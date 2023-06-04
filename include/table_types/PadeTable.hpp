/*
  LUT using [M/N] pade approximants with uniform sampling. Polynomial coefficients are calculated using
  Armadillo.

  Usage example using [4/3] approximants:
    PadeTable<4,3> look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - This class only works if TOUT and TIN can both be cast to double. 
    Armadillo Mat<T>'s `is_supported_elem_type<T>` will only let us do arithmetic
    with float or double (not even long double) and arma::field is useless.
    TODO disable constructor if TOUT cannot be cast to double

  - TODO add a way to build these tables with a pole on the left or right endpoints
  - table precomputes and stores any coefficients so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
  - Available template values are all M,N such that 0 < N <= M and M+N<=7
  - Template values where M < N are not supported
  - Requires both Armadillo and Boost version 1.71.0 or newer to generate
*/
#pragma once
#include "MetaTable.hpp"
#include "config.hpp" // FUNC_USE_BOOST, FUNC_USE_ARMADILLO
#include <stdexcept>
#include <cmath> //isinfite

#ifdef FUNC_USE_ARMADILLO
#include <armadillo>
#endif

namespace func {

static double constexpr fact[] = {1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0};

template <unsigned int M, unsigned int N, typename TIN, typename TOUT=TIN, GridTypes GT=GridTypes::UNIFORM>
class PadeTable final : public MetaTable<M+N+1,TIN,TOUT,GT>
{
  INHERIT_META(M+N+1,TIN,TOUT,GT);
public:
  // build the LUT from scratch or look in filename for an existing LUT
  PadeTable(const FunctionContainer<TIN,TOUT>& func_container, const LookupTableParameters<TIN>& par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<M+N+1,TIN,TOUT,GT>(jsonStats.empty() ? // use the default move constructor for MetaTable (probably not elided...)
      std::move(MetaTable<M+N+1,TIN,TOUT,GT>(func_container, par)) :
      std::move(MetaTable<M+N+1,TIN,TOUT,GT>(jsonStats)))
  {
#if !defined(FUNC_USE_BOOST) || !defined(FUNC_USE_ARMADILLO)
    /* This could theoretically be a compile time error; however, that will only stop us from registering this table (which is not useful!) */
    if(jsonStats.empty())
      throw std::invalid_argument("Error in func::PadeTable: Pade LUTs need both Armadillo and Boost to be generated");
#else
    if(!jsonStats.empty())
      return; // we already loaded the LUT from json in MetaTable

    using boost::math::differentiation::make_fvar;

    /* Base class default variables */
    m_name = grid_type_to_string<GT>() + "PadeTable<" + std::to_string(M) + "," + std::to_string(N) + ">";
    m_order = M+N+1;
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = sizeof(m_table[0]) * m_numTableEntries;

    auto boost_fun = func_container.template get_nth_func<M+N>();
    if(boost_fun == nullptr)
      throw std::invalid_argument(m_name + " needs the " + std::to_string(N+M) + "th derivative but this was not provided in FunctionContainer");

    m_grid.reset(new TIN[m_numTableEntries]);
    m_table.reset(new polynomial<TOUT,M+N+1>[m_numTableEntries]);
    FUNC_BUILDPAR
    for (unsigned int ii=0;ii<m_numTableEntries;++ii) {
      /* nonuniform grids are not supported for PadeTables */
      TIN x = m_minArg + ii*m_stepSize;
      TIN h = m_stepSize;
      m_grid[ii] = x;

      // build the matrix of taylor coefficients
      arma::Mat<double> T = arma::zeros<arma::Mat<double>>(M+N+1, N+1);
      const auto derivs = boost_fun(make_fvar<TIN,M+N>(x+0.5*h));
      //const auto derivs = boost_fun(make_fvar<TIN,M+N>(x));
      for(unsigned int i=0; i<M+N+1; i++)
        T(i,0) = static_cast<double>(derivs.derivative(i)*pow(h,i))/(fact[i]);

      /* copy the first column of T down and to the right  */
      for(unsigned int i=0; i<N+1; i++)
        T(arma::span(i,N+M), i) = T(arma::span(0,N+M-i), 0);

      /* find the coefficients of Q */
      arma::Mat<double> Q = arma::null(T.rows(M+1, M+N));

      /* Check if this Pade approximant behaves spurriously. That might be poles or poorly defined coefficients...
       * If the approximant is bad then lower the degree of Q until it is sufficiently not bad. */
      bool Q_bad = false;

      /* TODO Q might have more than 1 column even for analytic functions!!!
       * This is confounding because Pade approximants are supposed to be unique...????? */
      if(Q.n_elem != N+1){
        Q_bad = Q_bad || true;
        //std::cerr << "Pade approximant is bad and it doesn't make any mathematical sense...\n";
      }

      /* Check for roots of Q within the subinterval [0,h] by building a bracket.
       * Q(0)=1 is the positive endpoint so we need only check if Q is ever negative */
      for(unsigned int k=N; k>0; k--){
        /* scale Q such that its first entry is 1 */
        Q_bad = Q_bad || !(Q[0]); // first entry cannot be 0
        Q=Q/Q[0];

        // check Q at the right endpoint
        Q_bad = Q_bad || is_Q_negative(-static_cast<double>(h)/2,Q) || is_Q_negative(static_cast<double>(h)/2,Q);

        // Check Q for negativity at any of its extrema
        double desc = 0.0;
        switch(k){
          case 1:
            break;
          case 2:
            Q_bad = Q_bad || is_Q_negative(-Q[1]/(2.0*Q[2]),Q);
            break;
          case 3:
            desc = Q[2]*Q[2]-3*Q[1]*Q[3];
            Q_bad = Q_bad || (desc > 0.0 && (is_Q_negative(-Q[2]+sqrt(desc)/(3.0*Q[3]),Q) || is_Q_negative(-Q[2]-sqrt(desc)/(3.0*Q[3]),Q)));
            break;
        }

        // switch to the [M,k-1] pade approximant on this table interval
        if(Q_bad){
          Q = arma::zeros<arma::Mat<double>>(N+1);
          if(k == 1){
            Q[0] = 1.0; // use the degree M Taylor sum
          }else{
            Q.rows(0,k-1) = arma::null(T(arma::span(M+1, M+k-1), arma::span(0,k-1)));
          }
        }else break; // Q is chill
      }

      /* Compute the coefficients of P */
      arma::Col<double> P = T.rows(0,M)*Q;

      /* TODO we can save ourselves from subtracting dx by 0.5 in operator() if we get these 4 lines to work */
      //for(unsigned int k=0; k<N+1; k++)
      //  Q[k] = polynomial_diff(Q,-0.5,k)/boost::math::factorial<double>(k);
      //for(unsigned int k=0; k<M+1; k++)
      //  P[k] = polynomial_diff(P,-0.5,k)/boost::math::factorial<double>(k);

      // move these coefs into m_table
      for (unsigned int k=0; k<M+1; k++)
        m_table[ii].coefs[k] = static_cast<TOUT>(P[k]);

      for (unsigned int k=0; k<N; k++)
        m_table[ii].coefs[M+1+k] = static_cast<TOUT>(Q[k+1]); // ignore the first coef of Q b/c it's always 1.

    }
    // special case to make lut(tableMaxArg) work
    m_grid[m_numTableEntries-1] = m_tableMaxArg;
    m_table[m_numTableEntries-1].coefs[0] = func_container.standard_fun(m_tableMaxArg);
    for (unsigned int k=1; k<N+1; k++)
      m_table[m_numTableEntries-1].coefs[k] = 0;
#endif
  }

#if defined(FUNC_USE_BOOST) && defined(FUNC_USE_ARMADILLO)
  // TODO these should have a consistent signature
  /* check if x is within Q's domain and if Q(x) is negative */
  bool is_Q_negative(double x, arma::mat Q) {
    if(x < -m_stepSize/2 || m_stepSize/2 > x) return false;

    double sum = 0.0;
    for(int k=Q.n_elem; k>=0; k--)
      sum = Q[k-1] + x*sum;
    return (sum <= 0.0);
  };

  double polynomial_diff(arma::mat p, double x, unsigned s){
    double sum = 0.0;
    for(unsigned int k=p.n_elem; k>s; k--)
      sum = p[k-1]*permutation(k-1,s) + sum*x;
    return sum;
  }
#endif

  /* override operator() for rational functions
   * TODO can avoid subtracting dx by 0.5 if we shift P and Q after constructing them */
  TOUT operator()(TIN x) const final
  {
    unsigned int x0; TOUT dx;
    std::tie(x0,dx) = MetaTable<M+N+1,TIN,TOUT,GT>::template hash<GT>(x);
    dx -= 0.5;

    // general degree horners method, evaluated from the inside out.
    // TODO evaluate P and Q at the same time if possible
    TOUT P = dx*m_table[x0].coefs[M];
    for (int k=M-1; k>0; k--)
      P = dx*(m_table[x0].coefs[k] + P);
    P = P+m_table[x0].coefs[0];

    TOUT Q = dx*m_table[x0].coefs[M+N];
    for (int k=N-1; k>0; k--)
      Q = dx*(m_table[x0].coefs[M+k] + Q);
    Q = static_cast<TOUT>(1.0)+Q; // the constant term in Q will always be 1
    return P/Q;
  }
};

template <unsigned int M, unsigned int N, typename TIN, typename TOUT=TIN>
using UniformPadeTable = PadeTable<M,N,TIN,TOUT,GridTypes::UNIFORM>;
} // namespace func

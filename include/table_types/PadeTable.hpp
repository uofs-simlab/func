/*
 *  TODO FIX THE OPERATOR()
 *  - use armadillo to solve the uppertriangular system
 *  - 
 *
 *
 *
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
      return; // all our work is already done

    using boost::math::differentiation::make_fvar;

    /* Base class default variables */
    m_name = grid_type_to_string<GT>() + "PadeTable<" + std::to_string(N) + ">";
    m_order = M+N+1;
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = static_cast<unsigned>(sizeof(m_table[0]) * (m_numTableEntries));

    auto boost_fun = func_container.template get_nth_func<M+N>();
    if(boost_fun == nullptr)
      throw std::invalid_argument(m_name + " needs the " + std::to_string(N+M) + "th derivative but this was not provided in FunctionContainer");

    /* Allocate and set table */
    m_grid.reset(new TIN[m_numTableEntries]);
    m_table.reset(new polynomial<TOUT,M+N+1>[m_numTableEntries]);
    FUNC_BUILDPAR
    for (unsigned int ii=0;ii<m_numTableEntries;++ii) {
      // nonuniform grids are not supported for PadeTables
      TIN x = m_minArg + ii*m_stepSize;
      // grid points
      m_grid[ii] = x;

      // build the matrix of taylor coefficients
      arma::Mat<double> T = arma::zeros<arma::Mat<double>>(M+N+1, N+1);
      const auto derivs = boost_fun(make_fvar<TIN,M+N>(x));
      for(unsigned int i=0; i<M+N+1; i++)
        T(i,0) = static_cast<double>(derivs.derivative(i))/(fact[i]);

      // copy the first column of T down and to the right
      for(unsigned int i=0; i<N+1; i++)
        T(arma::span(i,N+M), i) = T(arma::span(0,N+M-i), 0);

      // find the coefficients of Q.
      arma::Mat<double> Q = arma::null(T.rows(M+1, M+N));
      bool Q_has_root = false; // TODO using brent_find_minima to check for a minimum will result in less source code

      /* TODO This can happen! Which is particularly confounding because Pade approximants are supposed to be unique...?? */
      if(Q.n_elem != N+1){
        Q = Q.col(0);
        Q_has_root = true;
        //throw std::runtime_error("Error in func::PadeTable:" + m_name + " is too poorly conditioned because the matrix Q has nullspace with dimension > 1");
      }

      // scale Q such that its first entry equals 1.
      Q_has_root = !(Q[0]); // first entry cannot be 0
      Q=Q/Q[0];

      //std::cout << Q << "\n";

      // find the coefficients of P
      arma::Col<double> P = T.rows(0,M)*Q;

      /* Check if the Q has any roots within the subinterval [0,m_stepSize]
       * by building a bracket (Q(0)=1 is the positive endpoint so we just need a negative endpoint).
       * If any roots exist, then lower the degree of Q. */
      for(unsigned int k=N; k>0; k--){
        // check Q at the subinterval endpoints
        Q_has_root = Q_has_root || Q_is_negative(static_cast<double>(-m_stepSize/2.0),Q,ii) || Q_is_negative(static_cast<double>(m_stepSize/2.0),Q,ii);

        // Check Q for negativity at any of its vertexes
        double desc = 0.0;
        if(!Q_has_root)
          switch(k){
            case 1:
              break;
            case 2:
              Q_has_root = Q_is_negative(-Q[1]/(2.0*Q[2]),Q,ii);
              break;
            case 3:
              desc = Q[2]*Q[2]-3*Q[1]*Q[3];
              Q_has_root = desc > 0.0 &&
                (Q_is_negative(-Q[2]+sqrt(desc)/(3*Q[3]),Q,ii) || Q_is_negative(-Q[2]+sqrt(desc)/(3*Q[3]),Q,ii));
              break;
          }

        // switch to using the [M,k-1] pade approximant on this table interval
        if(Q_has_root){
          if(k == 1){
            // just use a Taylor series
            Q = arma::zeros<arma::Mat<double>>(N+1);
            Q[0] = 1.0;
            P = T(arma::span(0,M),0);
          }else{
            Q[k] = 0.0;
            Q.rows(0,k-1) = arma::null(T(arma::span(M+1, M+k-1), arma::span(0,k-1)));
            Q_has_root = !(Q[0]);
            Q = Q/Q[0];
            P = T.rows(0,M)*Q;
          }
        }else // Q is free to go if it has no roots here
          break;
      }

      /* TODO is this exception ever thrown? */
      for(unsigned int i=1; i<M+N+1; i++)
        if(!std::isfinite(Q[i])) // check for NaN or +/-inf
          throw std::runtime_error(m_name + " is too poorly conditioned: coef " + std::to_string(i) + " of Q is " + std::to_string(Q[i]));

      // move these coefs into m_table
      for (unsigned int k=0; k<M+1; k++)
        m_table[ii].coefs[k] = static_cast<TOUT>(P[k]);

      for (unsigned int k=0; k<N; k++)
        m_table[ii].coefs[M+1+k] = static_cast<TOUT>(Q[k+1]); // ignore the first coef of Q b/c it's always 1.
    }
    // does not need special case grid entry!
#endif
  }

#if defined(FUNC_USE_BOOST) && defined(FUNC_USE_ARMADILLO)
    bool Q_is_negative(double x, arma::mat Q, unsigned int ii) {
      // Check if x is within this subinterval's range
      if(((ii == 0 && x < 0.0) || (ii == m_numTableEntries - 1 && x > 0.0)))
        return false;

      // compute Q(x) using horners method evaluating from the inside out
      double sum = 0;
      for(int k=N; k>0; k--)
        sum = x*(Q[k] + sum);
      sum += 1;
      // was Q(x) negative?
      return sum <= 0.0;
    };
#endif

  // override operator() from MetaTable so we work with rational functions instead
  TOUT operator()(TIN x) const final
  {
    // nondimensionalized x position
    TOUT dx  = (x-m_minArg);
    TOUT x1r = dx/m_stepSize+0.5;
    // index of previous table entry
    unsigned x1 = ((unsigned) x1r);
    dx -= x1*m_stepSize;

    // general degree horners method, evaluated from the inside out.
    // TODO evaluate P and Q at the same time if possible
    TOUT P = dx*m_table[x1].coefs[M];
    for (int k=M-1; k>0; k--)
      P = dx*(m_table[x1].coefs[k] + P);
    P = P+m_table[x1].coefs[0];

    TOUT Q = dx*m_table[x1].coefs[M+N];
    for (int k=N-1; k>0; k--)
      Q = dx*(m_table[x1].coefs[M+k] + Q);
    Q = 1+Q; // the constant term in Q will always be 1
    return P/Q;
  }
};

template <unsigned int M, unsigned int N, typename TIN, typename TOUT=TIN>
using UniformPadeTable = PadeTable<M,N,TIN,TOUT,GridTypes::UNIFORM>;
} // namespace func

/* Implementation of a UniformPrecomputed Lookup table with linear interpolation */
#include "UniformPadeTable.hpp"
#include <armadillo>
#include <iostream>
#include <stdexcept>
#include <cmath>

// The registration looks terrible so it's at the bottom of the class
using boost::math::differentiation::make_fvar;

static double const fact[] = {1,1,2,6,24,120,720,5040};

template <unsigned int M, unsigned int N>
UniformPadeTable<M,N>::UniformPadeTable(FunctionContainer *func_container, UniformLookupTableParameters par) : 
  UniformLookupTable(func_container, par)
{
  /* Base class default variables */
  m_name = "UniformPadeTable<" + std::to_string(M) + "," + std::to_string(N) + ">";
  m_order = M+N+1;  
  m_numTableEntries = m_numIntervals+1;
  m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;

  __IS_NULL(func_container->get_nth_func<M+N>()); // the least descriptive exception
  mp_boost_func = func_container->get_nth_func<M+N>();

  /* Allocate and set table */
  m_table.reset(new polynomial<M+N+1, M+N<4? 32:64>[m_numTableEntries]);
  for (int ii=0;ii<m_numIntervals;++ii) {
    const double x = m_minArg + ii*m_stepSize;
    // grid points
    m_grid[ii] = x;
    
    // build the matrix of taylor coefficients
    arma::mat T = arma::zeros(M+N+1, N+1);
    const auto derivs = (mp_boost_func)(make_fvar<double,M+N>(x));
    for(unsigned int i=0; i<M+N+1; i++)
      T(i,0) = derivs.derivative(i)/(fact[i]);

    // copy a shifted column 1 of T into the rest of the matrix
    for(unsigned int i=0; i<N+1; i++)
      T(arma::span(i,N+M), i) = T(arma::span(0,N+M-i), 0);

    // find the coefficients of Q.
    arma::mat Q = arma::null(T.rows(M+1, M+N));

    // scale Q such that its first entry equals 1.
    Q=Q/Q[0];
    if(!std::isfinite(Q[1]))
      std::cerr << "Q = " << Q << std::endl << std::endl;

    // find the coefficients of P
    arma::vec P = T.rows(0,M)*Q;

    // check if the roots of Q lie in this interval
    // this is a bit of a mess of math/logic. We just need to know if any roots
    // exist within \pm(m_stepSize/2) of (m_minArg + ii*m_stepSize)
    // with a special case at the table's endpoints
    auto Q_is_negative = [this, &Q, &ii](double x) -> bool {
      // Tell us if this point is within this subinterval's range
      if(((ii == 0 && x < 0.0) || (ii == m_numIntervals - 1 && x > 0.0)))
        return false;

      // compute Q(x)
      double sum = x*Q[N];
      for (int k=N-1; k>0; k--)
        sum = x*(Q[k] + sum);
      sum += 1;
      // Tell us if this point is negative
      return sum < 0.0;
    };

    //for(unsigned int k=N; k>0; k--){
    //  // check if Q is negative at the subinterval endpoints
    //  bool Q_has_root = Q_is_negative(-m_stepSize/2.0) || Q_is_negative(m_stepSize/2.0);
    //  // Check if Q is negative at any of its vertexes
    //  if(!Q_has_root)
    //    switch(k){
    //      case 1:
    //        break;
    //      case 2:
    //        Q_has_root  = Q_is_negative(-Q[1]/(2.0*Q[2]));
    //        break;
    //      case 3:
    //        double desc = Q[2]*Q[2]-3*Q[1]*Q[3];
    //        Q_has_root  = desc > 0.0 && (Q_is_negative(-Q[2]+sqrt(desc)/(3*Q[3])) || Q_is_negative(-Q[2]+sqrt(desc)/(3*Q[3])));
    //        break;
    //    }
    //  // switch to the [M,N-1] pade approximant on this table interval
    //  if(Q_has_root){
    //    if(k == 1){
    //      Q = arma::zeros(N+1);
    //      Q[0] = 1.0;
    //      P = T(arma::span(0,M),0);
    //    }else{
    //      Q[k] = 0.0;
    //      Q.rows(0,k-1) = arma::null(T(arma::span(M+1, M+k-1), arma::span(0,k-1));
    //      Q = Q/Q[0];
    //      P = T(arma::span(0,M-2),arma::span(0,k-1))*Q;
    //    }
    //  }else // Q is free to go if it has no roots here
    //    break;
    //}

    bool Q_has_root = false;
    double desc = 0.0;
    switch(N){
      case 3 :
        Q_has_root = Q_is_negative(-m_stepSize/2.0) || Q_is_negative(m_stepSize/2.0);
        desc = Q[2]*Q[2]-3*Q[1]*Q[3];
        Q_has_root = Q_has_root || (desc > 0.0 && 
            (Q_is_negative(-Q[2]+sqrt(desc)/(3*Q[3])) || Q_is_negative(-Q[2]+sqrt(desc)/(3*Q[3]))));
        if(Q_has_root){
          Q[3] = 0.0;
          Q.rows(0,2) = arma::null(T(arma::span(M+1, M+2), arma::span(0,2)));
          Q = Q/Q[0];
          P = T.rows(0,M)*Q;
        }else
          break;
      case 2 :
        Q_has_root = Q_is_negative(-m_stepSize/2.0) || Q_is_negative(m_stepSize/2.0);
        Q_has_root = Q_has_root || Q_is_negative(-Q[1]/(2.0*Q[2]));
        if(Q_has_root){
          Q[2] = 0.0;
          Q.rows(0,1) = arma::null(T(arma::span(M+1, M+1), arma::span(0,1)));
          Q = Q/Q[0];
          P = T.rows(0,M)*Q;
        }else
          break;
      case 1 :
        Q_has_root = Q_is_negative(-m_stepSize/2.0) || Q_is_negative(m_stepSize/2.0);
        if(Q_has_root){
          Q = arma::zeros(N+1);
          Q[0] = 1.0;
          P = T(arma::span(0,M),0);
        }
    }

    // move these coefs into m_table
    for (unsigned int k=0; k<M+1; k++)
      m_table[ii].coefs[k] = P[k];

    for (unsigned int k=0; k<N; k++)
      m_table[ii].coefs[M+1+k] = Q[k+1]; // ignore the first coef of Q b/c it's always 1.
  }
}

template <unsigned int M, unsigned int N>
double UniformPadeTable<M,N>::operator()(double x)
{ 
  // nondimensionalized x position
  double  dx  = (x-m_minArg);
  double  x1r = dx/m_stepSize+0.5;
  // index of previous table entry
  unsigned x1 = ((unsigned) x1r);
  dx -= x1*m_stepSize;
  
  // general degree horners method, evaluated from the inside out.
  double P = dx*m_table[x1].coefs[M];
  for (int k=M-1; k>0; k--)
    P = dx*(m_table[x1].coefs[k] + P);
  P = P+m_table[x1].coefs[0];

  double Q = dx*m_table[x1].coefs[M+N];
  for (int k=N-1; k>0; k--)
    Q = dx*(m_table[x1].coefs[M+k] + Q);
  Q = 1+Q;  // the constant term in Q will always be 1
  return P/Q;
}

// Template substitution happens way after the preprocessor does it's work so
// we'll register all the available template values this way
template<> REGISTER_ULUT_IMPL(UniformPadeTable<1,1>);
template<> REGISTER_ULUT_IMPL(UniformPadeTable<2,1>);
template<> REGISTER_ULUT_IMPL(UniformPadeTable<3,1>);
template<> REGISTER_ULUT_IMPL(UniformPadeTable<4,1>);
template<> REGISTER_ULUT_IMPL(UniformPadeTable<5,1>);
template<> REGISTER_ULUT_IMPL(UniformPadeTable<6,1>);

template<> REGISTER_ULUT_IMPL(UniformPadeTable<2,2>);
template<> REGISTER_ULUT_IMPL(UniformPadeTable<3,2>);
template<> REGISTER_ULUT_IMPL(UniformPadeTable<4,2>);
template<> REGISTER_ULUT_IMPL(UniformPadeTable<5,2>);

template<> REGISTER_ULUT_IMPL(UniformPadeTable<3,3>);
template<> REGISTER_ULUT_IMPL(UniformPadeTable<4,3>);

// declaration of the available values for the template values M,N
template class UniformPadeTable<1,1>;
template class UniformPadeTable<2,1>;
template class UniformPadeTable<3,1>;
template class UniformPadeTable<4,1>;
template class UniformPadeTable<5,1>;
template class UniformPadeTable<6,1>;

template class UniformPadeTable<2,2>;
template class UniformPadeTable<3,2>;
template class UniformPadeTable<4,2>;
template class UniformPadeTable<5,2>;

template class UniformPadeTable<3,3>;
template class UniformPadeTable<4,3>;

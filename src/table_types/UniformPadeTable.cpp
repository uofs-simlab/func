/* Implementation of a UniformPrecomputed Lookup table with linear interpolation */
#include "UniformPadeTable.hpp"
#include <armadillo>
#include <iostream>

// The registration looks terrible so it's at the bottom of the class

static double const fact[] = {1,1,2,6,24,120,720,5040};
double (EvaluationFunctor<double,double>::*derivs[7])(double);

template <unsigned int M, unsigned int N>
UniformPadeTable<M,N>::UniformPadeTable(EvaluationFunctor<double,double> *func, UniformLookupTableParameters par) : 
  UniformLookupTable(func, par)
{
  /* Base class default variables */
  m_name = "UniformPadeTable<" + std::to_string(M) + "," + std::to_string(N) + ">";
  m_order = M+N+1;  
  m_numTableEntries = (M+N+1)*(m_numIntervals+1);
  m_dataSize = (unsigned) sizeof(double) * m_numTableEntries;

  // assign the first 7 derivatives to the derivs array for easy enumeration
  derivs[0]=&EvaluationFunctor<double,double>::deriv;
  derivs[1]=&EvaluationFunctor<double,double>::deriv2;
  derivs[2]=&EvaluationFunctor<double,double>::deriv3;
  derivs[3]=&EvaluationFunctor<double,double>::deriv4;
  derivs[4]=&EvaluationFunctor<double,double>::deriv5;
  derivs[5]=&EvaluationFunctor<double,double>::deriv6;
  derivs[6]=&EvaluationFunctor<double,double>::deriv7;

  /* Allocate and set table */
  m_table.reset(new double[m_numTableEntries]);
  for (int ii=0;ii<m_numIntervals;++ii) {
    const double x = m_minArg + ii*m_stepSize;
    // grid points
    m_grid[ii] = x;
    
    // build the matrix of taylor coefficients
    arma::mat T = arma::zeros(M+N+1, N+1);
    T(0,0) = (*func)(x);
    for(unsigned int i=1; i<M+N+1; i++)
      T(i,0) = (func->*(derivs[i-1]))(x)/(fact[i]);
    
    // copy a shifted column 1 of T into the rest of the matrix
    for(unsigned int i=0; i<N+1; i++)
      T(arma::span(i,N+M), i) = T(arma::span(0,N+M-i), 0);

    // find the coefficients of Q.
    arma::mat Q = arma::null(T.rows(M+1, M+N)); 
    if(Q.n_cols > 1){ // Arma could throw an exception if Q isn't a vector but this is more descriptive
      throw "Pade table of order [" + std::to_string(M) + "/" + std::to_string(N) + "] does not exist.";
      return;
    }

    // scale Q such that its first entry equals 1.
    Q=Q/Q[0]; 

    // find the coefficients of P
    arma::vec P = T.rows(0,M)*Q;

    // move these coefs into m_table
    for (unsigned int k=0; k<M+1; k++)
      m_table[(M+N+1)*ii+k] = P[k];

    for (unsigned int k=0; k<N; k++)
      m_table[(M+N+1)*ii+M+1+k] = Q[k+1]; // ignore the first coef of Q b/c it's always 1.
  }
}

template <unsigned int M, unsigned int N>
double UniformPadeTable<M,N>::operator()(double x)
{ 
  // nondimensionalized x position
  double  dx  = (x-m_minArg);
  double  x1r = dx/m_stepSize+0.5;
  // index of previous table entry
  unsigned x1 = (M+N+1)*((unsigned) x1r);
  dx -= x1*m_stepSize/(M+N+1);
  
  // general degree horners method, evaluated from the inside out.
  double P = dx*m_table[x1+M];
  for (int k=M-1; k>0; k--)
    P = dx*(m_table[x1+k] + P);
  P = P+m_table[x1];

  double Q = dx*m_table[x1+M+N];
  for (int k=N-1; k>0; k--)
    Q = dx*(m_table[x1+M+k] + Q);
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
/* Implementation of a UniformPrecomputed Lookup table with linear interpolation */
#include "UniformArmadilloPrecomputedInterpolationTable.hpp"
#include <armadillo>

#define SOLVE_OPTS arma::solve_opts::refine

// Template substitution happens way after the preprocessor does it's work so
// we'll register all the available template values this way
template<>REGISTER_ULUT_IMPL(UniformArmadilloPrecomputedInterpolationTable<4>);
template<>REGISTER_ULUT_IMPL(UniformArmadilloPrecomputedInterpolationTable<5>);
template<>REGISTER_ULUT_IMPL(UniformArmadilloPrecomputedInterpolationTable<6>);
template<>REGISTER_ULUT_IMPL(UniformArmadilloPrecomputedInterpolationTable<7>);

template <unsigned int N>
UniformArmadilloPrecomputedInterpolationTable<N>::UniformArmadilloPrecomputedInterpolationTable(EvaluationFunctor<double,double> *func, UniformLookupTableParameters par) : 
  UniformLookupTable(func, par)
{
  /* Base class default variables */
  m_name = "UniformArmadilloPrecomputedInterpolationTable<" + std::to_string(N) + ">";
  m_order = N+1;  // take N as the degree of the polynomial interpolant which is of order N+1
  m_numTableEntries = m_numIntervals+1;
  m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;
   
  /* build the vandermonde system for finding the interpolating polynomial's coefficients */
  arma::mat Van = arma::ones(N+1, N+1);
  Van.col(1) = arma::linspace(0,1,N+1);
  for(unsigned int i=2; i<N+1; i++)
    Van.col(i) = Van.col(i-1) % Van.col(1); // the % does elementwise multiplication
  
  // LU factor the matrix we just built
  arma::mat L, U, P;
  arma::lu(L,U,P,Van);

  /* Allocate and set table */
  m_table.reset(new polynomial<N+1,64>[m_numTableEntries]);
  for (int ii=0;ii<m_numIntervals;++ii) {
    const double x = m_minArg + ii*m_stepSize;
    // grid points
    m_grid[ii] = x;
    
    // build the vector of coefficients from uniformly spaced function values
    arma::vec y = arma::linspace(x,x+m_stepSize,N+1);
    for (unsigned int k=0; k<N+1; k++)
      y[k] = (*mp_func)(y[k]);
    
    // make y the coefficients of the polynomial interpolant
    y = solve(trimatu(U), solve(trimatl(L), P*y));
    //y = arma::solve(Van, y, SOLVE_OPTS);
    
    // move this back into the m_table array
    for (unsigned int k=0; k<N+1; k++)
      m_table[ii].coefs[k] = y[k];
  }
}

template <unsigned int N>
double UniformArmadilloPrecomputedInterpolationTable<N>::operator()(double x)
{
  // nondimensionalized x position, scaled by step size
  double   dx = m_stepSize_inv*(x-this->m_minArg);
  // index of previous table entry
  unsigned x0  = (unsigned) dx;
  // value of table entries around x position
  dx -= x0;
  
  // general degree horners method, evaluated from the inside out.
  double sum = dx*m_table[x0].coefs[N];
  for (int k=N-1; k>0; k--)
    sum = dx*(m_table[x0].coefs[k] + sum);
  return m_table[x0].coefs[0]+sum;
}

// declaration of the available values for the template N
template class UniformArmadilloPrecomputedInterpolationTable<4>;
template class UniformArmadilloPrecomputedInterpolationTable<5>;
template class UniformArmadilloPrecomputedInterpolationTable<6>;
template class UniformArmadilloPrecomputedInterpolationTable<7>;

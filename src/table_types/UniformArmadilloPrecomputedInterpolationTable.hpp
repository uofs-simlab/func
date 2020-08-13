/*
  4th to 7th degree polynomial interpolation LUT with uniform sampling (precomputed coefficients using an
  Armadillo matrix)

  Usage example for a 4th degree interpolant:
    UniformArmadilloPrecomputedInterpolationTable<4> look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
  - the template implementation is only for N=4,5,6,7 (ie, available polynomial interpolation is 
  of degrees 4 up to degree 7)
*/
#pragma once
#include "UniformLookupTable.hpp"
#include <cmath> // ceil()

#define ARMA_USE_CXX11
#include <armadillo>
#define SOLVE_OPTS arma::solve_opts::none

template <typename IN_TYPE, typename OUT_TYPE, unsigned int N>
class UniformArmadilloPrecomputedInterpolationTable final : public UniformLookupTable<IN_TYPE,OUT_TYPE>
{
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE);

  REGISTER_LUT(UniformArmadilloPrecomputedInterpolationTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,N+1>[]> m_table;

public:
  UniformArmadilloPrecomputedInterpolationTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container,
      UniformLookupTableParameters<IN_TYPE> par) :
    UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, par)
  {
    /* Base class default variables */
    m_name = "UniformArmadilloPrecomputedInterpolationTable<" + std::to_string(N) + ">";
    m_order = N+1;  // take N as the degree of the polynomial interpolant which is of order N+1
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = (unsigned) sizeof(m_table[0]) * (m_numTableEntries);
     
    /* build the vandermonde system for finding the interpolating polynomial's coefficients */
    arma::mat Van = arma::ones(N+1, N+1);
    Van.col(1) = arma::linspace(0,1,N+1);
    for(unsigned int i=2; i<N+1; i++)
      Van.col(i) = Van.col(i-1) % Van.col(1); // the % does elementwise multiplication
    
    // LU factor the matrix we just built
    arma::mat L, U, P;
    arma::lu(L,U,P,Van);

    /* Allocate and set table */
    m_table.reset(new polynomial<OUT_TYPE,N+1>[m_numTableEntries]);
    for (int ii=0;ii<m_numIntervals;++ii) {
      const IN_TYPE x = m_minArg + ii*m_stepSize;
      // grid points
      m_grid[ii] = x;
      
      // build the vector of coefficients from uniformly spaced function values
      arma::vec y = arma::linspace(x,x+m_stepSize,N+1);
      for (unsigned int k=0; k<N+1; k++)
        y[k] = mp_func(y[k]);
      
      // make y the coefficients of the polynomial interpolant
      y = solve(trimatu(U), solve(trimatl(L), P*y));
      //y = arma::solve(Van, y, SOLVE_OPTS);
      
      // move this back into the m_table array
      for (unsigned int k=0; k<N+1; k++)
        m_table[ii].coefs[k] = y[k];
    }
  }

  OUT_TYPE operator()(IN_TYPE x) override
  {
    // nondimensionalized x position, scaled by step size
    OUT_TYPE dx = (OUT_TYPE) m_stepSize_inv*(x-m_minArg);
    // index of previous table entry
    unsigned x0  = (unsigned) dx;
    // value of table entries around x position
    dx -= x0;
    
    // general degree horners method, evaluated from the inside out.
    OUT_TYPE sum = dx*m_table[x0].coefs[N];
    for (int k=N-1; k>0; k--)
      sum = dx*(m_table[x0].coefs[k] + sum);
    return m_table[x0].coefs[0]+sum;
  }
};

// Template substitution happens way after the preprocessor does it's work so
// we'll register all the available template values this way
REGISTER_TEMPLATED_DOUBLE_AND_FLOAT_LUT_IMPLS(UniformArmadilloPrecomputedInterpolationTable,4);
REGISTER_TEMPLATED_DOUBLE_AND_FLOAT_LUT_IMPLS(UniformArmadilloPrecomputedInterpolationTable,5);
REGISTER_TEMPLATED_DOUBLE_AND_FLOAT_LUT_IMPLS(UniformArmadilloPrecomputedInterpolationTable,6);
REGISTER_TEMPLATED_DOUBLE_AND_FLOAT_LUT_IMPLS(UniformArmadilloPrecomputedInterpolationTable,7);

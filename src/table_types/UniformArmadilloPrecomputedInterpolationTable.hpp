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
#include "UniformLookupTable.hpp"
#define ARMA_USE_CXX11
#include <armadillo>
#define SOLVE_OPTS arma::solve_opts::none

template <typename IN_TYPE, typename OUT_TYPE, unsigned int N>
class UniformArmadilloPrecomputedInterpolationTable final : public UniformLookupTable<IN_TYPE,OUT_TYPE>
{
  // Template substitution happens way after the preprocessor does it's work so
  // we'll register all the available template values this way
  REGISTER_LUT(UniformArmadilloPrecomputedInterpolationTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,N+1,64>[]> m_table;

public:
  UniformArmadilloPrecomputedInterpolationTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container,
      UniformLookupTableParameters<IN_TYPE> par) :
    UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, par)
  {
    /* Base class default variables */
    this->m_name = "UniformArmadilloPrecomputedInterpolationTable<" + std::to_string(N) + ">";
    this->m_order = N+1;  // take N as the degree of the polynomial interpolant which is of order N+1
    this->m_numTableEntries = this->m_numIntervals+1;
    this->m_dataSize = (unsigned) sizeof(m_table[0]) * this->m_numTableEntries;
     
    /* build the vandermonde system for finding the interpolating polynomial's coefficients */
    arma::mat Van = arma::ones(N+1, N+1);
    Van.col(1) = arma::linspace(0,1,N+1);
    for(unsigned int i=2; i<N+1; i++)
      Van.col(i) = Van.col(i-1) % Van.col(1); // the % does elementwise multiplication
    
    // LU factor the matrix we just built
    arma::mat L, U, P;
    arma::lu(L,U,P,Van);

    /* Allocate and set table */
    m_table.reset(new polynomial<OUT_TYPE,N+1,64>[this->m_numTableEntries]);
    for (int ii=0;ii<this->m_numIntervals;++ii) {
      const IN_TYPE x = this->m_minArg + ii*this->m_stepSize;
      // grid points
      this->m_grid[ii] = x;
      
      // build the vector of coefficients from uniformly spaced function values
      arma::vec y = arma::linspace(x,x+this->m_stepSize,N+1);
      for (unsigned int k=0; k<N+1; k++)
        y[k] = this->mp_func(y[k]);
      
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
    IN_TYPE dx = this->m_stepSize_inv*(x-this->m_minArg);
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

register_double_and_float_types(UniformArmadilloPrecomputedInterpolationTable,4)
register_double_and_float_types(UniformArmadilloPrecomputedInterpolationTable,5)
register_double_and_float_types(UniformArmadilloPrecomputedInterpolationTable,6)
register_double_and_float_types(UniformArmadilloPrecomputedInterpolationTable,7)

/* Implementation of a UniformPrecomputed Lookup table with linear interpolation */
#include "UniformQuarticPrecomputedInterpolationTable.hpp"
#include <armadillo>
#include <iostream>

#define IMPL_NAME UniformQuarticPrecomputedInterpolationTable
REGISTER_ULUT_IMPL(IMPL_NAME);

arma::mat calc_vander(){
  arma::mat A = arma::ones(5, 5);
  A.col(1) = arma::linspace(0,1,5);
  for(unsigned int i=2; i<5; i++)
    A.col(i) = A.col(i-1) % A.col(1); // % does elementwise multiplication
  return A;
}

static arma::mat const vander = calc_vander();

UniformQuarticPrecomputedInterpolationTable::UniformQuarticPrecomputedInterpolationTable(EvaluationFunctor<double,double> *func, UniformLookupTableParameters par) : 
  UniformLookupTable(func, par)
{
  /* Base class default variables */
  m_name = STR(IMPL_NAME);
  m_order = 5;
  m_numTableEntries = 5*(m_numIntervals+1);
  m_dataSize = (unsigned) sizeof(double) * m_numTableEntries;
  
  /* Allocate and set table */
  m_table.reset(new double[m_numTableEntries]);
  for (int ii=0;ii<m_numIntervals;++ii) {
    const double x = m_minArg + ii*m_stepSize;
    // grid points
    m_grid[ii] = x;
    // polynomial coefficients
    arma::vec y=arma::linspace(x,x+m_stepSize,5);
    y[0]=(*mp_func)(y[0]);
    y[1]=(*mp_func)(y[1]);
    y[2]=(*mp_func)(y[2]);
    y[3]=(*mp_func)(y[3]);
    y[4]=(*mp_func)(y[4]);
    
    y=solve(vander,y);

    m_table[4*ii]   = y[0];
    m_table[4*ii+1] = y[1];
    m_table[4*ii+2] = y[2];
    m_table[4*ii+3] = y[3];
    m_table[4*ii+4] = y[4];
  }
}

double UniformQuarticPrecomputedInterpolationTable::operator()(double x)
{
  // nondimensionalized x position, scaled by step size
  double   dx = m_stepSize_inv*(x-this->m_minArg);
  // index of previous table entry
  unsigned x0  = (unsigned) dx;
  // value of table entries around x position
  dx -= x0;
  x0 *= 5;
  // quartic horners method
  return m_table[x0]+dx*(m_table[x0+1]+dx*(m_table[x0+2]+dx*(m_table[x0+3]+dx*m_table[x0+4])));
}

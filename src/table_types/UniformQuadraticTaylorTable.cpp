/* Implementation of a Uniform Lookup table with linear interpolation */
#include "UniformQuadraticTaylorTable.hpp"

#define IMPL_NAME UniformQuadraticTaylorTable
REGISTER_ULUT_IMPL(IMPL_NAME);

UniformQuadraticTaylorTable::UniformQuadraticTaylorTable(EvaluationFunctor<double,double> *func, UniformLookupTableParameters par) : UniformLookupTable(func, par)
{
  /* Base class default variables */
  m_name = STR(IMPL_NAME);
  m_order = 3;
  m_numTableEntries = 3*m_numIntervals;
  m_dataSize = (unsigned) sizeof(double) * m_numTableEntries;

  /* Allocate and set table */
  m_table.reset(new double[m_numTableEntries]);
  for (int ii=0;ii<m_numIntervals;++ii) {
    double x = (m_minArg + ii*m_stepSize);
    m_grid[ii] = x;
    m_table[3*ii] = (*func)(x);
    m_table[3*ii+1] = func->deriv(x);
    m_table[3*ii+2] = func->deriv2(x)/2;
  }
}

double UniformQuadraticTaylorTable::operator()(double x)
{
  // nondimensionalized x position
  double  dx  = (x-m_minArg);
  double  x1r = dx/m_stepSize+0.5;
  // index of previous table entry
  unsigned x1 = 3*((unsigned) x1r);
  dx -= x1*m_stepSize/3;
  // linear Taylor series from grid point below x
  return m_table[x1]+dx*(m_table[x1+1]+dx*m_table[x1+2]);
}

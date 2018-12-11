/* Implementation of a Uniform Lookup table with linear interpolation */
#include "UniformLinearTaylorTable.hpp"

#define IMPL_NAME UniformLinearTaylorTable
REGISTER_ULUT_IMPL(IMPL_NAME);

UniformLinearTaylorTable::UniformLinearTaylorTable(EvaluationFunctor<double,double> *func, UniformLookupTableParameters par) : UniformLookupTable(func, par)
{
  /* Base class variables */
  m_name = STR(IMPL_NAME);
  m_order = 2;
  m_numTableEntries = 2*m_numIntervals;
  m_dataSize = (unsigned) sizeof(double) * m_numTableEntries;

  /* Allocate and set table */
  m_table.reset(new double[m_numTableEntries]);
  for (int ii=0; ii<m_numIntervals; ++ii) {
    const double x = (m_minArg + ii*m_stepSize);
    m_grid[ii]      = x;
    m_table[2*ii]   = (*mp_func)(x);
    m_table[2*ii+1] = mp_func->deriv(x);
  }
}

double UniformLinearTaylorTable::operator()(double x)
{
  // nondimensionalized x position
  double  dx  = (x-m_minArg);
  double  x1r = dx/m_stepSize+0.5;
  // index of previous table entry
  unsigned x1 = 2*((unsigned) x1r);
  dx -= 0.5*x1*m_stepSize;
  // linear Taylor series from grid point below x
  return m_table[x1]+m_table[x1+1]*dx;
}

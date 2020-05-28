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
  m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;

  /* Allocate and set table */
  m_table.reset(new polynomial<2,16>[m_numTableEntries]);
  for (int ii=0; ii<m_numIntervals; ++ii) {
    const double x = (m_minArg + ii*m_stepSize);
    m_grid[ii]      = x;
    m_table[ii].coefs[0]   = (*mp_func)(x);
    m_table[ii].coefs[1] = mp_func->deriv(x);
  }
}

double UniformLinearTaylorTable::operator()(double x)
{
  // nondimensionalized x position
  double  dx  = (x-m_minArg);
  double  x1r = dx/m_stepSize+0.5;
  // index of previous table entry
  unsigned x1 = (unsigned) x1r;
  dx -= x1*m_stepSize;
  // linear Taylor series from grid point below x
  return m_table[x1].coefs[0]+m_table[x1].coefs[1]*dx;
}

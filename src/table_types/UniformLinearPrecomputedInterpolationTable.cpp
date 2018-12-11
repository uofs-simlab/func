/* Implementation of a UniformPrecomputed Lookup table with linear interpolation */
#include "UniformLinearPrecomputedInterpolationTable.hpp"

#define IMPL_NAME UniformLinearPrecomputedInterpolationTable
REGISTER_ULUT_IMPL(IMPL_NAME);

UniformLinearPrecomputedInterpolationTable::UniformLinearPrecomputedInterpolationTable(EvaluationFunctor<double,double> *func, UniformLookupTableParameters par) : UniformLookupTable(func, par)
{

  /* Base class default variables */
  m_name = STR(IMPL_NAME);
  m_order = 2;
  m_numTableEntries = 2*m_numIntervals+2;
  m_dataSize = (unsigned) sizeof(double) * m_numTableEntries;

  /* Allocate and set table */
  m_table.reset(new double[m_numTableEntries]);
  for (int ii=0;ii<m_numIntervals;++ii) {
    double x = m_minArg + ii*m_stepSize;
    m_grid[ii] = x;
    m_table[2*ii]   = (*mp_func)(x);
    x = m_minArg + (ii+1)*m_stepSize;
    m_table[2*ii+1] = (*mp_func)(x) - m_table[2*ii];
  }
}

double UniformLinearPrecomputedInterpolationTable::operator()(double x)
{
  // nondimensionalized x position, scaled by step size
  double   dx = m_stepSize_inv*(x-m_minArg);
  // index of previous table entry
  unsigned x0  = (unsigned) dx;
  dx -= x0;
  x0 *= 2;
  // linear interpolation
  return m_table[x0]+dx*m_table[x0+1];
}

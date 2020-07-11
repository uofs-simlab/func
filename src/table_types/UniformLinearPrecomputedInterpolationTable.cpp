/* Implementation of a UniformPrecomputed Lookup table with linear interpolation */
#include "UniformLinearPrecomputedInterpolationTable.hpp"

#define IMPL_NAME UniformLinearPrecomputedInterpolationTable
REGISTER_ULUT_IMPL(IMPL_NAME);

UniformLinearPrecomputedInterpolationTable::UniformLinearPrecomputedInterpolationTable(FunctionContainer *func_container, UniformLookupTableParameters par) :
  UniformLookupTable(func_container, par)
{

  /* Base class default variables */
  m_name = STR(IMPL_NAME);
  m_order = 2;
  m_numTableEntries = m_numIntervals+1;
  m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;

  /* Allocate and set table */
  m_table.reset(new polynomial<2,16>[m_numTableEntries]);
  for (int ii=0;ii<m_numIntervals;++ii) {
    double x = m_minArg + ii*m_stepSize;
    m_grid[ii] = x;
    m_table[ii].coefs[0] = (mp_func)(x);
    x = m_minArg + (ii+1)*m_stepSize;
    m_table[ii].coefs[1] = (mp_func)(x) - m_table[ii].coefs[0];
  }
}

double UniformLinearPrecomputedInterpolationTable::operator()(double x)
{
  // nondimensionalized x position, scaled by step size
  double   dx = m_stepSize_inv*(x-m_minArg);
  // index of previous table entry
  unsigned x0  = (unsigned) dx;
  dx -= x0;
  // linear interpolation
  return m_table[x0].coefs[0]+dx*m_table[x0].coefs[1];
}

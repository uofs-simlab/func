/* Implementation of a Uniform Lookup table with linear interpolation */
#include "NonUniformLinearInterpolationTable.hpp"

#define IMPL_NAME NonUniformLinearInterpolationTable
REGISTER_ULUT_IMPL(IMPL_NAME);

NonUniformLinearInterpolationTable::NonUniformLinearInterpolationTable(FunctionContainer *func_container, UniformLookupTableParameters par) :
  NonUniformLookupTable(func_container, par)
{
  /* Base class variables */
  m_name = STR(IMPL_NAME);
  m_order = 2;
  m_numTableEntries = m_numIntervals;
  m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;

  /* Allocate and set table */
  m_table.reset(new polynomial<1,8>[m_numTableEntries]);
  for (int ii=0; ii<m_numTableEntries; ++ii) {
    const double x = m_minArg + m_g(ii/(m_numTableEntries-1))*(m_maxArg-m_minArg);
    m_grid[ii]  = x;
    m_table[ii].coefs[0] = mp_func(x);
  }
}

double NonUniformLinearInterpolationTable::operator()(double x)
{
  double   dx = (m_numTableEntries-1)*m_g_inv((x-m_minArg)/(m_maxArg-m_minArg));
  // index of previous table entry
  unsigned x1  = (unsigned) dx;
  // value of table entries around x position
  dx -= x1;
  double   y1  = m_table[x1].coefs[0];
  double   y2  = m_table[x1+1].coefs[0];
  // linear interpolation
  return y1+dx*(y2-y1);
}

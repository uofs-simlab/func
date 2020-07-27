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
    const double x = m_minArg + m_transferFunction->g(ii/(double)(m_numTableEntries-1))*(m_maxArg-m_minArg);
    m_grid[ii]  = x;
    m_table[ii].coefs[0] = mp_func(x);
  }
}

double NonUniformLinearInterpolationTable::operator()(double x)
{
  // TODO simplify I hope
  double x_idx = (unsigned) (m_numTableEntries-1)*m_transferFunction->g_inv((x-m_minArg)/(m_maxArg-m_minArg));
  double h     = m_grid[x_idx+1] - m_grid[x_idx];
  double dx    = (x - m_grid[x_idx])/h;

  // value of table entries around x position
  double   y1  = m_table[x_idx].coefs[0];
  double   y2  = m_table[x_idx+1].coefs[0];
  // linear interpolation
  return y1+dx*(y2-y1);
}

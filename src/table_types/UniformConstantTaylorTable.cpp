/* Implementation of a Uniform Lookup table with linear interpolation */
#include "UniformConstantTaylorTable.hpp"

#define IMPL_NAME UniformConstantTaylorTable
REGISTER_ULUT_IMPL(IMPL_NAME);

UniformConstantTaylorTable::UniformConstantTaylorTable(FunctionContainer *func_container, UniformLookupTableParameters par) : 
  UniformLookupTable(func_container, par)
{

  /* Base class default variables */
  m_name = STR(IMPL_NAME);
  m_order = 1;
  m_numTableEntries = m_numIntervals;
  m_dataSize = (unsigned) sizeof(double) * m_numTableEntries;

  /* Allocate and set table */
  m_table.reset(new double[m_numTableEntries]);
  for (int ii=0; ii<m_numTableEntries; ++ii) {
    double x = m_minArg + ii*m_stepSize;
    m_grid[ii] = x;
    m_table[ii] = (*mp_func)(x);
  }
}

double UniformConstantTaylorTable::operator()(double x)
{
  // Constant interpolation from table point immediately below x
  return m_table[(unsigned)((x-m_minArg)/m_stepSize+0.5)];
}

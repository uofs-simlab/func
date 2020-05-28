/* Implementation of a Uniform Lookup table with linear interpolation */
#include "UniformLinearInterpolationTable.hpp"

#define IMPL_NAME UniformLinearInterpolationTable
REGISTER_ULUT_IMPL(IMPL_NAME);

UniformLinearInterpolationTable::UniformLinearInterpolationTable(EvaluationFunctor<double,double> *func, UniformLookupTableParameters par) : UniformLookupTable(func, par)
{
  /* Base class variables */
  m_name = STR(IMPL_NAME);
  m_order = 2;
  m_numTableEntries = m_numIntervals; // should this be m_numIntervals+1?
  m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;

  /* Allocate and set table */
  m_table.reset(new polynomial<1,8>[m_numTableEntries]);
  for (int ii=0; ii<m_numTableEntries; ++ii) {
    const double x = m_minArg + ii*m_stepSize;
    m_grid[ii]  = x;
    m_table[ii].coefs[0] = (*mp_func)(x);
  }
}

double UniformLinearInterpolationTable::operator()(double x)
{
  // nondimensionalized x position, scaled by step size
  double   dx = (x-this->m_minArg)/m_stepSize;
  // index of previous table entry
  unsigned x1  = (unsigned) dx;
  // value of table entries around x position
  dx -= x1;
  double   y1  = m_table[x1].coefs[0];
  double   y2  = m_table[x1+1].coefs[0];
  // linear interpolation
  return y1+dx*(y2-y1);
}

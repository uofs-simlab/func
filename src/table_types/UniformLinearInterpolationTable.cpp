/* Implementation of a Uniform Lookup table with linear interpolation */
#include "UniformLinearInterpolationTable.hpp"

#define IMPL_NAME UniformLinearInterpolationTable
REGISTER_ULUT_IMPL(IMPL_NAME);

UniformLinearInterpolationTable::UniformLinearInterpolationTable(EvaluationFunctor<double,double> *func, UniformLookupTableParameters par) : UniformLookupTable(func, par)
{
  /* Base class variables */
  m_name = STR(IMPL_NAME);
  m_order = 2;
  m_numTableEntries = m_numIntervals;
  m_dataSize = (unsigned) sizeof(double) * m_numTableEntries;

  /* Allocate and set table */
  m_table.reset(new double[m_numTableEntries]);
  for (int ii=0; ii<m_numTableEntries; ++ii) {
    const double x    = m_minArg + ii*m_stepSize;
    m_grid[ii]  = x;
    m_table[ii] = (*mp_func)(x);
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
  double   y1  = m_table[x1];
  double   y2  = m_table[x1+1];
  // linear interpolation
  return y1+dx*(y2-y1);
}

/* Implementation of a UniformPrecomputed Lookup table with linear interpolation */
#include "UniformCubicPrecomputedInterpolationTable.hpp"

#define IMPL_NAME UniformCubicPrecomputedInterpolationTable
REGISTER_ULUT_IMPL(IMPL_NAME);

UniformCubicPrecomputedInterpolationTable::UniformCubicPrecomputedInterpolationTable(EvaluationFunctor<double,double> *func, UniformLookupTableParameters par) : UniformLookupTable(func, par)
{

  /* Base class default variables */
  m_name = STR(IMPL_NAME);
  m_order = 4;
  m_numTableEntries = m_numIntervals+1;
  m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;

  /* Allocate and set table */
  m_table.reset(new polynomial<4,32>[m_numTableEntries]);
  for (int ii=0;ii<m_numIntervals;++ii) {
    const double x = m_minArg + ii*m_stepSize;
    // grid points
    m_grid[ii] = x;
    // polynomial coefficients
    const double y0 = (*mp_func)(x);
    const double y1 = (*mp_func)(x+m_stepSize/3);
    const double y2 = (*mp_func)(x+2*m_stepSize/3);
    const double y3 = (*mp_func)(x+m_stepSize);
    m_table[ii].coefs[0]   = y0;
    m_table[ii].coefs[1] = -11*y0/2+9*y1-9*y2/2+y3;
    m_table[ii].coefs[2] = 9*y0-45*y1/2+18*y2-9*y3/2;
    m_table[ii].coefs[3] = -9*y0/2+27*y1/2-27*y2/2+9*y3/2;
  }
}

double UniformCubicPrecomputedInterpolationTable::operator()(double x)
{
  // nondimensionalized x position, scaled by step size
  double   dx = m_stepSize_inv*(x-this->m_minArg);
  // index of previous table entry
  // unsigned x0  = (unsigned) floor(dx);
  unsigned x0  = (unsigned) dx;
  // value of table entries around x position
  dx -= x0;
  // cubic interpolation
  return m_table[x0].coefs[0]+dx*(m_table[x0].coefs[1]+dx*(m_table[x0].coefs[2]+dx*m_table[x0].coefs[3]));
}

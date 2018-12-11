/* Implementation of a UniformPrecomputed Lookup table with linear interpolation */
#include "UniformQuadraticPrecomputedInterpolationTable.hpp"

#define IMPL_NAME UniformQuadraticPrecomputedInterpolationTable
REGISTER_ULUT_IMPL(IMPL_NAME);


UniformQuadraticPrecomputedInterpolationTable::UniformQuadraticPrecomputedInterpolationTable(EvaluationFunctor<double,double> *func, UniformLookupTableParameters par) : UniformLookupTable(func, par)
{

  /* Base class default variables */
  m_name = STR(IMPL_NAME);
  m_order = 3;
  m_numTableEntries = 3*(m_numIntervals+1);
  m_dataSize = (unsigned) sizeof(double) * m_numTableEntries;

  /* Allocate and set table */
  m_table.reset(new double[m_numTableEntries]);
  for (int ii=0;ii<m_numIntervals;++ii) {
    const double x = m_minArg + ii*m_stepSize;
    // grid points
    m_grid[ii] = x;
    // polynomial coefficients
    const double y0 = (*mp_func)(x);
    const double y1 = (*mp_func)(x+m_stepSize/2);
    const double y2 = (*mp_func)(x+m_stepSize);
    m_table[3*ii]   = y0;
    m_table[3*ii+1] = -3*y0+4*y1-y2;
    m_table[3*ii+2] = 2*y0+-4*y1+2*y2;
  }
}

double UniformQuadraticPrecomputedInterpolationTable::operator()(double x)
{
  // nondimensionalized x position, scaled by step size
  double   dx = m_stepSize_inv*(x-this->m_minArg);
  // index of previous table entry
  // unsigned x0  = (unsigned) floor(dx);
  unsigned x0  = (unsigned) dx;
  // value of table entries around x position
  dx -= x0;
  x0 *= 3;
  // quadratic interpolation
  return m_table[x0]+dx*(m_table[x0+1]+dx*m_table[x0+2]);
}

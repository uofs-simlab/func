/* Implementation of a UniformPrecomputed Lookup table with linear interpolation */
#include "UniformCubicHermiteTable.hpp"
#include <boost/math/differentiation/autodiff.hpp>

#define IMPL_NAME UniformCubicHermiteTable
REGISTER_ULUT_IMPL(IMPL_NAME);

UniformCubicHermiteTable::UniformCubicHermiteTable(FunctionContainer *func_container, UniformLookupTableParameters par) :
  UniformLookupTable(func_container, par)
{
  using boost::math::differentiation::make_fvar;
  /* Base class default variables */
  m_name = STR(IMPL_NAME);
  m_order = 4;
  m_numTableEntries = m_numIntervals+1;
  m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;

  __IS_NULL(func_container->autodiff1_func);
  mp_boost_func = func_container->autodiff1_func;

  /* Allocate and set table */
  m_table.reset(new polynomial<4,32>[m_numTableEntries]);
  for (int ii=0; ii<m_numIntervals; ++ii) {
    const double x = m_minArg + ii*m_stepSize;
    m_grid[ii] = x;

    const auto derivs0 = (mp_boost_func)(make_fvar<double,1>(x));
    const double y0    = derivs0.derivative(0);
    const double m0    = derivs0.derivative(1);
    const auto derivs1 = (mp_boost_func)(make_fvar<double,1>(x+m_stepSize));
    const double y1    = derivs1.derivative(0);
    const double m1    = derivs1.derivative(1);

    m_table[ii].coefs[0] = y0;
    m_table[ii].coefs[1] = m_stepSize*m0;
    m_table[ii].coefs[2] = -3*y0+3*y1-(2*m0+m1)*m_stepSize;
    m_table[ii].coefs[3] = 2*y0-2*y1+(m0+m1)*m_stepSize;
  }
}

double UniformCubicHermiteTable::operator()(double x)
{
  // nondimensionalized x position, scaled by step size
  double   dx = m_stepSize_inv*(x-m_minArg);
  // index of previous table entry
  unsigned x0  = (unsigned) dx;
  dx -= x0;
  // linear interpolation
  return m_table[x0].coefs[0]+dx*(m_table[x0].coefs[1]+dx*(m_table[x0].coefs[2]+dx*m_table[x0].coefs[3]));
}

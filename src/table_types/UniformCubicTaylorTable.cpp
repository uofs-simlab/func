/* Implementation of a Uniform Lookup table with linear interpolation */
#include "UniformCubicTaylorTable.hpp"

using boost::math::differentiation::make_fvar;

#define IMPL_NAME UniformCubicTaylorTable
REGISTER_ULUT_IMPL(IMPL_NAME);

UniformCubicTaylorTable::UniformCubicTaylorTable(FunctionContainer *func_container, UniformLookupTableParameters par) :
  UniformLookupTable(func_container, par)
{
  /* Base class default variables */
  m_name = STR(IMPL_NAME);
  m_order = 4;
  m_numTableEntries = m_numIntervals;
  m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;

  __IS_NULLPTR(func_container->fvar3_func);
  mp_boost_func = func_container->fvar3_func;

  /* Allocate and set table */
  m_table.reset(new polynomial<4,32>[m_numTableEntries]);
  for (int ii=0;ii<m_numIntervals;++ii) {
    double x = (m_minArg + ii*m_stepSize);
    m_grid[ii] = x;
    auto const derivs = (*mp_boost_func)(make_fvar<double,3>(x));

    m_table[ii].coefs[0] = derivs.derivative(0);
    m_table[ii].coefs[1] = derivs.derivative(1);
    m_table[ii].coefs[2] = derivs.derivative(2)/2;
    m_table[ii].coefs[3] = derivs.derivative(3)/6;
  }
}

double UniformCubicTaylorTable::operator()(double x)
{
  // nondimensionalized x position
  double  dx  = (x-m_minArg);
  double  x1r = dx/m_stepSize+0.5;
  // index of previous table entry
  unsigned x1 = (unsigned) x1r;
  dx -= x1*m_stepSize;
  // Cubic Taylor series from grid point below x
  return m_table[x1].coefs[0]+dx*(m_table[x1].coefs[1]+dx*(m_table[x1].coefs[2]+dx*m_table[x1].coefs[3]));
}

/* Implementation of a Uniform Lookup table with linear interpolation */
#include "UniformLinearTaylorTable.hpp"

#define IMPL_NAME UniformLinearTaylorTable
REGISTER_ULUT_IMPL(IMPL_NAME);

UniformLinearTaylorTable::UniformLinearTaylorTable(FunctionContainer *func_container, UniformLookupTableParameters par) : 
  UniformLookupTable(func_container, par)
{
  using boost::math::differentiation::make_fvar;

  /* Base class variables */
  m_name = STR(IMPL_NAME);
  m_order = 2;
  m_numTableEntries = m_numIntervals;
  m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;

  __IS_NULL(func_container->autodiff1_func);
  mp_boost_func = func_container->autodiff1_func;

  /* Allocate and set table */
  m_table.reset(new polynomial<2,16>[m_numTableEntries]);
  for (int ii=0; ii<m_numIntervals; ++ii) {
    const double x = m_minArg + ii*m_stepSize;
    m_grid[ii]     = x;
    // get every derivative up to the first
    auto const derivs = (mp_boost_func)(make_fvar<double,1>(x));
    m_table[ii].coefs[0] = derivs.derivative(0);
    m_table[ii].coefs[1] = derivs.derivative(1);
  }
}

double UniformLinearTaylorTable::operator()(double x)
{
  // nondimensionalized x position
  double  dx  = (x-m_minArg);
  double  x1r = dx/m_stepSize+0.5;
  // index of previous table entry
  unsigned x1 = (unsigned) x1r;
  dx -= x1*m_stepSize;
  // linear Taylor series from grid point below x
  return m_table[x1].coefs[0]+m_table[x1].coefs[1]*dx;
}

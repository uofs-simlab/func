/* Implementation of a Uniform Lookup table with linear interpolation */
#include "UniformQuadraticTaylorTable.hpp"

using boost::math::differentiation::make_fvar;

#define IMPL_NAME UniformQuadraticTaylorTable
REGISTER_ULUT_IMPL_DIFF(2,IMPL_NAME);

UniformQuadraticTaylorTable::UniformQuadraticTaylorTable(EvaluationFunctor<autodiff_fvar<double,2>,autodiff_fvar<double,2>> *func, UniformLookupTableParameters par) :
  UniformAutoDiffTable(func, par)
{
  /* Base class default variables */
  m_name = STR(IMPL_NAME);
  m_order = 3;
  m_numTableEntries = m_numIntervals;
  m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;

  /* Allocate and set table */
  m_table.reset(new polynomial<3,32>[m_numTableEntries]);
  for (int ii=0;ii<m_numIntervals;++ii) {
    double x = (m_minArg + ii*m_stepSize);
    m_grid[ii] = x;
    auto const y = (*mp_boost_func)(make_fvar<double,2>(x));
    m_table[ii].coefs[0] = y.derivative(0);
    m_table[ii].coefs[1] = y.derivative(1);
    m_table[ii].coefs[2] = y.derivative(2)/2;
  }
}

double UniformQuadraticTaylorTable::operator()(double x)
{
  // nondimensionalized x position
  double  dx  = (x-m_minArg);
  double  x1r = dx/m_stepSize+0.5;
  // index of previous table entry
  unsigned x1 = (unsigned) x1r;
  dx -= x1*m_stepSize;
  // Quadratic Taylor series from grid point below x
  return m_table[x1].coefs[0]+dx*(m_table[x1].coefs[1]+dx*m_table[x1].coefs[2]);
}

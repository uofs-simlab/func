/*
  Cubic Interpolation LUT with uniform sampling (precomputed coefficients)

  Usage example:
    UniformCubicPrecomputedInterpolationTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#include "UniformLookupTable.hpp"

template <typename IN_TYPE, typename OUT_TYPE>
class UniformCubicPrecomputedInterpolationTable final : public UniformLookupTable<IN_TYPE,OUT_TYPE>
{
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE);
  REGISTER_LUT(UniformCubicPrecomputedInterpolationTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,4>[]> m_table;

public:
  UniformCubicPrecomputedInterpolationTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, UniformLookupTableParameters<IN_TYPE> par) :
    UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, par)
  {
    /* Base class default variables */
    m_name = STR(UniformCubicPrecomputedInterpolationTable);
    m_order = 4;
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;

    /* Allocate and set table */
    m_table.reset(new polynomial<OUT_TYPE,4>[m_numTableEntries]);
    for (int ii=0;ii<m_numIntervals;++ii) {
      const IN_TYPE x = m_minArg + ii*m_stepSize;
      // grid points
      m_grid[ii] = x;
      // polynomial coefficients
      const OUT_TYPE y0 = (mp_func)(x);
      const OUT_TYPE y1 = (mp_func)(x+m_stepSize/3);
      const OUT_TYPE y2 = (mp_func)(x+2*m_stepSize/3);
      const OUT_TYPE y3 = (mp_func)(x+m_stepSize);
      m_table[ii].coefs[0] = y0;
      m_table[ii].coefs[1] = -11*y0/2+9*y1-9*y2/2+y3;
      m_table[ii].coefs[2] = 9*y0-45*y1/2+18*y2-9*y3/2;
      m_table[ii].coefs[3] = -9*y0/2+27*y1/2-27*y2/2+9*y3/2;
    }
  }

  OUT_TYPE operator()(IN_TYPE x) override
  {
    // nondimensionalized x position, scaled by step size
    OUT_TYPE  dx = m_stepSize_inv*(x-m_minArg);
    // index of previous table entry
    // unsigned x0  = (unsigned) floor(dx);
    unsigned x0  = (unsigned) dx;
    // value of table entries around x position
    dx -= x0;
    // cubic interpolation
    return m_table[x0].coefs[0]+dx*(m_table[x0].coefs[1]+dx*(m_table[x0].coefs[2]+dx*m_table[x0].coefs[3]));
  }
};

REGISTER_DOUBLE_AND_FLOAT_LUT_IMPLS(UniformCubicPrecomputedInterpolationTable);

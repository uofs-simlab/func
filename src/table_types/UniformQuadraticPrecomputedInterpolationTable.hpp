/*
  Quadratic Interpolation LUT with uniform sampling (precomputed coefficients)

  Usage example:
    UniformQuadraticPrecomputedInterpolationTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "UniformLookupTable.hpp"

template <typename IN_TYPE, typename OUT_TYPE>
class UniformQuadraticPrecomputedInterpolationTable final : public UniformLookupTable<IN_TYPE,OUT_TYPE>
{
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE);
  REGISTER_LUT(UniformQuadraticPrecomputedInterpolationTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,3>[]> m_table;

public:
  UniformQuadraticPrecomputedInterpolationTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, UniformLookupTableParameters<IN_TYPE> par) :
    UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, par)
  {
    /* Base class default variables */
    m_name = STR(UniformQuadraticPrecomputedInterpolationTable);
    m_order = 3;
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;

    /* Allocate and set table */
    m_table.reset(new polynomial<OUT_TYPE,3>[m_numTableEntries]);
    for (int ii=0;ii<m_numIntervals;++ii) {
      const IN_TYPE x = m_minArg + ii*m_stepSize;
      // grid points
      m_grid[ii] = x;
      // polynomial coefficients
      const OUT_TYPE y0  = (mp_func)(x);
      const OUT_TYPE y1  = (mp_func)(x+m_stepSize/2);
      const OUT_TYPE y2  = (mp_func)(x+m_stepSize);
      m_table[ii].coefs[0] = y0;
      m_table[ii].coefs[1] = -3*y0+4*y1-y2;
      m_table[ii].coefs[2] = 2*y0+-4*y1+2*y2;
    }
  }

  OUT_TYPE operator()(IN_TYPE x) override
  {
    // nondimensionalized x position, scaled by step size
    OUT_TYPE dx = m_stepSize_inv*(x-m_minArg);
    // index of previous table entry
    unsigned x0  = (unsigned) dx;
    // value of table entries around x position
    dx -= x0;
    // quadratic interpolation
    return m_table[x0].coefs[0]+dx*(m_table[x0].coefs[1]+dx*m_table[x0].coefs[2]);
  }
};

REGISTER_DOUBLE_AND_FLOAT_LUT_IMPLS(UniformQuadraticPrecomputedInterpolationTable);

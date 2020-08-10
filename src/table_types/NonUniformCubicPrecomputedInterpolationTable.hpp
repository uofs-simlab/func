/*
  Cubic polynomial Interpolation LUT with nonuniform sampling

  Usage example:
    NonUniformCubicPrecomputedInterpolationTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "NonUniformLookupTable.hpp"

template <typename IN_TYPE, typename OUT_TYPE>
class NonUniformCubicPrecomputedInterpolationTable final : public NonUniformLookupTable<IN_TYPE,OUT_TYPE>
{
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE);
  INHERIT_NONUNIFORM_LUT(IN_TYPE,OUT_TYPE);

  REGISTER_LUT(NonUniformCubicPrecomputedInterpolationTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,4>[]> m_table;
public:
  NonUniformCubicPrecomputedInterpolationTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container,
      UniformLookupTableParameters<IN_TYPE> par) :
    NonUniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, par)
  {
    /* Base class variables */
    m_name = STR(NonUniformCubicPrecomputedInterpolationTable);
    m_order = 4;
    m_numTableEntries = m_numIntervals + 1;
    m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;
    
    /* Allocate and set table */
    m_table.reset(new polynomial<OUT_TYPE,4>[m_numTableEntries]);
    for (int ii=0; ii<m_numIntervals; ++ii) {
      const IN_TYPE x = m_minArg + m_transferFunction->g((IN_TYPE) (ii/(double)(m_numIntervals-1)))*(m_maxArg-m_minArg);
      // the local stepsize of our nonuniform grid
      const IN_TYPE h = m_minArg + m_transferFunction->g((IN_TYPE) ((ii+1)/(double)(m_numIntervals-1)))*(m_maxArg-m_minArg) - x;
      // grid points
      m_grid[ii] = x;
      // polynomial coefficients
      const OUT_TYPE y0 = mp_func(x);
      const OUT_TYPE y1 = mp_func(x+h/3);
      const OUT_TYPE y2 = mp_func(x+2*h/3);
      const OUT_TYPE y3 = mp_func(x+h);
      m_table[ii].coefs[0] = y0;
      m_table[ii].coefs[1] = -11*y0/2+9*y1-9*y2/2+y3;
      m_table[ii].coefs[2] = 9*y0-45*y1/2+18*y2-9*y3/2;
      m_table[ii].coefs[3] = -9*y0/2+27*y1/2-27*y2/2+9*y3/2;
    }
  }

  OUT_TYPE operator()(IN_TYPE x) override
  {
    // find the subinterval x lives in
    unsigned x_idx = (unsigned) (m_numIntervals-1)*m_transferFunction->g_inv((x-m_minArg)/(m_maxArg-m_minArg));
    //if(x < m_grid[x_idx]-0.00000005 || m_grid[x_idx+1] < x)
    //  std::cerr << "The hash thinks " << x << " is in [" << m_grid[x_idx] << "," << m_grid[x_idx+1] << ")" << std::endl;

    // find where x is in that subinterval
    double h  = m_grid[x_idx+1] - m_grid[x_idx];
    double dx = (x - m_grid[x_idx])/h;

    // cubic interpolation
    return m_table[x_idx].coefs[0]+dx*(m_table[x_idx].coefs[1]+dx*(m_table[x_idx].coefs[2]+dx*m_table[x_idx].coefs[3]));
  }
};

REGISTER_DOUBLE_AND_FLOAT_LUT_IMPLS(NonUniformCubicPrecomputedInterpolationTable);

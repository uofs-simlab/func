/*
  Linear Interpolation LUT with nonuniform sampling

  Usage example:
    NonUniformLinearInterpolationTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "NonUniformLookupTable.hpp"

template <typename IN_TYPE, typename OUT_TYPE>
class NonUniformLinearInterpolationTable final : public NonUniformLookupTable<IN_TYPE,OUT_TYPE>
{
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE);
  INHERIT_NONUNIFORM_LUT(IN_TYPE,OUT_TYPE);

  REGISTER_LUT(NonUniformLinearInterpolationTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,1>[]> m_table;
public:
  NonUniformLinearInterpolationTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container,
      UniformLookupTableParameters<IN_TYPE> par) :
    NonUniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, par)
  {
    /* Base class variables */
    m_name = STR(NonUniformLinearInterpolationTable);
    m_order = 2;
    m_numTableEntries = m_numIntervals;
    m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;
    
    /* Allocate and set table */
    m_table.reset(new polynomial<OUT_TYPE,1>[m_numTableEntries]);
    for (int ii=0; ii<m_numIntervals; ++ii) {
      const IN_TYPE x = m_minArg + m_transferFunction->g((IN_TYPE) (ii/(double)(m_numIntervals-1)))*(m_maxArg-m_minArg);
      m_grid[ii]  = x;
      m_table[ii].coefs[0] = mp_func(x);
    }
  }

  OUT_TYPE operator()(IN_TYPE x) override
  {
    // find the subinterval x lives in
    unsigned x_idx = (unsigned) (m_numTableEntries-1)*m_transferFunction->g_inv((x-m_minArg)/(m_maxArg-m_minArg));
    //if(x < m_grid[x_idx]-0.00000005 || m_grid[x_idx+1] < x)
    //  std::cerr << "The hash thinks " << x << " is in [" << m_grid[x_idx] << "," << m_grid[x_idx+1] << ")" << std::endl;

    // find where x is in that subinterval
    IN_TYPE h   = m_grid[x_idx+1] - m_grid[x_idx];
    OUT_TYPE dx = (OUT_TYPE) (x - m_grid[x_idx])/h;

    // value of table entries around x position
    double   y1  = m_table[x_idx].coefs[0];
    double   y2  = m_table[x_idx+1].coefs[0];
    // linear interpolation
    return y1+dx*(y2-y1);
  }
};

REGISTER_DOUBLE_AND_FLOAT_LUT_IMPLS(NonUniformLinearInterpolationTable);

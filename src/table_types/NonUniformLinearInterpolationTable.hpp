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

template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE, class TRANSFER_FUNC_TYPE = TransferFunctionSinh<IN_TYPE>>
class NonUniformLinearInterpolationTable final : public NonUniformLookupTable<IN_TYPE,OUT_TYPE,TRANSFER_FUNC_TYPE>
{
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE);
  INHERIT_NONUNIFORM_LUT(IN_TYPE,OUT_TYPE);

  REGISTER_LUT(NonUniformLinearInterpolationTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,1>[]> m_table;
public:
  NonUniformLinearInterpolationTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container,
      UniformLookupTableParameters<IN_TYPE> par) :
    NonUniformLookupTable<IN_TYPE,OUT_TYPE,TRANSFER_FUNC_TYPE>(func_container, par)
  {
    /* Base class variables */
    m_name = STR(NonUniformLinearInterpolationTable);
    m_order = 2;
    m_numTableEntries = m_numIntervals;
    m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;
    
    /* Allocate and set table */
    m_table.reset(new polynomial<OUT_TYPE,1>[m_numTableEntries]);
    for (int ii=0; ii<m_numIntervals; ++ii) {
      // transform the previously used uniform grid to a nonuniform grid
      const IN_TYPE x = m_transferFunction->g(m_minArg + ii*m_stepSize);
      m_grid[ii]  = x;
      m_table[ii].coefs[0] = mp_func(x);
    }
  }

  OUT_TYPE operator()(IN_TYPE x) override
  {
    // set x_idx = floor((g_inv(x)-m_minArg)/m_stepSize)
    // where each of the above member vars are encoded into g_inv
    // Note: the fractional part of g_inv can't be used as
    // dx b/c distances work differently in this nonuniform grid
    //unsigned int x_idx = (unsigned int) (m_transferFunction->g_inv(x));

    // find where x is in that subinterval
    //IN_TYPE  h  = m_grid[x_idx+1] - m_grid[x_idx];
    //IN_TYPE dx = (OUT_TYPE) (x - m_grid[x_idx])/h;

    OUT_TYPE dx = m_transferFunction->g_inv(x);
    unsigned x_idx = (unsigned) dx;
    dx -= x_idx;

    // value of table entries around x position
    OUT_TYPE y1  = m_table[x_idx].coefs[0];
    OUT_TYPE y2  = m_table[x_idx+1].coefs[0];
    // linear interpolation
    return y1+dx*(y2-y1);
  }
};

REGISTER_NONUNIFORM_IMPL(NonUniformLinearInterpolationTable,double,double,TransferFunctionSinh<double>);

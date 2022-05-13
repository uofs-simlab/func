/*
  Constant Taylor LUT with uniform sampling

  Usage example:
    UniformConstantTaylorTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "UniformLookupTable.hpp"

template <typename IN_TYPE, typename OUT_TYPE>
class UniformConstantTaylorTable final : public UniformLookupTable<IN_TYPE,OUT_TYPE>
{
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE);

  REGISTER_LUT(UniformConstantTaylorTable);

  __attribute__((aligned)) std::unique_ptr<OUT_TYPE[]> m_table;
public:
  UniformConstantTaylorTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, UniformLookupTableParameters<IN_TYPE> par) :
    UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, par)
  {
    /* Base class default variables */
    m_name = STR(UniformConstantTaylorTable);
    m_order = 1;
    m_numTableEntries = m_numIntervals;
    m_dataSize = (unsigned) sizeof(OUT_TYPE) * m_numTableEntries;

    /* Allocate and set table */
    m_table.reset(new OUT_TYPE[m_numTableEntries]);
    for (int ii=0; ii<m_numTableEntries; ++ii) {
      IN_TYPE x = m_minArg + ii*m_stepSize;
      m_grid[ii] = x;
      m_table[ii] = mp_func(x);
    }
  }

  /* Constant interpolation from table point immediately below x */
  OUT_TYPE operator()(IN_TYPE x) override
  {
    return m_table[(unsigned)((x-m_minArg)/m_stepSize+0.5)];
  }
};

REGISTER_DOUBLE_AND_FLOAT_LUT_IMPLS(UniformConstantTaylorTable);

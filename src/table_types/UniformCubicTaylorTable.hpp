/*
  Cubic Taylor LUT with uniform sampling

  Usage example:
    UniformCubicTaylorTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "UniformLookupTable.hpp"

template <typename IN_TYPE, typename OUT_TYPE>
class UniformCubicTaylorTable final : public UniformLookupTable<IN_TYPE,OUT_TYPE>
{
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE);
  REGISTER_LUT(UniformCubicTaylorTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,4>[]> m_table;
  std::function<fvar<OUT_TYPE,3>(fvar<OUT_TYPE,3>)> mp_boost_func;

public:
  UniformCubicTaylorTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, UniformLookupTableParameters<IN_TYPE> par) :
    UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, par)
  {
    using boost::math::differentiation::make_fvar;

    /* Base class default variables */
    m_name = STR(UniformCubicTaylorTable);
    m_order = 4;
    m_numTableEntries = m_numIntervals;
    m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;

    __IS_NULL(func_container->autodiff3_func);
    mp_boost_func = func_container->autodiff3_func;

    /* Allocate and set table */
    m_table.reset(new polynomial<OUT_TYPE,4>[m_numTableEntries]);
    for (int ii=0;ii<m_numIntervals;++ii) {
      IN_TYPE x = (m_minArg + ii*m_stepSize);
      m_grid[ii] = x;
      auto const derivs = (mp_boost_func)(make_fvar<IN_TYPE,3>(x));

      m_table[ii].coefs[0] = derivs.derivative(0);
      m_table[ii].coefs[1] = derivs.derivative(1);
      m_table[ii].coefs[2] = derivs.derivative(2)/2;
      m_table[ii].coefs[3] = derivs.derivative(3)/6;
    }
  }

  OUT_TYPE operator()(IN_TYPE x) override
  {
    // nondimensionalized x position
    OUT_TYPE dx  = (x-m_minArg);
    OUT_TYPE x1r = dx/m_stepSize+0.5;
    // index of previous table entry
    unsigned  x1 = (unsigned) x1r;
    dx -= x1*m_stepSize;
    // Cubic Taylor series from grid point below x
    return m_table[x1].coefs[0]+dx*(m_table[x1].coefs[1]+dx*(m_table[x1].coefs[2]+dx*m_table[x1].coefs[3])); 
  }

  std::function<fvar<OUT_TYPE,3>(fvar<OUT_TYPE,3>)> boost_function(){ return mp_boost_func; }
};

REGISTER_DOUBLE_AND_FLOAT_LUT_IMPLS(UniformCubicTaylorTable);

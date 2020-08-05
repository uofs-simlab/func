/*
  Quadratic Taylor LUT with uniform sampling

  Usage example:
    UniformQuadraticTaylorTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "UniformLookupTable.hpp"

template <typename IN_TYPE, typename OUT_TYPE>
class UniformQuadraticTaylorTable final : public UniformLookupTable<IN_TYPE,OUT_TYPE>
{
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE);
  REGISTER_LUT(UniformQuadraticTaylorTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,3>[]> m_table;
  std::function<fvar<OUT_TYPE,2>(fvar<IN_TYPE,2>)> mp_boost_func;

public:
  UniformQuadraticTaylorTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, UniformLookupTableParameters<IN_TYPE> par) :
    UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, par)
  {
    using boost::math::differentiation::make_fvar;

    /* Base class default variables */
    m_name = STR(UniformQuadraticTaylorTable);
    m_order = 3;
    m_numTableEntries = m_numIntervals;
    m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;

    __IS_NULL(func_container->autodiff2_func);
    mp_boost_func = func_container->autodiff2_func;

    /* Allocate and set table */
    m_table.reset(new polynomial<OUT_TYPE,3>[m_numTableEntries]);
    for (int ii=0;ii<m_numIntervals;++ii) {
      IN_TYPE x = (m_minArg + ii*m_stepSize);
      m_grid[ii] = x;
      auto const derivs = (mp_boost_func)(make_fvar<IN_TYPE,2>(x));
      m_table[ii].coefs[0] = derivs.derivative(0);
      m_table[ii].coefs[1] = derivs.derivative(1);
      m_table[ii].coefs[2] = derivs.derivative(2)/2;
    }
  }

  OUT_TYPE operator()(IN_TYPE x) override
  {
    // nondimensionalized x position
    OUT_TYPE dx  = (x-m_minArg);
    OUT_TYPE x1r = dx/m_stepSize+0.5;
    // index of previous table entry
    unsigned x1 = (unsigned) x1r;
    dx -= x1*m_stepSize;
    // Quadratic Taylor series from grid point below x
    return m_table[x1].coefs[0]+dx*(m_table[x1].coefs[1]+dx*m_table[x1].coefs[2]);
  }

  std::function<fvar<OUT_TYPE,2>(fvar<IN_TYPE,2>)> boost_function(){ return mp_boost_func; }
};

REGISTER_DOUBLE_AND_FLOAT_LUT_IMPLS(UniformQuadraticTaylorTable);

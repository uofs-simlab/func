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
  REGISTER_ULUT(UniformCubicTaylorTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,4,32>[]> m_table;
  std::function<fvar<OUT_TYPE,3>(fvar<OUT_TYPE,3>)> mp_boost_func;
public:
  UniformCubicTaylorTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, UniformLookupTableParameters<IN_TYPE> par);
  OUT_TYPE operator()(IN_TYPE x) override;
  std::function<fvar<OUT_TYPE,3>(fvar<OUT_TYPE,3>)> boost_function(){ return mp_boost_func; }
};

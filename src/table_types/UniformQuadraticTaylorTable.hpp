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
  REGISTER_ULUT(UniformQuadraticTaylorTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,3,32>[]> m_table;
  std::function<fvar<OUT_TYPE,2>(fvar<IN_TYPE,2>)> mp_boost_func;
public:
  UniformQuadraticTaylorTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container,
      UniformLookupTableParameters<IN_TYPE> par);
  OUT_TYPE operator()(IN_TYPE x) override;
  std::function<fvar<OUT_TYPE,2>(fvar<IN_TYPE,2>)> boost_function(){ return mp_boost_func; }
};

/*
  Linear Taylor LUT with uniform sampling

  Usage example:
    UniformLinearTaylorTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "UniformLookupTable.hpp"

class UniformLinearTaylorTable final : public UniformLookupTable
{
  REGISTER_ULUT(UniformLinearTaylorTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<2,16>[]> m_table;
  std::function<fvar1(fvar1)> mp_boost_func;
public:
  UniformLinearTaylorTable(FunctionContainer *func_container, UniformLookupTableParameters par);
  double operator()(double x) override;
  std::function<fvar1(fvar1)> boost_function(){ return mp_boost_func; }
};

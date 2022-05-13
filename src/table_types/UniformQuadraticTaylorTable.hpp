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

class UniformQuadraticTaylorTable final : public UniformLookupTable
{
  REGISTER_ULUT(UniformQuadraticTaylorTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<3,32>[]> m_table;
  std::function<fvar2(fvar2)> mp_boost_func;
public:
  UniformQuadraticTaylorTable(FunctionContainer *func_container, UniformLookupTableParameters par);
  double operator()(double x) override;
  std::function<fvar2(fvar2)> boost_function(){ return mp_boost_func; }
};

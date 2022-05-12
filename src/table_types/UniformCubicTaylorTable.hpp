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

class UniformCubicTaylorTable final : public UniformLookupTable
{
  REGISTER_ULUT(UniformCubicTaylorTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<4,32>[]> m_table;
  EvaluationFunctor<fvar3,fvar3> *mp_boost_func;
public:
  UniformCubicTaylorTable(FunctionContainer *func_container, UniformLookupTableParameters par);
  double operator()(double x) override;
  EvaluationFunctor<fvar3,fvar3> *boost_function(){ return mp_boost_func; }
};

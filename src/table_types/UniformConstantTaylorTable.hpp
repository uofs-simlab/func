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

class UniformConstantTaylorTable final : public UniformLookupTable
{
  REGISTER_ULUT(UniformConstantTaylorTable);

public:
  UniformConstantTaylorTable(EvaluationFunctor<double,double> *func, UniformLookupTableParameters par);
  double operator()(double x) override;
};

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

#include <memory>

class UniformQuadraticTaylorTable final : public UniformLookupTable
{
  REGISTER_ULUT(UniformQuadraticTaylorTable);

public:
  UniformQuadraticTaylorTable(EvaluationFunctor<double,double> *func, UniformLookupTableParameters par);
  double operator()(double x) override;
};

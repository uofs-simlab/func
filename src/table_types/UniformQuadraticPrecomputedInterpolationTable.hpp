/*
  Quadratic Interpolation LUT with uniform sampling (precomputed coefficients)

  Usage example:
    UniformQuadraticPrecomputedInterpolationTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "UniformLookupTable.hpp"

#include <cmath>

class UniformQuadraticPrecomputedInterpolationTable final : public UniformLookupTable
{
  REGISTER_ULUT(UniformQuadraticPrecomputedInterpolationTable);

public:
  UniformQuadraticPrecomputedInterpolationTable(EvaluationFunctor<double,double> *func, UniformLookupTableParameters par);
  double operator()(double x) override;
};

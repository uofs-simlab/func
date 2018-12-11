#pragma once
#include "EvaluationFunctor.hpp"
#include <cmath>
#define FUNC(x) (7.7-13.0287*log(x))
#define FUNCNAME "(7.7-13.0287*log(x))"

class MyFunction final : public EvaluationFunctor<double,double>
{
public:
  double operator()(double x) override { return FUNC(x); }
  double deriv(double x) override
  {
    return -13.0287/x;
  }
  double deriv2(double x) override
  {
    return 13.0287/x/x;
  }
  double deriv3(double x) override
  {
    return -26.0574/x/x/x;
  }
};

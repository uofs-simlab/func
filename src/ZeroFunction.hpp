#pragma once
#include "EvaluationFunctor.hpp"
#define FUNC(x) 0.0
#define FUNCNAME "Zero"
class ZeroFunction final : public EvaluationFunctor<double,double>
{
public:
  double operator()(double x) override { return FUNC(X); }
  double deriv(double x) override
  {
    return 0.0;
  }
  double deriv2(double x) override
  {
    return 0.0;
  }
  double deriv3(double x) override
  {
    return 0.0;
  }
};

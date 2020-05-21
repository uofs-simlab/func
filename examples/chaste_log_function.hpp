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
  double deriv4(double x) override
  {
    return 78.1722/x/x/x/x;
  }
  double deriv5(double x) override
  {
    return -312.6888/x/x/x/x/x;
  }
  double deriv6(double x) override
  {
    return 1563.444/x/x/x/x/x/x;
  }
  double deriv7(double x) override
  {
    return -9380.664/x/x/x/x/x/x/x;
  }
};

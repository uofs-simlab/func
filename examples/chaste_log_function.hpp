#pragma once
#include "EvaluationFunctor.hpp"
#include <cmath>
#define FUNC(x) (7.7-13.0287*log(x))
#define FUNCNAME "(7.7-13.0287*log(x))"

template <typename T>
class MyFunction final : public EvaluationFunctor<T,T>
{
public:
  T operator()(T x) override { return FUNC(x); }
  T deriv(T x) override
  {
    return -13.0287/x;
  }
  T deriv2(T x) override
  {
    return 13.0287/x/x;
  }
  T deriv3(T x) override
  {
    return -26.0574/x/x/x;
  }
  T deriv4(T x) override
  {
    return 78.1722/x/x/x/x;
  }
  T deriv5(T x) override
  {
    return -312.6888/x/x/x/x/x;
  }
  T deriv6(T x) override
  {
    return 1563.444/x/x/x/x/x/x;
  }
  T deriv7(T x) override
  {
    return -9380.664/x/x/x/x/x/x/x;
  }
};

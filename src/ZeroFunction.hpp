#pragma once
#include "EvaluationFunctor.hpp"
#define FUNC(x) 0.0
#define FUNCNAME "Zero"

template <typename T>
class ZeroFunction final : public EvaluationFunctor<T,T>
{
public:
  T operator()(T x) override { return FUNC(X); }
  T deriv(T x) override
  {
    return 0.0;
  }
  T deriv2(T x) override
  {
    return 0.0;
  }
  T deriv3(T x) override
  {
    return 0.0;
  }
};

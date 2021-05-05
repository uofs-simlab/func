#pragma once
#include "EvaluationFunctor.hpp"
#include <cmath>
#define FUNC(x) (7.7*x*exp(x)-13.0287*x*x*log(x)-x +13.0*x*x*x)
#define FUNCNAME "(7.7*x*exp(x)-13.0287*x*x*log(x)-x +13.0*x*x*x)"

template <typename T>
class MyFunction final : public EvaluationFunctor<T,T>
{
public:
  T operator()(T x) override { return FUNC(x); }
};

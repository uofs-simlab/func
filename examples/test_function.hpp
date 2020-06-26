#pragma once
#include "EvaluationFunctor.hpp"
#include <cmath>
#define FUNC(x) (7.7*x*exp(x)-13.0287*x*x*log(x)-x +13.0*x*x*x)
#define FUNCNAME "(7.7*x*exp(x)-13.0287*x*x*log(x)-x +13.0*x*x*x)"

class MyFunction final : public EvaluationFunctor<boost_fvar,boost_fvar>
{
public:
  boost_fvar operator()(boost_fvar x) override { return FUNC(x); }
  boost_fvar deriv(boost_fvar x) override
  {
    return 7.7*exp(x)*(1.0+x) - 13.0287*(2.0*x*log(x)+x) - 1.0 + 39.0*x*x;
  }
  boost_fvar deriv2(boost_fvar x) override
  {
    return 7.7*exp(x)*(2.0+x) - 13.0287*(2.0*log(x)+3.0) + 78.0*x;
  }
  boost_fvar deriv3(boost_fvar x) override
  {
    return 7.7*exp(x)*(3.0+x) - 13.0287*(2.0/x) + 78.0;
  }
};

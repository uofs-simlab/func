/*
  Virtual base functor for evaluating a mathematical function

  - needed for combining evaluation and derivative(s) evaluation so that
    Taylor and Interpolation tables can have the same interface
*/
#pragma once

// TODO use a better exception type
#define NOT_DEFINED(ORDER)			\
  throw ORDER" derivative not defined.";	\
  return 0;

template <class IN_TYPE, class OUT_TYPE>
class EvaluationFunctor
{
public:
  virtual OUT_TYPE operator()(IN_TYPE x) = 0;
  /* Raise an exception in these by default? */
  virtual OUT_TYPE deriv(IN_TYPE x){ NOT_DEFINED("1st"); };
  virtual OUT_TYPE deriv2(IN_TYPE x){ NOT_DEFINED("2nd"); };
  virtual OUT_TYPE deriv3(IN_TYPE x){ NOT_DEFINED("3rd"); };
  virtual OUT_TYPE deriv4(IN_TYPE x){ NOT_DEFINED("4th"); };
  virtual OUT_TYPE deriv5(IN_TYPE x){ NOT_DEFINED("5th"); };
};

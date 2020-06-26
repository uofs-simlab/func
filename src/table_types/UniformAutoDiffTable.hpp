/*
  Templated header only abstract class for uniform LUTs 
  using boosts automatic differentiation for construction
*/
#pragma once
#include "UniformLookupTable.hpp"

#include <memory>
#include <map>
#include <vector>
#include <functional>
#include <boost/math/differentiation/autodiff.hpp>

using boost::math::differentiation::autodiff_fvar;

/* This can be reworked when we move to std::function */
template <unsigned int N>
class FuncWrapper : public EvaluationFunctor<double,double>
{
  public:
    // TODO test if this can be simplified with friend + remove warning???
    EvaluationFunctor<autodiff_fvar<double,N>,autodiff_fvar<double,N>> *fvar_EvalFunctor;
    FuncWrapper(EvaluationFunctor<autodiff_fvar<double,N>,autodiff_fvar<double,N>> *new_EvalFunctor)
      { fvar_EvalFunctor = new_EvalFunctor; }
    double operator()(double x) override { return (double) (*fvar_EvalFunctor)(x); }
};

template <unsigned int N>
class UniformAutoDiffTable : public UniformLookupTable
{
  friend class func_wrapper;

protected:
  EvaluationFunctor<autodiff_fvar<double,N>, autodiff_fvar<double,N>> *mp_boost_func;

public:
  UniformAutoDiffTable(EvaluationFunctor<autodiff_fvar<double,N>,autodiff_fvar<double,N>> *func, UniformLookupTableParameters par) :
    mp_boost_func(func), UniformLookupTable(new FuncWrapper<N>(mp_boost_func),par){}

  virtual ~UniformAutoDiffTable(){};

  /* public access of protected data */
  EvaluationFunctor<autodiff_fvar<double,N>,autodiff_fvar<double,N>> *fvar_function(){ return mp_boost_func; }
};

/*
  Templated header only abstract class for uniform LUTs 
  using boosts automatic differentiation for construction
*/
#pragma once
#include "UniformLookupTable.hpp"

using boost::math::differentiation::autodiff_fvar;

template <unsigned int N>
class UniformAutoDiffTable : public UniformLookupTable
{
protected:
  EvaluationFunctor<autodiff_fvar<double,N>, autodiff_fvar<double,N>> *mp_boost_func;

public:
  UniformAutoDiffTable(FunctionContainer *func_container, UniformLookupTableParameters par) :
    UniformLookupTable(func_container, par)
  {
    // select, check for null, and set appropriate boost function with modification to FunctionContainer
  }

  virtual ~UniformAutoDiffTable(){};

  /* public access of protected data */
  EvaluationFunctor<autodiff_fvar<double,N>,autodiff_fvar<double,N>> *fvar_function(){ return mp_boost_func; }
};

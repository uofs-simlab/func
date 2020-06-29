/*
 * will soon rename this file to FunctionContainer
  Virtual base functor for evaluating a mathematical function along with
  the data structure used to pass this 

  - needed for combining evaluation and derivative(s) evaluation so that
    Taylor and Interpolation tables can have the same interface
  - easy (though obtuse) to use if mathematical functions are templated
  - copy and paste the following example code into a new file and use
    the command :%s/new_name/foo/g in vim or (TODO something else) in emacs 
    to write your own function.

  Example usage:
  #include EvaluationFunctor.hpp
  // this could just be a function when we move to std::function
  template <typename T>
  struct foo : EvaluationFunctor
  {
    T operator()(T x){ return x; }
  };

  int main(){
    // these 'new' memory allocations can be removed when we move to std::function
    FunctionContainer foo_container{new foo<double>, new foo<fvar1>,
        new foo<fvar2>, new foo<fvar3>, new foo<fvar4>, 
        new foo<fvar5>, new foo<fvar6>, new foo<fvar7>};
    // or just
    // FunctionContainer foo_container;
    // foo_container.double_func = new foo<double>;
    // if you're not using any Taylor/Pade/Hermite tables
    // (the compiler will warn you)
    return 0;
  }
*/
#pragma once
#include <stdexcept>
#include <memory>
#include <boost/math/differentiation/autodiff.hpp>

using boost::math::differentiation::autodiff_fvar;
typedef autodiff_fvar<double,1> fvar1;
typedef autodiff_fvar<double,2> fvar2;
typedef autodiff_fvar<double,3> fvar3;
typedef autodiff_fvar<double,4> fvar4;
typedef autodiff_fvar<double,5> fvar5;
typedef autodiff_fvar<double,6> fvar6;
typedef autodiff_fvar<double,7> fvar7;

/* Used by each table type to check if the required function
 * type has been provided.
 * TODO decide whether or not to make this null for -DNDEBUG */
#define __IS_NULLPTR(VAR) \
    if(VAR == nullptr)     \
      throw std::invalid_argument(#VAR " is a nullptr")

template <class IN_TYPE, class OUT_TYPE>
class EvaluationFunctor
{
public:
  virtual OUT_TYPE operator()(IN_TYPE x) = 0;
};


class FunctionContainer
{
  // create a set of structs so we can specify function signatures
  // (ie, number of differentiaions needed) with an index
  template<unsigned int N>
  struct func_type{
    using type = EvaluationFunctor<autodiff_fvar<double,N>,autodiff_fvar<double,N>>*;
  };
  // special case for N=0
  template<>
  struct func_type<0>{
    using type = EvaluationFunctor<double,double>*;
  };

public: 
  // will become 'std::function<type(type)> type_func'
  EvaluationFunctor<double,double> *double_func;
  EvaluationFunctor<fvar1,fvar1>   *fvar1_func;
  EvaluationFunctor<fvar2,fvar2>   *fvar2_func;
  EvaluationFunctor<fvar3,fvar3>   *fvar3_func;
  EvaluationFunctor<fvar4,fvar4>   *fvar4_func;
  EvaluationFunctor<fvar5,fvar5>   *fvar5_func;
  EvaluationFunctor<fvar6,fvar6>   *fvar6_func;
  EvaluationFunctor<fvar7,fvar7>   *fvar7_func;

  // provide a method for accessing each member function uniformly
  // call as 'get_nth_func<N>()' for N auto derivatives
  template<int x> typename func_type<x>::type get_nth_func() 
  { 
    throw std::out_of_range("Template value must be in 0<=N<=7");
  };
  template<> typename func_type<0>::type get_nth_func<0>() { return double_func; };
  template<> typename func_type<1>::type get_nth_func<1>() { return fvar1_func;  };
  template<> typename func_type<2>::type get_nth_func<2>() { return fvar2_func;  };
  template<> typename func_type<3>::type get_nth_func<3>() { return fvar3_func;  };
  template<> typename func_type<4>::type get_nth_func<4>() { return fvar4_func;  };
  template<> typename func_type<5>::type get_nth_func<5>() { return fvar5_func;  };
  template<> typename func_type<6>::type get_nth_func<6>() { return fvar6_func;  };
  template<> typename func_type<7>::type get_nth_func<7>() { return fvar7_func;  };

  // some constructors
  FunctionContainer(){}
 
  FunctionContainer(EvaluationFunctor<double,double> *func,
      EvaluationFunctor<fvar1,fvar1> *func1, EvaluationFunctor<fvar2,fvar2> *func2,
      EvaluationFunctor<fvar3,fvar3> *func3, EvaluationFunctor<fvar4,fvar4> *func4,
      EvaluationFunctor<fvar5,fvar5> *func5, EvaluationFunctor<fvar6,fvar6> *func6,
      EvaluationFunctor<fvar7,fvar7> *func7) :
    double_func(func), fvar1_func(func1),
    fvar2_func(func2), fvar3_func(func3),
    fvar4_func(func4), fvar5_func(func5),
    fvar6_func(func6), fvar7_func(func7) {}
};

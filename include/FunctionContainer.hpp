/*
  A data structure used to pass mathematical functions to FunC tables.

  Notes:
  - needed for passing functions to tables. Tables using boost's
    automatic differentiation (such as Taylor, Hermite, Pade, and every NonUniformLUT)
    take advantage of the adVar[1-7] functions.
  - easy to use if your mathematical function is easy to template
  - copy and paste the following example code into a new file and use
    the command :%s/foo/new_name/g in vim or (TODO something else) in emacs
    to rename the example to your own function.
  - Boils down to a std::function wrapper if Boost's version is lower than 1.71
  - the readability of this class has been sacrificed in order to have get_nth_func
    return the ith order autodiff functions based on an index

  Example usage:
  #include FunctionContainer.hpp
  template <typename T>
  T foo(T x){ return x; }

  int main(){
    FunctionContainer<double,double> foo_container {FUNC_SET_F(foo,double)};
    // or if it's inconvenient to template your function
    FunctionContainer foo_container2 {foo<double>};
    // you just can't use any of the lookup tables that need
    // automatic differentiation (FunC will throw an exception if you do)
    return 0;
  }
*/
#pragma once
#include <stdexcept>
#include <functional>
#include <boost/math/differentiation/autodiff.hpp>
#include "config.hpp" // FUNC_USE_BOOST

#ifdef FUNC_USE_BOOST
/* Let the user quickly define their function container with a macro */
#define FUNC_SET_F_ONE_TYPE(F,TYPE)                                \
  F<TYPE>, F<adVar<TYPE,1>>, F<adVar<TYPE,2>>, F<adVar<TYPE,3>>, \
  F<adVar<TYPE,4>>, F<adVar<TYPE,5>>, F<adVar<TYPE,6>>, F<adVar<TYPE,7>>

#define FUNC_SET_F_TWO_TYPE(F,IN_TYPE,OUT_TYPE)                                  \
  F<IN_TYPE,OUT_TYPE>, F<adVar<IN_TYPE,1>,adVar<OUT_TYPE,1>>,                 \
  F<adVar<IN_TYPE,2>,adVar<OUT_TYPE,2>>, F<adVar<IN_TYPE,3>,adVar<OUT_TYPE,3>>, \
  F<adVar<IN_TYPE,4>,adVar<OUT_TYPE,4>>, F<adVar<IN_TYPE,5>,adVar<OUT_TYPE,5>>, \
  F<adVar<IN_TYPE,6>,adVar<OUT_TYPE,6>>, F<adVar<IN_TYPE,7>,adVar<OUT_TYPE,7>>

#define FUNC_GET_MACRO_FUNCTION_CONTAINER(_1,_2,_3,NAME,...) NAME
// Call with FUNC_SET_F(foo,template-type...)
#define FUNC_SET_F(...) FUNC_GET_MACRO_FUNCTION_CONTAINER(__VA_ARGS__, FUNC_SET_F_TWO_TYPE, FUNC_SET_F_ONE_TYPE, )(__VA_ARGS__)

namespace func {
// setup the automatically differentiable variable
using boost::math::differentiation::autodiff_fvar;
template <typename T, unsigned int N>
using adVar = autodiff_fvar<T,N>;

// create a set of structs so we can specify
// FunctionContainer::get_nth_func's return type with an index
template<typename IN_TYPE, typename OUT_TYPE, unsigned int N>
struct nth_differentiable{
  using type = std::function<adVar<OUT_TYPE,N>(adVar<IN_TYPE,N>)>;
};
/* special case for N=0 */
template<typename IN_TYPE, typename OUT_TYPE>
struct nth_differentiable<IN_TYPE,OUT_TYPE,0>{
  using type = std::function<OUT_TYPE(IN_TYPE)>;
};

#else
#define FUNC_SET_F(F,TYPE) F<TYPE>
#endif // FUNC_USE_BOOST

template<typename IN_TYPE, typename OUT_TYPE = IN_TYPE>
class FunctionContainer
{
#ifdef FUNC_USE_BOOST
  template<unsigned int N>
  using func_type = nth_differentiable<IN_TYPE,OUT_TYPE,N>;
  // 2 parter for providing a way to access each member function with a number.
  // overload func_type 9 different ways to get a function that seemingly does different
  // things based on its template value.
  template<unsigned int N> typename func_type<N>::type get_nth_func(func_type<N>)
  {
    throw std::out_of_range("Template value must be in 0<=N<=7");
  }
  typename func_type<0>::type get_nth_func(func_type<0>){ return standard_func;  };
  typename func_type<1>::type get_nth_func(func_type<1>){ return autodiff1_func; };
  typename func_type<2>::type get_nth_func(func_type<2>){ return autodiff2_func; };
  typename func_type<3>::type get_nth_func(func_type<3>){ return autodiff3_func; };
  typename func_type<4>::type get_nth_func(func_type<4>){ return autodiff4_func; };
  typename func_type<5>::type get_nth_func(func_type<5>){ return autodiff5_func; };
  typename func_type<6>::type get_nth_func(func_type<6>){ return autodiff6_func; };
  typename func_type<7>::type get_nth_func(func_type<7>){ return autodiff7_func; };

public:
  // call as func_container->template get_nth_func<N>() to get the member
  // function that is differentiated N times for each function call.
  template<unsigned int N>
  typename func_type<N>::type get_nth_func(){ return get_nth_func(func_type<N>()); }
#endif // FUNC_USE_BOOST
public:

  std::function<OUT_TYPE(IN_TYPE)> standard_func;
#ifdef FUNC_USE_BOOST
  std::function<adVar<OUT_TYPE,1>(adVar<IN_TYPE,1>)> autodiff1_func;
  std::function<adVar<OUT_TYPE,2>(adVar<IN_TYPE,2>)> autodiff2_func;
  std::function<adVar<OUT_TYPE,3>(adVar<IN_TYPE,3>)> autodiff3_func;
  std::function<adVar<OUT_TYPE,4>(adVar<IN_TYPE,4>)> autodiff4_func;
  std::function<adVar<OUT_TYPE,5>(adVar<IN_TYPE,5>)> autodiff5_func;
  std::function<adVar<OUT_TYPE,6>(adVar<IN_TYPE,6>)> autodiff6_func;
  std::function<adVar<OUT_TYPE,7>(adVar<IN_TYPE,7>)> autodiff7_func;
#endif // FUNC_USE_BOOST

  // some constructors
  FunctionContainer(){}

  FunctionContainer(std::function<OUT_TYPE(IN_TYPE)> func) :
    standard_func(func) {};

#ifdef FUNC_USE_BOOST
  FunctionContainer(std::function<OUT_TYPE(IN_TYPE)>   func,
      std::function<adVar<OUT_TYPE,1>(adVar<IN_TYPE,1>)> func1,
      std::function<adVar<OUT_TYPE,2>(adVar<IN_TYPE,2>)> func2,
      std::function<adVar<OUT_TYPE,3>(adVar<IN_TYPE,3>)> func3,
      std::function<adVar<OUT_TYPE,4>(adVar<IN_TYPE,4>)> func4,
      std::function<adVar<OUT_TYPE,5>(adVar<IN_TYPE,5>)> func5,
      std::function<adVar<OUT_TYPE,6>(adVar<IN_TYPE,6>)> func6,
      std::function<adVar<OUT_TYPE,7>(adVar<IN_TYPE,7>)> func7) :
    standard_func(func),   autodiff1_func(func1),
    autodiff2_func(func2), autodiff3_func(func3),
    autodiff4_func(func4), autodiff5_func(func5),
    autodiff6_func(func6), autodiff7_func(func7) {}
#endif // FUNC_USE_BOOST
};
} // namespace func

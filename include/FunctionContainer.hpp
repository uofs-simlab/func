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
#include "config.hpp" // FUNC_USE_BOOST

#ifdef FUNC_USE_BOOST
#include <boost/math/differentiation/autodiff.hpp>
// two helper macros to make FUNC_SET_F appear like a variadic macro
#define FUNC_SET_F_ONE_TYPE(F,TYPE)               \
  F<TYPE>, F<func::adVar<TYPE,1>>,                \
  F<func::adVar<TYPE,2>>, F<func::adVar<TYPE,3>>, \
  F<func::adVar<TYPE,4>>, F<func::adVar<TYPE,5>>, \
  F<func::adVar<TYPE,6>>, F<func::adVar<TYPE,7>>

#define FUNC_SET_F_TWO_TYPE(F,TIN,TOUT)                                                 \
  F<TIN,TOUT>, F<func::adVar<TIN,1>,func::adVar<TOUT,1>>,                               \
  F<func::adVar<TIN,2>,func::adVar<TOUT,2>>, F<func::adVar<TIN,3>,func::adVar<TOUT,3>>, \
  F<func::adVar<TIN,4>,func::adVar<TOUT,4>>, F<func::adVar<TIN,5>,func::adVar<TOUT,5>>, \
  F<func::adVar<TIN,6>,func::adVar<TOUT,6>>, F<func::adVar<TIN,7>,func::adVar<TOUT,7>>

namespace func {
// convenient typename for Boost's forward mode autodiff variable
using boost::math::differentiation::autodiff_fvar;
template <typename T, unsigned int N>
using adVar = autodiff_fvar<T,N>;

// create a set of structs so we can specify
// FunctionContainer::get_nth_func's return type with an index
template<typename TIN, typename TOUT, unsigned int N>
struct nth_differentiable{
  using type = std::function<adVar<TOUT,N>(adVar<TIN,N>)>;
};
/* special case for N=0 */
template<typename TIN, typename TOUT>
struct nth_differentiable<TIN,TOUT,0>{
  using type = std::function<TOUT(TIN)>;
};

#else
#define FUNC_SET_F_ONE_TYPE(F,TYPE)     F<TYPE>
#define FUNC_SET_F_TWO_TYPE(F,TIN,TOUT) F<TIN,TOUT>

namespace func {
#endif // FUNC_USE_BOOST

#define FUNC_GET_MACRO_FUNCTION_CONTAINER(_1,_2,_3,NAME,...) NAME
// Call with FUNC_SET_F(foo,template-type...)
#define FUNC_SET_F(...) FUNC_GET_MACRO_FUNCTION_CONTAINER(__VA_ARGS__, FUNC_SET_F_TWO_TYPE, FUNC_SET_F_ONE_TYPE, )(__VA_ARGS__)


template<typename TIN, typename TOUT = TIN>
class FunctionContainer
{
#ifdef FUNC_USE_BOOST
  template<unsigned int N>
  using fun_type = nth_differentiable<TIN,TOUT,N>;
  // 2 parter for providing a way to access each member function with a number.
  // overload fun_type 9 different ways to get a function that seemingly does different
  // things based on its template value.
  template<unsigned int N> typename fun_type<N>::type get_nth_func(fun_type<N>) const
  {
    throw std::out_of_range("Template value must be in 0<=N<=7");
  }
  typename fun_type<0>::type get_nth_func(fun_type<0>) const { return standard_fun;  };
  typename fun_type<1>::type get_nth_func(fun_type<1>) const { return autodiff1_fun; };
  typename fun_type<2>::type get_nth_func(fun_type<2>) const { return autodiff2_fun; };
  typename fun_type<3>::type get_nth_func(fun_type<3>) const { return autodiff3_fun; };
  typename fun_type<4>::type get_nth_func(fun_type<4>) const { return autodiff4_fun; };
  typename fun_type<5>::type get_nth_func(fun_type<5>) const { return autodiff5_fun; };
  typename fun_type<6>::type get_nth_func(fun_type<6>) const { return autodiff6_fun; };
  typename fun_type<7>::type get_nth_func(fun_type<7>) const { return autodiff7_fun; };

public:
  // call as func_container->template get_nth_func<N>() to get the member
  // function that is differentiated N times for each function call.
  template<unsigned int N>
  typename fun_type<N>::type get_nth_func() const { return get_nth_func(fun_type<N>()); }
#endif // FUNC_USE_BOOST
public:

  std::function<TOUT(TIN)> standard_fun;
#ifdef FUNC_USE_BOOST
  std::function<adVar<TOUT,1>(adVar<TIN,1>)> autodiff1_fun;
  std::function<adVar<TOUT,2>(adVar<TIN,2>)> autodiff2_fun;
  std::function<adVar<TOUT,3>(adVar<TIN,3>)> autodiff3_fun;
  std::function<adVar<TOUT,4>(adVar<TIN,4>)> autodiff4_fun;
  std::function<adVar<TOUT,5>(adVar<TIN,5>)> autodiff5_fun;
  std::function<adVar<TOUT,6>(adVar<TIN,6>)> autodiff6_fun;
  std::function<adVar<TOUT,7>(adVar<TIN,7>)> autodiff7_fun;
#endif // FUNC_USE_BOOST

  // some constructors
  FunctionContainer(){}

  FunctionContainer(std::function<TOUT(TIN)> fun) :
    standard_fun(fun) {};

#ifdef FUNC_USE_BOOST
  FunctionContainer(std::function<TOUT(TIN)>   fun,
      std::function<adVar<TOUT,1>(adVar<TIN,1>)> fun1,
      std::function<adVar<TOUT,2>(adVar<TIN,2>)> fun2,
      std::function<adVar<TOUT,3>(adVar<TIN,3>)> fun3,
      std::function<adVar<TOUT,4>(adVar<TIN,4>)> fun4,
      std::function<adVar<TOUT,5>(adVar<TIN,5>)> fun5,
      std::function<adVar<TOUT,6>(adVar<TIN,6>)> fun6,
      std::function<adVar<TOUT,7>(adVar<TIN,7>)> fun7) :
    standard_fun(fun),   autodiff1_fun(fun1),
    autodiff2_fun(fun2), autodiff3_fun(fun3),
    autodiff4_fun(fun4), autodiff5_fun(fun5),
    autodiff6_fun(fun6), autodiff7_fun(fun7) {}
#endif // FUNC_USE_BOOST
};

} // namespace func

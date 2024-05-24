/** \file FunctionContainer.hpp
 * \brief Define a FunctionContainer with the (necessary) convenience function
 * `get_nth_func()`, the "variadic" macro `#FUNC_SET_F`.
 * \note Boost's autodiff.cpp is only included if CMake found a new enough version of Boost */
#pragma once
#include <stdexcept>
#include <functional>
#include "config.hpp" // FUNC_USE_BOOST

#ifdef FUNC_USE_BOOST
#include <boost/math/differentiation/autodiff.hpp>

// two helper macros to make FUNC_SET_F work like a variadic macro
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
/** convenient shorthand for Boost's forward mode autodiff variable */
using boost::math::differentiation::autodiff_fvar;
template <typename T, unsigned int N>
using adVar = autodiff_fvar<T,N>;

/** \brief These structs provide an indexed typedef for Boost's autodiff_fvar type.
 *
 * The advantage of using `nth_differentiable` instead of a normal typedef is
 * that we can define the function overloads for FunctionContainer::get_nth_func
 * (the return type must be specified with an index). 
 * \note Does not exist in FunC if Boost's autodiff is not available */
template<typename TIN, typename TOUT, unsigned int N>
struct nth_differentiable{
  using type = std::function<adVar<TOUT,N>(adVar<TIN,N>)>;
};
/** \brief Base case for `nth_differentiable` when N=0 */
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

/** \brief Macro to list out FunctionContainer constructor arguments.
 * Call as `FUNC_SET_F(foo,template-type...)` where
 * - foo is your function name, and
 * - template-type is as list of the 1 or 2 types it maps between. */
#define FUNC_SET_F(...) FUNC_GET_MACRO_FUNCTION_CONTAINER(__VA_ARGS__, FUNC_SET_F_TWO_TYPE, FUNC_SET_F_ONE_TYPE, )(__VA_ARGS__)


/**
  \brief Wrapper for `std::function<TOUT(TIN)>` and some optional `std::functions`
  of Boost's automatic differentiation type

  Used to pass mathematical functions to FunC's LUTs.

  Notes:
  - Only the LUTs that need derivatives (Taylor, Hermite, Pade, and every
    NonUniformLUT) need Boost's adVar[1-7] functions.
  - The automatic differentiation requires the mathematical function is templated
    on some abstract type
  - Autodiff was introduced in Boost 1.71. This class reduces to a simple std::function
    wrapper if Boost is not available or is too old.
  - Most of the machinery is necessary to use `get_nth_func<N>` which returns
    the ith order autodiff functions based on an index

  - copy and paste the following example code into a new file and
    rename the example to your own function with
    1. :%s/foo/new_name/g in Vim or
    2. (TODO whatever it is) in emacs.

  Example:
  \code{.cpp}
  #include FunctionContainer.hpp
  template <typename T>
  T foo(T x){ return x; }

  #define TYPE double

  int main(){
    FunctionContainer<TYPE> foo_container {FUNC_SET_F(foo,TYPE)};

    // Use this version if it's inconvenient to template your function:
    FunctionContainer<TYPE> foo_container2 {foo<TYPE>};
    // The only downside is then you can't generate LUTs that need
    // automatic differentiation (runtime exception if you try)
    return 0;
  }
  \endcode */
template<typename TIN, typename TOUT = TIN>
class FunctionContainer
{
#ifdef FUNC_USE_BOOST
  template<unsigned int N>
  using fun_type = nth_differentiable<TIN,TOUT,N>;

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
  /** \brief return the Boost autodiff function that automatically
   * differentiates the user's function N times.
   *
   * \note call as func_container->template get_nth_func<N>() to get the member
   * \note Only supports up to 7 differentiations.
   * \note Implemented with 9 different overloads of fun_type<N> */
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

  FunctionContainer(){}
  FunctionContainer(std::function<TOUT(TIN)> fun) :
    standard_fun(fun) {};

#ifdef FUNC_USE_BOOST
  /** \brief Users should not use this function directly. See the example usage
   * for FunctionContainer and the macro `FUNC_SET_F(...)` */
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

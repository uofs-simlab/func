/*
  A data structure used to pass mathematical functions to FunC tables.

  Notes:
  - needed for passing functions to tables. Tables using boost's 
    automatic differentiation (such as Taylor, Hermite, and Pade) 
    take advantage of the fvar[1-9] functions.
  - easy (though obtuse) to use if mathematical functions are templated
  - copy and paste the following example code into a new file and use
    the command :%s/foo/new_name/g in vim or (TODO something else) in emacs 
    to rename the example to your own function.
  - some of this class is illegible

  Example usage:
  #include FunctionContainer.hpp
  template <typename T>
  T foo(T x){ return x; }

  int main(){
    FunctionContainer foo_container{SET_F(foo)};
    return 0;
  }
*/
#pragma once
#include <stdexcept>
#include <functional>
#include <boost/math/differentiation/autodiff.hpp>

// some typedefs
using boost::math::differentiation::autodiff_fvar;
// will these typedefs still work if we template this class?
typedef autodiff_fvar<double,1> fvar1;
typedef autodiff_fvar<double,2> fvar2;
typedef autodiff_fvar<double,3> fvar3;
typedef autodiff_fvar<double,4> fvar4;
typedef autodiff_fvar<double,5> fvar5;
typedef autodiff_fvar<double,6> fvar6;
typedef autodiff_fvar<double,7> fvar7;

/* Used by each table type to check if the required function
 * type has been provided.
 * TODO decide whether or not to make this null op for -DNDEBUG */
#define __IS_NULL(VAR) \
  if(VAR == NULL)      \
    throw std::invalid_argument(#VAR " is NULL")

#define SET_F(F) \
  F<double>, F<fvar1>,F<fvar2>, F<fvar3>, F<fvar4>, F<fvar5>, F<fvar6>, F<fvar7>

// create a set of structs so we can specify 
// FunctionContainer::get_nth_func's return type with an index
template<unsigned int N>
struct func_type{
  using type = std::function<autodiff_fvar<double,N>(autodiff_fvar<double,N>)>;
};
// special case for N=0
template<>
struct func_type<0>{
  using type = std::function<double(double)>;
};

class FunctionContainer
{
private:
  // 2 parter for providing a way to access each member function with a number.
  // overload func_type 9 different ways to get a function that seemingly does different
  // things based on its template value.
  template<unsigned int N> typename func_type<N>::type get_nth_func(func_type<N>)
  { 
    throw std::out_of_range("Template value must be in 0<=N<=7");
  };
  typename func_type<0>::type get_nth_func(func_type<0>) { return double_func; };
  typename func_type<1>::type get_nth_func(func_type<1>) { return fvar1_func;  };
  typename func_type<2>::type get_nth_func(func_type<2>) { return fvar2_func;  };
  typename func_type<3>::type get_nth_func(func_type<3>) { return fvar3_func;  };
  typename func_type<4>::type get_nth_func(func_type<4>) { return fvar4_func;  };
  typename func_type<5>::type get_nth_func(func_type<5>) { return fvar5_func;  };
  typename func_type<6>::type get_nth_func(func_type<6>) { return fvar6_func;  };
  typename func_type<7>::type get_nth_func(func_type<7>) { return fvar7_func;  };

public: 
  // call as get_nth_func<N>() to get the member function that is differentiated N 
  // times for each function call
  template<unsigned int N>
  typename func_type<N>::type get_nth_func(){ return get_nth_func(func_type<N>()); }

  std::function<double(double)> double_func;
  std::function<fvar1(fvar1)>   fvar1_func;
  std::function<fvar2(fvar2)>   fvar2_func;
  std::function<fvar3(fvar3)>   fvar3_func;
  std::function<fvar4(fvar4)>   fvar4_func;
  std::function<fvar5(fvar5)>   fvar5_func;
  std::function<fvar6(fvar6)>   fvar6_func;
  std::function<fvar7(fvar7)>   fvar7_func;
    
  // some constructors
  FunctionContainer(){}
 
  FunctionContainer(std::function<double(double)> func,
      std::function<fvar1(fvar1)> func1, std::function<fvar2(fvar2)> func2,
      std::function<fvar3(fvar3)> func3, std::function<fvar4(fvar4)> func4,
      std::function<fvar5(fvar5)> func5, std::function<fvar6(fvar6)> func6,
      std::function<fvar7(fvar7)> func7) :
    double_func(func), fvar1_func(func1),
    fvar2_func(func2), fvar3_func(func3),
    fvar4_func(func4), fvar5_func(func5),
    fvar6_func(func6), fvar7_func(func7) {}
};

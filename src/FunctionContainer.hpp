/*
  A data structure used to pass mathematical functions to FunC tables.

  - needed for passing functions to tables. Tables using boost's 
    automatic differentiation (such as Taylor, Hermite, and Pade) 
    take advantage of the fvar[1-9] functions.
  - easy (though obtuse) to use if mathematical functions are templated
  - copy and paste the following example code into a new file and use
    the command :%s/foo/new_name/g in vim or (TODO something else) in emacs 
    to rename the example to your own function.

  Example usage:
  #include FunctionContainer.hpp
  template <typename T>
  T foo(T x){ return x; }

  int main(){
    FunctionContainer foo_container{foo<double>, foo<fvar1>,
        foo<fvar2>, foo<fvar3>, foo<fvar4>, 
        foo<fvar5>, foo<fvar6>, foo<fvar7>};
    // or just
    // FunctionContainer foo_container;
    // foo_container.double_func = foo<double>;
    // if you're not using any Taylor/Pade/Hermite tables
    // (You'll get an exception if something needed is not set)
    return 0;
  }
*/
#pragma once
#include <stdexcept>
#include <functional>
#include <boost/math/differentiation/autodiff.hpp>

// some typedefs
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
#define __IS_NULL(VAR) \
    if(VAR == NULL)     \
      throw std::invalid_argument(#VAR " is NULL")

class FunctionContainer
{
  // create a set of structs so we can specify function signatures
  // (ie, number of differentiaions needed) with an index
  template<unsigned int N>
  struct func_type{
    using type = std::function<autodiff_fvar<double,N>(autodiff_fvar<double,N>)>;
  };
  // special case for N=0
  template<>
  struct func_type<0>{
    using type = std::function<double(double)>;
  };

public: 
  std::function<double(double)> double_func;
  std::function<fvar1(fvar1)>   fvar1_func;
  std::function<fvar2(fvar2)>   fvar2_func;
  std::function<fvar3(fvar3)>   fvar3_func;
  std::function<fvar4(fvar4)>   fvar4_func;
  std::function<fvar5(fvar5)>   fvar5_func;
  std::function<fvar6(fvar6)>   fvar6_func;
  std::function<fvar7(fvar7)>   fvar7_func;
  
  // provide a method for accessing each member function uniformly.
  // Call as 'get_nth_func<N>()' for N derivatives to be calculated for 
  // each function call
  template<unsigned int N> typename func_type<N>::type get_nth_func() 
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

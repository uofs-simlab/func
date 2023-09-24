
#pragma once
#ifndef FUNC_NO_CXX17
#include "LookupTable.hpp"
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

namespace func {

template <unsigned int N, typename TERR, class LUT, typename F>
TERR metric(LUT L, F f){
  TERR max_err = 0;
  if constexpr (N == 0) {
    using std::abs;
    auto fval = f();
    max_err = -abs(L - fval)/(static_cast<TERR>(1.0) + abs(fval));
  } else {
    using namespace boost::math::tools;
    /* get number of binary bits in mantissa */
    int bits = std::numeric_limits<TERR>::digits/2; // effective maximum for brent_find_minima
    boost::uintmax_t max_it = 52;
    #pragma omp parallel for
    for(unsigned ii=0; ii<L.num_subintervals(); ii++){
      auto intEndPoints = L.bounds_of_subinterval(ii);
      TERR x = static_cast<TERR>(boost::math::float_next(intEndPoints.first));
      TERR xtop = static_cast<TERR>(boost::math::float_prior(intEndPoints.second));
      if(intEndPoints.second > L.max_arg())
        xtop = static_cast<TERR>(L.max_arg());

      std::pair<TERR,TERR> r = brent_find_minima(
          [&L,&f](const TERR& z){ 
            auto z2 = static_cast<typename decltype(L)::input_type>(z);
            return metric<TERR,N-1>(L(z2), [&z2,&f](auto... args){ return f(z2, args...); }); },
          x, xtop, bits, max_it
      );

      #pragma omp critical
      {
        max_err = std::min(max_err, r.second);
      }
    }
  }
  return max_err;
}

/* convenience functions for defining LUTs of LUTs */
template <unsigned int N, typename TIN, typename TOUT, template <typename...> class classname>
struct curriedLUT {
  using type = classname<TIN, typename curriedLUT<N-1,TIN,TOUT,classname>::type>;
};

template<typename TIN, typename TOUT, template <typename...> class classname>
struct curriedLUT<0,TIN,TOUT,classname> {
  using type = classname<TIN,TOUT>;
};

/* TODO:
 * - Is not compatible with function derivatives (needed for nonuniform partition & Taylor tables). Is that impossible to support anyways???
 * - classname should be variadic! Is this possible?????
 * Call as func::ndimLUT<ndim,type1,type2,luttype>(f, params) */
template<unsigned int N, typename TIN, typename TOUT, template <typename...> class classname, class F, typename... TIN2>
constexpr typename curriedLUT<N-1,TIN,TOUT,classname>::type ndimLUT(F f, const std::vector<LookupTableParameters<TIN>>& params, TIN2... other) {
  auto H = sizeof...(other);
  if constexpr (N == 1) {
    return typename curriedLUT<0  ,TIN,TOUT,classname>::type({[&f,other...](TIN x){ return f(other..., x); }}, params[H]);
  } else {
    return typename curriedLUT<N-1,TIN,TOUT,classname>::type({[&f,&params,other...](TIN x){ return ndimLUT<N-1,TIN,TOUT,classname>(f, params, other..., x); }}, params[H]);
  }
}

}

#endif

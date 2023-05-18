/*
  Accepts an array of polynomial coefficients.
  We only mandate that a TransferFunction is a monotonically increasing
  cubic polynomial, and the table hash must be baked into those coefficients.
  For example, it initializes itself the identity transfer function which
  is the linear polynomial with coefs {-minArg/stepSize,1/stepSize,0,0}.

  If this class is given a function container then it will use
  Boost's automatic differentiation to generate a cubic monotone Hermite
  interpolating polynomial which will be a decent TransferFunction.

  Given a function f and its domain [a,b], build a pair of functions
  g: [a,b] -> [a,b] and it's inverse g^{-1}.
  g is formally defined as
    g(x) = a + \frac{b-a}{c}\int_a^x\frac{dt}{\sqrt{1+[f'(t)]^2}}.
  where c = \int_a^b\frac{dt}{\sqrt{1+[f'(t)]^2}}.
  g_inv has to be fast so we approximate is as a monotone Hermite
  cubic polynomial and then rebuilt g as g_inv^{-1} using
  a Newton's method/ bisection mix.

Notes:
  - We seem to get better grids if f' is largest near the endpoints [a,b].
    If f' is largest near the middle of an interval (eg e^{-x^2} when a<-3 and b>3)
    then that information is largely ignored.
  - We need Boost version 1.71.0 or newer to generate a candidate transfer function.
    Boost is not required if TransferFunctionSinh is given pre-generated coefficients.
 */
#pragma once
#include "config.hpp" // FUNC_USE_BOOST
#include <cmath> // sqrt
#include <limits> // std::numeric_limits<T>::digits
#include <array> // std::array
#include <utility> // std::pair
#include <type_traits>

#include "FunctionContainer.hpp"

#ifdef FUNC_USE_BOOST
#define BOOST_MATH_GAUSS_NO_COMPUTE_ON_DEMAND
#include <boost/math/quadrature/gauss_kronrod.hpp> // gauss_kronrod::integrate
#include <boost/math/tools/roots.hpp> // newton_raphson_iterate
#endif

namespace func {

template <typename TIN>
class TransferFunctionSinh
{
  // aligned array of polynomial coefs approximating g_inv
  __attribute__((aligned)) std::array<TIN,4> m_inv_coefs = {{0,0,0,0}};
  /* This min, max must be the same as the corresponding table's min and
  max resp. (though note that the table max is not necessarily equal to
  the function's max arg) */
  //TIN minArg, tableMaxArg;
  //TIN stepSize;
public:
  /* Set m_inv_coefs equal to a vector that is (presumably) either the identity or came from a json file */
  TransferFunctionSinh(const std::array<TIN,4>& inv_coefs) { m_inv_coefs = inv_coefs; }

  TransferFunctionSinh() = default;
  /* initialize the identity transfer function */
  //TransferFunctionSinh(TIN minArg, TIN tableMaxArg, TIN stepSize) :
  //  TransferFunctionSinh<TIN>(minArg, tableMaxArg, stepSize, {-minArg/stepSize,1/stepSize,0,0}) {}

  /* Build the coefficients in g_inv 
   * TODO do SFINAE because TOUT need only be defined for sqrt() (ie TOUT cannot be LookupTable<...>).
   * So, it doesn't necessarily need to be accepted by std::is_arithmetic */
  template<typename TOUT>
           //typename std::enable_if<std::is_arithmetic<TOUT>::value, bool>::type = true>
  TransferFunctionSinh(const FunctionContainer<TIN,TOUT>& fc, TIN minArg, TIN tableMaxArg, TIN stepSize) :
  {
#ifndef FUNC_USE_BOOST
    /* cause a compile time error precisely when this constructor is called without Boost: */
    static_assert(sizeof(TIN) != sizeof(TIN), "Cannot generate a nonuniform grid without Boost verion 1.71.0 or higher");
#else
    using boost::math::quadrature::gauss_kronrod;
    using boost::math::differentiation::make_fvar;

    if(fc->autodiff1_func == nullptr)
      throw std::invalid_argument("Error in func::TransferFunction. 1st derivative of function is needed to generate nonuniform grids but is null");

    // build a function to return the first derivative of f
    std::function<TOUT(TIN)> f_prime = [&fc](TIN x) -> TOUT {
      return (fc->autodiff1_func)(make_fvar<TIN,1>(x)).derivative(1);
    };

    // build the integrand. Notice that it is strictly positive
    std::function<TIN(TIN)> integrand = [&f_prime](TIN x) -> TIN {
      return static_cast<TIN>(1.0)/static_cast<TIN>(sqrt(static_cast<TOUT>(1.0) + f_prime(x)*f_prime(x)));
    };

    TIN a = minArg;
    TIN b = tableMaxArg;

    // integrate over [a,b] using adaptive quadrature & a default tol of sqrt(epsilon).
    TIN c = gauss_kronrod<TIN, 15>::integrate(integrand, a, b);

    // rescale the integrand to map [a,b] -> [a,b]
    std::function<TIN(TIN)> g_prime = [&integrand,a,b,c](TIN x) -> TIN {
      return (b-a) * integrand(x)/c;
    };

    // now g(a)=a and g(b)=b so
    // g_inv_prime(a) = 1/g_prime(a) and g_inv_prime(b) = 1/g_prime(b)
    TIN m0 = 1/g_prime(a);
    TIN m1 = 1/g_prime(b);

    /* ensure monotonicity of the Hermite interpolating polynomial.
     * TODO this doesn't capture every condition of cubic monotone polynomials
     * Note that m0,m1 >= 0 */
    m0 = (m0 > 3) ? 3 : m0;
    m1 = (m1 > 3) ? 3 : m1;

    // Compute the Hermite interpolating polynomial p for g_inv satisfying
    // p(a)=a, p(b)=b, p'(a)=m0, p'(b)=m1
    // (the symbolic expression was computed with Matlab)
    m_inv_coefs[0] = (a*b*(a + b - a*m1 - b*m0))/(a - b)/(a - b);
    m_inv_coefs[1] = (a*a*m1 - 6*a*b + b*b*m0 + 2*a*b*m0 + 2*a*b*m1)/(a - b)/(a - b);
    m_inv_coefs[2] = -(a*m0 - 3*b - 3*a + 2*a*m1 + 2*b*m0 + b*m1)/(a - b)/(a - b);
    m_inv_coefs[3] = (m0 + m1 - 2)/(a - b)/(a - b);

    /* build the version of g_inv we'll actually use by encoding the
       underlying table's hash into the transfer function eval.
       This way, the only cost of using a transfer function is to lookup 4 stack allocated numbers */
    m_inv_coefs[0] = m_inv_coefs[0] - minArg;
    for(unsigned int i=0; i<4; i++)
      m_inv_coefs[i] = m_inv_coefs[i] / stepSize;
#endif
  }

  TIN g_inv(TIN x)
  {
    // Horner's method
    TIN sum = 0;
    for (int k=3; k>0; k--)
      sum = x*(m_inv_coefs[k] + sum);
    return sum + m_inv_coefs[0];
  }

  TIN g_inv_prime(TIN x)
  {
    // Horner's method
    TIN sum = 0;
    for (int k=3; k>1; k--)
      sum = x*(k*m_inv_coefs[k] + sum);
    return sum + m_inv_coefs[1];
  }

  /* Use newton-raphson_iteration on g_inv. Recall that each coef in g was divided by h
   * and we subtracted by minArg */
  TIN g(TIN x)
  {
    // This will have at least 0.9*std::numeric_limits<TIN>::digits digits of accuracy
    boost::uintmax_t maxit = 55;
    return boost::math::tools::newton_raphson_iterate(
        [this,x](const TIN z){ return std::make_tuple(g_inv(z) - (x-minArg)/stepSize, g_inv_prime(z));},
        x,minArg, tableMaxArg, 0.9*std::numeric_limits<TIN>::digits, maxit);
  }

  // public access to private vars
  void print_details(std::ostream& out)
  {
    out << "degree 3 monotone Hermite interpolation with polynomial: \n";
    out << std::to_string(m_inv_coefs[3]) << "x^3 + " << std::to_string(m_inv_coefs[2]) << "x^2 + " <<
      std::to_string(m_inv_coefs[1]) << "x + " << std::to_string(m_inv_coefs[0]) << std::endl;
  }

  std::array<TIN,4> get_coefs() const { return m_inv_coefs; }
  std::pair<TIN,TIN> arg_bounds_of_interval(){ return std::make_pair(minArg, tableMaxArg); }
};
}

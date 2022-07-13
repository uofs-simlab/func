/*
  Accepts an array of polynomial coefficients.
  We only mandate that a TransferFunction is a monotonically increasing
  cubic polynomial, and the table hash must be baked into those coefficients.
  For example, it initializes itself the identity transfer function which
  is the linear polynomial with coefs {-m_minArg/m_stepSize,1/m_stepSize,0,0}.

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

#include "FunctionContainer.hpp"

#ifdef FUNC_USE_BOOST
#define BOOST_MATH_GAUSS_NO_COMPUTE_ON_DEMAND
#include <boost/math/quadrature/gauss_kronrod.hpp> // gauss_kronrod::integrate
#include <boost/math/tools/roots.hpp> // newton_raphson_iterate
#endif

namespace func {

template <typename IN_TYPE>
class TransferFunctionSinh
{
  /* This min, max must be the same as the corresponding table's min and
  max resp. (though note that the table max is not necessarily equal to
  the function's max arg) */
  IN_TYPE m_minArg, m_tableMaxArg;
  IN_TYPE m_stepSize;

  // aligned array of polynomial coefs approximating g_inv
  __attribute__((aligned)) std::array<IN_TYPE,4> m_inv_coefs = {{0,0,0,0}};

public:
  /* Set m_inv_coefs equal to a vector that is (presumably) either the identity or came from a json file */
  TransferFunctionSinh(IN_TYPE minArg, IN_TYPE tableMaxArg, IN_TYPE stepSize, std::array<IN_TYPE,4> inv_coefs) :
    m_minArg(minArg), m_tableMaxArg(tableMaxArg), m_stepSize(stepSize) { m_inv_coefs = inv_coefs; }

  TransferFunctionSinh() = default;
  /* initialize the identity transfer function */
  //TransferFunctionSinh(IN_TYPE minArg, IN_TYPE tableMaxArg, IN_TYPE stepSize) :
  //  TransferFunctionSinh<IN_TYPE>(minArg, tableMaxArg, stepSize, {-minArg/stepSize,1/stepSize,0,0}) {}

  /* Build the coefficients in g_inv */
  template<typename OUT_TYPE>
  TransferFunctionSinh(FunctionContainer<IN_TYPE,OUT_TYPE> *fc,
      IN_TYPE minArg, IN_TYPE tableMaxArg, IN_TYPE stepSize) :
    m_minArg(minArg), m_tableMaxArg(tableMaxArg), m_stepSize(stepSize)
  {
#ifndef FUNC_USE_BOOST
    // Template code is only compiled if the template is instantiated
    // so this will cause a compile time error only when this
    // constructor is called without Boost available:
    static_assert(sizeof(IN_TYPE) != sizeof(IN_TYPE), "Cannot generate a nonuniform grid without Boost verion 1.71.0 or higher");
#else
    using boost::math::quadrature::gauss_kronrod;
    using boost::math::differentiation::make_fvar;

    if(fc->autodiff1_func == nullptr)
      throw std::invalid_argument("Error in func::TransferFunction. 1st derivative of function is needed to generate nonuniform grids but is null");

    // build a function to return the first derivative of f
    std::function<OUT_TYPE(IN_TYPE)> f_prime = [fc](IN_TYPE x) -> OUT_TYPE {
      return (fc->autodiff1_func)(make_fvar<IN_TYPE,1>(x)).derivative(1);
    };

    // build the integrand. Notice that it is strictly positive
    std::function<IN_TYPE(IN_TYPE)> integrand = [f_prime](IN_TYPE x) -> IN_TYPE {
          return 1/((IN_TYPE) sqrt(1 + f_prime(x)*f_prime(x)));
        };

    IN_TYPE a = m_minArg;
    IN_TYPE b = m_tableMaxArg;

    // integrate over [a,b] using adaptive quadrature & a default tol of sqrt(epsilon).
    IN_TYPE c = gauss_kronrod<IN_TYPE, 15>::integrate(integrand, a, b);

    // rescale the integrand to map [a,b] -> [a,b]
    std::function<IN_TYPE(IN_TYPE)> g_prime = [integrand,a,b,c](IN_TYPE x) -> IN_TYPE
    {
      return (b-a) * integrand(x)/c;
    };

    // now g(a)=a and g(b)=b so
    // g_inv_prime(a) = 1/g_prime(a) and g_inv_prime(b) = 1/g_prime(b)
    IN_TYPE m0 = 1/g_prime(a);
    IN_TYPE m1 = 1/g_prime(b);

    // ensure monotonicity of the Hermite interpolating polynomial.
    // We already know m0,m1 >= 0
    m0 = (m0 > 3) ? 3 : m0;
    m1 = (m1 > 3) ? 3 : m1;

    // Compute the Hermite interpolating polynomial p for g_inv satisfying
    // p(a)=a, p(b)=b, p'(a)=m0, p'(b)=m1
    // (the symbolic expression was computed with Matlab)
    m_inv_coefs[0] = (a*b*(a + b - a*m1 - b*m0))/(a - b)/(a - b);
    m_inv_coefs[1] = (a*a*m1 - 6*a*b + b*b*m0 + 2*a*b*m0 + 2*a*b*m1)/(a - b)/(a - b);
    m_inv_coefs[2] = -(a*m0 - 3*b - 3*a + 2*a*m1 + 2*b*m0 + b*m1)/(a - b)/(a - b);
    m_inv_coefs[3] = (m0 + m1 - 2)/(a - b)/(a - b);

    /* build the real version of g_inv by encoding the
       underlying table's hash into the transfer function eval. This
       will obfuscate our code, but now the only price of using
       nonuniform tables is an indirection */
    m_inv_coefs[0] = m_inv_coefs[0] - m_minArg;
    for(unsigned int i=0; i<4; i++)
      m_inv_coefs[i] = m_inv_coefs[i] / stepSize;
#endif
  }

  IN_TYPE g_inv(IN_TYPE x)
  {
    // Horner's method
    IN_TYPE sum = 0;
    for (int k=3; k>0; k--)
      sum = x*(m_inv_coefs[k] + sum);
    return sum + m_inv_coefs[0];
  }

  IN_TYPE g_inv_prime(IN_TYPE x)
  {
    // Horner's method
    IN_TYPE sum = 0;
    for (int k=3; k>1; k--)
      sum = x*(k*m_inv_coefs[k] + sum);
    return sum + m_inv_coefs[1];
  }

  /* Use newton-raphson_iteration on g_inv. Recall that each coef in g was divided by h
   * and we subtracted by m_minArg */
  IN_TYPE g(IN_TYPE x)
  {
    // This will have at least 0.9*std::numeric_limits<IN_TYPE>::digits digits of accuracy
    boost::uintmax_t maxit = 55;
    return boost::math::tools::newton_raphson_iterate(
        [this,x](const IN_TYPE z){ return std::make_tuple(g_inv(z) - (x-m_minArg)/m_stepSize, g_inv_prime(z));},
        x,m_minArg, m_tableMaxArg, 0.9*std::numeric_limits<IN_TYPE>::digits, maxit);
  }

  // public access to private vars
  void print_details(std::ostream& out)
  {
    out << "degree 3 monotone Hermite interpolation with polynomial: \n";
    out << std::to_string(m_inv_coefs[3]) << "x^3 + " << std::to_string(m_inv_coefs[2]) << "x^2 + " <<
      std::to_string(m_inv_coefs[1]) << "x + " << std::to_string(m_inv_coefs[0]) << std::endl;
  }

  std::array<IN_TYPE,4> get_coefs() const { return m_inv_coefs; }

  std::pair<IN_TYPE,IN_TYPE> arg_bounds_of_interval(){ return std::make_pair(m_minArg, m_tableMaxArg); }
};
}

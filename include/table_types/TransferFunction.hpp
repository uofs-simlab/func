/*
  A TransferFunction transforms a uniformly spaced partition of $[a,b]$ into a nonuniform partition of $[a,b]$.
  For efficiency, we require a Transfer function is simply an increasing cubic polynomial such that p(a)=0, p(b)=b/stepSize.
  To help compete against uniform lookup tables, part of the operator() must be baked into those coefficients (hence the appearance of stepsize).
  Can be constructed with a std::array<TIN,4> of polynomial coefficients.
  For example, the identity transfer function which is the linear polynomial with coefs {-m_minArg/m_stepSize,1/m_stepSize,0,0}.

  When given a FunctionContainer with a defined first derivative then this will construct a valid TransferFunction.
  Let f be a function defined over [a,b]. Then, construct S: [a,b] -> [a,b] as
    S(x) = a + \frac{b-a}{c}\int_a^x\frac{1}{\sqrt{1+[f'(t)]^2}} dt.
  where c = \int_a^b\frac{1}{\sqrt{1+[f'(t)]^2}} dt. This way, S(a)=a, S(b)=b.
  Computing S^{-1} has to be fast so we approximate is as a monotone Hermite
  cubic polynomial and then rebuild S as (S^{-1})^{-1} using
  Boost's newton_raphson_iterate.

Notes:
  - We need Boost version 1.71.0 or newer to generate a candidate transfer function.
    Boost is not required if TransferFunction is given coefficients.
  - We get better results if f' is largest near the endpoints [a,b].
    If f' is largest near the middle of an interval (eg e^{-x^2} when a<-3 and b>3)
    then the grid remains uniform. This is an issue with the way we approximate S!

  TODO Currently, transfer functions are basically useless if $f'$ is not extreme near the endpoints of its domain.
  (b/c we use cubic hermite interpolation at the endpoints).
  
  Note that a cubic polynomial
    p(x) = a_0 + a_1x + a_2x^2 + a_3x^3
  is monotone over R if and only if a_2^2 < 3a_1a_3. Maybe we can compute a polynomial minimizing
    int_a^b |f(t)-p(t)| dt
  such that
    p(a) = a, p(b) = b, and a_2^2 < 3a_1a_3.
  That would be a much better general purpose solution.
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
class TransferFunction
{
  TIN m_minArg, m_tableMaxArg, m_stepSize;
  /* This class is a polynomial approximating g inverse.
   * The identity transfer function is {-m_minArg/m_stepSize,1/m_stepSize,0,0} */
  __attribute__((aligned)) std::array<TIN,4> m_inverse_coefs = {{0,0,0,0}};
public:
  /* Set m_inverse_coefs equal to a vector that is (presumably) either the identity or came from a json file */
  TransferFunction(const std::array<TIN,4>& inv_coefs) { m_inverse_coefs = inv_coefs; }

  TransferFunction() = default;

  /* Build the coefficients in g_inv 
   * TODO can we use SFINAE to check if sqrt(TOUT) is defined? (ie TOUT cannot be LookupTable<...>) */
  template<typename TOUT>
  TransferFunction(const FunctionContainer<TIN,TOUT>& fc, TIN minArg, TIN tableMaxArg, TIN stepSize) : 
    m_minArg(minArg), m_tableMaxArg(tableMaxArg), m_stepSize(stepSize) {
#ifndef FUNC_USE_BOOST
    /* cause a compile time error because this constructor should never be called without Boost available */
    static_assert(sizeof(TIN) != sizeof(TIN), "Cannot generate a nonuniform grid without Boost verion 1.71.0 or higher");
#else
    using boost::math::quadrature::gauss_kronrod;
    using boost::math::differentiation::make_fvar;

    if(fc.autodiff1_fun == nullptr)
      throw std::invalid_argument("Error in func::TransferFunction. 1st derivative of function is needed to generate nonuniform grids but is null");

    /* std::function to return the first derivative of f */
    std::function<TOUT(TIN)> f_prime = [&fc](TIN x) -> TOUT {
      return (fc.autodiff1_fun)(make_fvar<TIN,1>(x)).derivative(1);
    };

    /* build the integrand. Notice that it is strictly positive */
    std::function<TIN(TIN)> integrand = [&f_prime](TIN x) -> TIN {
      return static_cast<TIN>(1.0)/static_cast<TIN>(sqrt(static_cast<TOUT>(1.0) + f_prime(x)*f_prime(x)));
    };

    TIN a = m_minArg;
    TIN b = m_tableMaxArg;

    /* integrate over [a,b] using adaptive quadrature & default tol of sqrt(epsilon_machine). */
    TIN c = gauss_kronrod<TIN, 15>::integrate(integrand, a, b);

    /* rescale the integrand to map [a,b] -> [a,b] */
    std::function<TIN(TIN)> g_prime = [&integrand,a,b,c](TIN x) -> TIN {
      return (b-a) * integrand(x)/c;
    };

    /* g(a)=a and g(b)=b so
     * g_inv_prime(a) = 1/g_prime(a) and g_inv_prime(b) = 1/g_prime(b) */
    TIN m0 = 1/g_prime(a);
    TIN m1 = 1/g_prime(b);

    /* Ensure monotonicity of the Hermite interpolating polynomial. Note that m0,m1 >= 0.
     * TODO this is sufficient to ensure p is monotone, but it's not necessary. Figure out how to allow m0>3 or m1>3. */
    m0 = (m0 > 3) ? 3 : m0;
    m1 = (m1 > 3) ? 3 : m1;

    /* Compute the Hermite interpolating polynomial p for g_inv satisfying
     * p(a)=a, p(b)=b, p'(a)=m0, p'(b)=m1
     * (this symbolic expression was computed by Matlab)
     * TODO maybe we can use this to get equations for the slopes m0 and m1? */
    m_inverse_coefs[0] = (a*b*(a + b - a*m1 - b*m0))/(a - b)/(a - b);
    m_inverse_coefs[1] = (a*a*m1 - 6*a*b + b*b*m0 + 2*a*b*m0 + 2*a*b*m1)/(a - b)/(a - b);
    m_inverse_coefs[2] = -(a*m0 - 3*b - 3*a + 2*a*m1 + 2*b*m0 + b*m1)/(a - b)/(a - b);
    m_inverse_coefs[3] = (m0 + m1 - 2)/(a - b)/(a - b);

    /* build the version of g_inv we'll actually use by encoding the
       underlying table's hash into the transfer function eval.
       This way, the only cost of using a transfer function is to lookup 4 stack allocated numbers */
    m_inverse_coefs[0] = m_inverse_coefs[0] - m_minArg;
    for(unsigned int i=0; i<4; i++)
      m_inverse_coefs[i] = m_inverse_coefs[i] / m_stepSize;
#endif
  }

  /* Horner's method */
  TIN inverse(TIN x) const {
    TIN sum = static_cast<TIN>(0);
    for (int k=3; k>0; k--)
      sum = x*(m_inverse_coefs[k] + sum);
    return sum + m_inverse_coefs[0];
  }
  
  /* Horner's method */
  TIN inverse_diff(TIN x) {
    TIN sum = static_cast<TIN>(0);
    for (int k=3; k>1; k--)
      sum = x*(k*m_inverse_coefs[k] + sum);
    return sum + m_inverse_coefs[1];
  }

  /* Use newton-raphson_iteration on p. Recall that each coef in
   * g was divided by h and we subtracted by m_minArg */
  TIN operator()(TIN x)
  {
    // This will have at least 0.9*std::numeric_limits<TIN>::digits digits of accuracy
    boost::uintmax_t maxit = 55;
    return boost::math::tools::newton_raphson_iterate(
        [this,x](const TIN z){ return std::make_tuple(inverse(z) - (x-m_minArg)/m_stepSize, inverse_diff(z));},
        x,m_minArg, m_tableMaxArg, 0.9*std::numeric_limits<TIN>::digits, maxit);
  }

  std::array<TIN,4> get_coefs() const { return m_inverse_coefs; }
  TIN min_arg() const { return m_minArg; }
  TIN max_arg() const { return m_tableMaxArg; }
};

/* print basic info about a TransferFunction */
template <typename TIN>
std::ostream& operator<<(std::ostream& out, const TransferFunction<TIN>& F){
  auto coefs = F.get_coefs();
  out << "degree 3 monotone Hermite interpolation with polynomial: \n";
  out << std::to_string(coefs[3]) << "x^3 + " << std::to_string(coefs[2]) << "x^2 + " <<
    std::to_string(coefs[1]) << "x + " << std::to_string(coefs[0]) << ". Defined over [" << F.min_arg() << "," << F.max_arg() << "].\n";
}

/* wraps operator<< */
template <typename TIN>
inline std::string to_string(const TransferFunction<TIN>& F) {
  std::ostringstream ss;
  ss << F;
  return ss.str();
}

}

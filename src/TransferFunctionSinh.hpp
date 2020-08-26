/*
  Given a function f and its domain [a,b], build a pair of functions
  g: [a,b] -> [a,b] and it's inverse g^{-1}.
  g is formally defined as
    g(x) = a + \frac{b-a}{c}\int_a^x\frac{dt}{\sqrt{1+[f'(t)]^2}}.
  where c = \int_a^b\frac{dt}{\sqrt{1+[f'(t)]^2}}

Notes:
  - To keep evaluations of g^{-1} quick, it's being approximated using 
  some form of inverse polynomial interpolation. This approximation ends
  up not being accurate enough to be useful in most cases so we'll change g
  to be the accurate inverse of our approximation to g^{-1} (that is, we'll be
  deepfrying g after approximating g^{-1}).
  - NUM_COEFS defines the number of coefficients used for approximating g^{-1}.

  Experimental atm.
  - TODO Make this class accept an arbitrary function which determines g from f
  - TODO for the functions that use derivative info, we could pass a slow g^{-1}
  and use g^{-1}' = (g'(g^{-1}))^{-1}
  - TODO quantify the conditioning of any given implementation of g^{-1}
  - TODO if g'=0 where we need 1/g' then swap out that method as well
  - TODO it would be nice if this class didn't neeed Armadillo to operate
 */
#pragma once
#include "config.hpp" // FUNC_USE_BOOST_AUTODIFF, FUNC_USE_ARMADILLO
#include <cmath> // sqrt
#include <string>
#include <limits> // std::numeric_limits<T>::max
#include <functional> // std::function

#ifndef FUNC_USE_BOOST_AUTODIFF
#error "TransferFunctionSinh needs boost version >= 1.71"
#endif

#ifndef FUNC_USE_ARMADILLO
#error "TransferFunctionSinh needs Armadillo"
#endif

#define ARMA_USE_CXX11
#include <armadillo>
#define FUNC_TRANSFER_FUNCTION_SOLVE_OPTS arma::solve_opts::refine

#include "TransferFunctionInterface.hpp"
#include "FunctionContainer.hpp"

#define BOOST_MATH_GAUSS_NO_COMPUTE_ON_DEMAND
#include <boost/math/quadrature/gauss_kronrod.hpp> // gauss_kronrod::integrate
#include <boost/math/tools/toms748_solve.hpp> // toms748_solve

static double tol = 1e-8;

template <typename IN_TYPE, unsigned int NUM_COEFS = 4>
class TransferFunctionSinh final : public TransferFunctionInterface<IN_TYPE>
{
  IMPLEMENT_TRANSFER_FUNCTION_INTERFACE(IN_TYPE);

  /* --- member vars --- */
  std::function<IN_TYPE(IN_TYPE)> m_g;
  // m_g will be approximated using Newton's method, hence the following
  // std::functions (even though m_g_inv isn't used for the member function g_inv())
  std::function<IN_TYPE(IN_TYPE)> m_g_inv;
  std::function<IN_TYPE(IN_TYPE)> m_g_inv_prime;

  __attribute__((aligned)) std::array<IN_TYPE,NUM_COEFS> m_inv_coefs;

  /* --- Private Functions for approximating g^{-1} --- */
  /* Approximate g^{-1} using inverse polynomial interpolation */
  static std::array<IN_TYPE,NUM_COEFS> inverse_poly_interp(
      std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp,
      IN_TYPE a, IN_TYPE b);

  /* Approximate g^{-1} using inverse polynomial interpolation
   * and specifying the slopes of inner points */
  static std::array<IN_TYPE,NUM_COEFS> inverse_poly_interior_slopes_interp(
      std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp,
      IN_TYPE a, IN_TYPE b);

  /* Approximate g^{-1} using inverse polynomial interpolation
   * and specifying slopes at function endpoints */
  static std::array<IN_TYPE,NUM_COEFS> inverse_hermite_interp(
      std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp,
      IN_TYPE a, IN_TYPE b);

  /* approximate g_inv by projecting g onto the corresponding polynomial space. */
  static std::array<IN_TYPE,NUM_COEFS> inverse_polynomial_projection(
    std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp,
    IN_TYPE a, IN_TYPE b);

  /* --- Private Helper functions --- */
  /* Make a std::array of coefs into a polynomial evaluated using horners
   * Note: std::array does deep copy's by default */
  template<unsigned long N>
  std::function<IN_TYPE(IN_TYPE)> make_horners(std::array<IN_TYPE,N> coefs);

  /* Use Newton's method to build g's inverse function,
   * or use bisection if things go south */
  static std::function<IN_TYPE(IN_TYPE)> newtons_inv(
      std::function<IN_TYPE(IN_TYPE)> const& g, std::function<IN_TYPE(IN_TYPE)> const& gp,
      IN_TYPE a, IN_TYPE b);

  /* fill a vector with N linearly spaced points with respect to g in [a,b].
   * assumes g is monotone on [a,b], g(a)=a, and g(b)=b */
  static arma::vec gspace(unsigned int N, std::function<IN_TYPE(IN_TYPE)> g,
      std::function<IN_TYPE(IN_TYPE)> gp, IN_TYPE a, IN_TYPE b);

public:
  /* public constructor */
  template<typename OUT_TYPE>
  TransferFunctionSinh(FunctionContainer<IN_TYPE,OUT_TYPE> *fc,
      IN_TYPE minArg, IN_TYPE maxArg, IN_TYPE stepSize);

  IN_TYPE g(IN_TYPE x) override { return m_g(x); }
  IN_TYPE g_inv(IN_TYPE x) override
  {
    IN_TYPE sum = x*m_inv_coefs[NUM_COEFS-1];
    for (int k=NUM_COEFS-2; k>0; k--)
      sum = x*(m_inv_coefs[k] + sum);
    return sum + m_inv_coefs[0];
  }

  void print_details(std::ostream& out) override
  {
    out << NUM_COEFS;
  }
};

template <typename IN_TYPE, unsigned int NUM_COEFS>
template <typename OUT_TYPE>
inline TransferFunctionSinh<IN_TYPE,NUM_COEFS>::TransferFunctionSinh(
    FunctionContainer<IN_TYPE,OUT_TYPE> *fc, IN_TYPE minArg, IN_TYPE maxArg, IN_TYPE stepSize) :
  TransferFunctionInterface<IN_TYPE>(fc,minArg,maxArg,stepSize)
{
  using boost::math::quadrature::gauss_kronrod;
  using boost::math::differentiation::make_fvar;

  // build a function to return the first derivative of f
  std::function<OUT_TYPE(IN_TYPE)> f_prime = [fc](IN_TYPE x) -> OUT_TYPE {
    return (fc->autodiff1_func)(make_fvar<IN_TYPE,1>(x)).derivative(1);
  };

  // build our transfer function, starting with the integrand.
  // This will be the definition of g described in the documentation.
  m_g = [f_prime](IN_TYPE x) -> IN_TYPE {
    return 1/((IN_TYPE) sqrt(1 + f_prime(x)*f_prime(x)));
  };

  // perform adaptive quadrature with a default tol of sqrt(epsilon).
  // Used to scale the actual g such that g(b)=b
  IN_TYPE c = gauss_kronrod<IN_TYPE, 15>::integrate(m_g, m_minArg, m_maxArg);

  // build m_g = a + \frac{b-a}{c}\int_a^x\frac{dt}{\sqrt{1+[f'(t)]^2}} which maps from [a,b] -> [a,b]
  // Note: m_g will be changed in a moment!
  m_g = [this,f_prime,c](IN_TYPE x) -> IN_TYPE {
    if(x <= m_minArg) return m_minArg; // boost gets upset if we don't do this
    return m_minArg + (m_maxArg - m_minArg)*gauss_kronrod<IN_TYPE, 15>::integrate(
        [f_prime](IN_TYPE t) -> IN_TYPE { return 1 / ((IN_TYPE)sqrt(1 + f_prime(t)*f_prime(t))); },
        m_minArg, x) / c;
  };

  // build m_g_prime. Only used for generating g_inv in a bit
  std::function<IN_TYPE(IN_TYPE)> g_prime = [this, f_prime, c](IN_TYPE x) -> IN_TYPE
  {
    return (m_maxArg-m_minArg) / (IN_TYPE) sqrt(1 + f_prime(x)*f_prime(x)) / c;
  };

  /* We'll build g_inv by using some form of inverse polynomial interpolantion.
     This is the experimental part, so it's modular. All we're checking for is some
     approximation of g_inv that is monotone. Currently there's no error est. */
  std::function<IN_TYPE(IN_TYPE)> formal_g_inv;

  // an array of each approximation method
  // rather than do an error est. we'll just hope the theoretically
  // better approximations are well conditioned
  const std::array<std::function<std::array<IN_TYPE,NUM_COEFS>(
      std::function<IN_TYPE(IN_TYPE)>,std::function<IN_TYPE(IN_TYPE)>, IN_TYPE, IN_TYPE
      )>, 2> approx_methods {inverse_poly_interior_slopes_interp, inverse_poly_interp};

  bool is_terrible = true; // assume every estimate is terrible unless proven otherwise
  for(unsigned int k=0; k<approx_methods.size(); k++){
    m_inv_coefs = approx_methods[k](m_g, g_prime, m_minArg, m_maxArg);

    // lets scope out the quality of this polynomial with a quick copy
    formal_g_inv = make_horners(m_inv_coefs);

    /* check if this version of g_inv is any good and make a fuss if it's terrible */
    // base conditions
    if(abs(formal_g_inv(m_minArg) - m_minArg) > tol || abs(formal_g_inv(m_maxArg) - m_maxArg) > tol){
      continue; // this estimate is terrible
    }
    // check g_inv at N linearly spaced points
    unsigned int const N = 50;
    for(unsigned int i=1; i<=N; i++){
      // check for monotonicity
      if(formal_g_inv((i-1)/(IN_TYPE)N) > formal_g_inv(i/(IN_TYPE)N))
        continue; // this estimate is also terrible
      // TODO error estimate
      //std::function<IN_TYPE(IN_TYPE)> slow_g_inv = newtons_inv(m_g, g_prime, m_minArg, m_maxArg);
    }
    is_terrible = false; // this estimate is at least passable
    break;
  }

  if(is_terrible)
    throw std::range_error("Every available polynomial approximation of the given transfer function "
        "using " + std::to_string(NUM_COEFS) + " coefficients is too poorly conditioned");

  // compute formal_g_inv's derivative
  std::array<IN_TYPE,NUM_COEFS-1> inv_prime_coefs;
  for(unsigned int j=1; j < NUM_COEFS; j++)
    inv_prime_coefs[j-1] = j*m_inv_coefs[j];
  m_g_inv_prime = make_horners(inv_prime_coefs);

  /* build the real version of g_inv by encoding the
     underlying table's hash into the transfer function eval. This
     will obfuscate our code, but now the only price of using
     nonuniform tables is an indirection */
  m_inv_coefs[0] = m_inv_coefs[0] - m_minArg;
  for(unsigned int i=0; i<NUM_COEFS; i++)
    m_inv_coefs[i] = m_inv_coefs[i] / stepSize;

  /* Now that we have a fast approx to g_inv, we'll make it more "accurate" by
     setting our original g to g_inv_inv. This Newton's method is the reason why
     we computed all those derivatives earlier */
  m_g = newtons_inv(formal_g_inv, m_g_inv_prime, m_minArg, m_maxArg);
}


template <typename IN_TYPE, unsigned int NUM_COEFS>
inline std::array<IN_TYPE,NUM_COEFS> TransferFunctionSinh<IN_TYPE,NUM_COEFS>::inverse_poly_interp(
    std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp, IN_TYPE a, IN_TYPE b)
{
  using arma::span;
  // check if this is possible
  if(NUM_COEFS < 2)
    throw std::invalid_argument("inverse_poly_coefs() must"
      " produce at least 2 polynomial coefficients");

  // generate vandermonde system
  arma::mat A = arma::ones(NUM_COEFS,NUM_COEFS);
  A(span(0,NUM_COEFS-1), 1)=arma::linspace(a,b,NUM_COEFS);
  for(unsigned int i=2; i<NUM_COEFS; i++)
    A(span(0,NUM_COEFS-1),i) = A(span(0,NUM_COEFS-1),i-1) % A(span(0,NUM_COEFS-1),1);

  // generate solution vector
  arma::vec y = arma::ones(NUM_COEFS);
  // note: gspace returns a vector of length NUM_COEFS whos points are all
  // between 0 and 1
  y.rows(0,NUM_COEFS-1) = gspace(NUM_COEFS, g, gp, a, b);

  y = arma::solve(A,y,FUNC_TRANSFER_FUNCTION_SOLVE_OPTS);

  // move from arma's vector to a std::array<IN_TYPE,NUM_COEFS>
  auto coefs = std::array<IN_TYPE,NUM_COEFS>();
  for(unsigned int i = 0; i < NUM_COEFS; i++)
    coefs[i] = y[i];

  return coefs;
}

/* approximate g_inv with inverse polynomial interpolation and by specifying
 * slopes at interior points as 1/g_prime */
template <typename IN_TYPE, unsigned int NUM_COEFS>
inline std::array<IN_TYPE,NUM_COEFS> TransferFunctionSinh<IN_TYPE,NUM_COEFS>::inverse_poly_interior_slopes_interp(
    std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp, IN_TYPE a, IN_TYPE b)
{
  using arma::span;
  // generate vandermonde system
  // assert the user is asking for an even number of coefficients
  if(NUM_COEFS % 2 != 0)
    throw std::invalid_argument("NUM_COEFS is odd"
      " but inverse_poly_interior_slopes_coefs() can only produce an"
      " even number of polynomial coefficients");

  unsigned int M = NUM_COEFS/2+1; // number of unique points being sampled from
  arma::mat A = arma::ones(NUM_COEFS,NUM_COEFS);

  A(span(0,M-1), 1)=arma::linspace(a,b,M);
  for(unsigned int i=2; i<NUM_COEFS; i++)
    A(span(0,M-1),i) = A(span(0,M-1),i-1) % A(span(0,M-1),1);

  // set the bottom half to derivative values
  A(span(M,NUM_COEFS-1),0) = arma::zeros(M-2,1);
  for(unsigned int i=1; i<NUM_COEFS; i++)
    A(span(M,NUM_COEFS-1),i) = i*A(span(1,M-2),i-1);

  // generate solution vector
  arma::vec y = arma::ones(NUM_COEFS);
  y.rows(0,M-1) = gspace(M, g, gp, a, b);
  // using 1/g'(x_i) as an approximation for the slopes of the inverse functions
  for(int i=1; i<M-1; i++){
    y[M-1+i] = 1.0/gp(y[i]); // requires f(y[i])\neq\pm\infy
  }

  y = arma::solve(A,y,FUNC_TRANSFER_FUNCTION_SOLVE_OPTS);

  // move from arma's vector to a std::array<IN_TYPE,NUM_COEFS>
  auto coefs = std::array<IN_TYPE,NUM_COEFS>();
  for(unsigned int i = 0; i < NUM_COEFS; i++)
    coefs[i] = y[i];

  return coefs;
}

template <typename IN_TYPE, unsigned int NUM_COEFS>
inline std::array<IN_TYPE,NUM_COEFS> TransferFunctionSinh<IN_TYPE,NUM_COEFS>::inverse_polynomial_projection(
    std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp, IN_TYPE a, IN_TYPE b)
{
  // move from arma's vector to a std::array<IN_TYPE,NUM_COEFS>
  auto coefs = std::array<IN_TYPE,NUM_COEFS>();
  for(unsigned int i = 0; i < NUM_COEFS; i++)
    coefs[i] = 1;
  
  return coefs;
}

/* approximate g_inv with inverse hermite interpolation. ie specify
 * slopes at the endpoints as 1/g_prime */
template <typename IN_TYPE, unsigned int NUM_COEFS>
inline std::array<IN_TYPE,NUM_COEFS> TransferFunctionSinh<IN_TYPE,NUM_COEFS>::inverse_hermite_interp(
    std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp, IN_TYPE a, IN_TYPE b)
{
  // TODO this is outdated and I don't there are
  // many situations where this will be useful...
  using arma::span;
  // generate vandermonde system
  // assert the user is asking for at least 4 coefs
  if(NUM_COEFS < 4)
      throw std::invalid_argument("NUM_COEFS is too small."
      " inverse_hermite_interp() can only produce 4 or more coefs");
  unsigned int M = NUM_COEFS-2; // number of unique points being sampled from
  arma::mat A = arma::ones(NUM_COEFS,NUM_COEFS);

  A(span(0,M-1), 1)=arma::linspace(0,1,M);
  for(unsigned int i=2; i<NUM_COEFS; i++)
    A(span(0,M-1),i) = A(span(0,M-1),i-1) % A(span(0,M-1),1);

  // set the bottom half to derivative values.
  A(span(M,NUM_COEFS-1),0) = arma::zeros(2,1);
  for(unsigned int i=1; i<NUM_COEFS; i++){
    A(M,i) = i*A(0,i-1);
    A(M+1,i) = i*A(M-1,i-1);
  }

  // generate solution vector
  arma::vec y = arma::ones(NUM_COEFS);
  y.rows(0,M-1) = gspace(M, g, gp, a, b);
  // using 1/g'(x_i) as an approximation for the slopes of the inverse functions
  y[M] = 1.0/gp(y[0]); // requires f(y[i])\neq\pm\infy
  y[M+1] = 1.0/gp(y[M-1]); // requires f(y[i])\neq\pm\infy

  y = arma::solve(A,y,FUNC_TRANSFER_FUNCTION_SOLVE_OPTS);

  // move from arma's vector to a std::array<IN_TYPE,NUM_COEFS>
  auto coefs = std::array<IN_TYPE,NUM_COEFS>(new IN_TYPE[NUM_COEFS]);
  for(unsigned int i = 0; i < NUM_COEFS; i++)
    coefs[i] = y[i];
  
  return coefs;
}

template <typename IN_TYPE, unsigned int NUM_COEFS>
template <unsigned long N>
inline std::function<IN_TYPE(IN_TYPE)> TransferFunctionSinh<IN_TYPE,NUM_COEFS>::make_horners(
    std::array<IN_TYPE,N> coefs)
{
  // We need the coefs by reference so that make_horners doesn't
  // return a function with a dangling stack value. But the lambda takes
  // a 
  return [coefs](IN_TYPE x) -> IN_TYPE {
    IN_TYPE sum = x*coefs[N-1];
    for (unsigned int k = N-2; k > 0; k--)
      sum = x*(coefs[k] + sum);
    return sum + coefs[0];
  };
}

template <typename IN_TYPE, unsigned int NUM_COEFS>
inline std::function<IN_TYPE(IN_TYPE)> TransferFunctionSinh<IN_TYPE,NUM_COEFS>::newtons_inv(
    std::function<IN_TYPE(IN_TYPE)> const& g, std::function<IN_TYPE(IN_TYPE)> const& gp, IN_TYPE a, IN_TYPE b)
{
  return [g, gp, a, b](IN_TYPE z) -> IN_TYPE {
    using boost::math::tools::eps_tolerance;
    using boost::math::tools::toms748_solve;

    const unsigned int MAX_NEWTON_IT = 20; // might be a bit large
    boost::uintmax_t const MAX_BISECTION = 54;

    // Find what g maps to x
    IN_TYPE x = z;
    IN_TYPE x0;
    unsigned int NEWTON_IT = 0;
    // Use Newton's method if we can get away with it. o.w. resort to bisection
    do{
      NEWTON_IT += 1;
      x0 = x;
      if(gp == NULL || gp(x) == 0.0 || NEWTON_IT > MAX_NEWTON_IT){
        // if g prime is undefined, do a hard switch to bisection
        boost::uintmax_t MAX_IT = MAX_BISECTION;
        x0 = x = toms748_solve(
            [g, &z](IN_TYPE h) -> IN_TYPE { return g(h) - z; }, // shift g
            a, b, a - z, b - z, eps_tolerance<IN_TYPE>(), MAX_IT).first;
      }else{
        x = x-(g(x)-z)/gp(x);
      }
    }while(abs(x0-x) > tol);

    return x;
  };
}

template <typename IN_TYPE, unsigned int NUM_COEFS>
inline arma::vec TransferFunctionSinh<IN_TYPE,NUM_COEFS>::gspace(unsigned int N,
    std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp, IN_TYPE a, IN_TYPE b)
{
  using boost::math::tools::eps_tolerance;
  using boost::math::tools::toms748_solve;

  arma::vec const linear_pts = arma::linspace(a,b,N);
  arma::vec v = arma::ones(N,1);
  v[0] = a;
  const unsigned int MAX_NEWTON_IT = 20;
  boost::uintmax_t const MAX_BISECTION = 54;

  for(int i=1; i<N-1; i++){
    // Find what g maps to linear_pts[i]
    IN_TYPE x0;
    IN_TYPE x = linear_pts[i];
    unsigned int NEWTON_IT = 0;
    // Use Newton's method if we can get away with it. o.w. resort to bisection
    do{
      NEWTON_IT += 1;
      x0 = x;
      if(gp == NULL || gp(x) == 0.0 || x < a || x > b || NEWTON_IT > MAX_NEWTON_IT){
        // if g prime is undefined, do a hard switch to bisection
        boost::uintmax_t MAX_IT = MAX_BISECTION;
        x0 = x = toms748_solve(
            [g, linear_pts, &i](IN_TYPE z) -> IN_TYPE { return g(z) - (IN_TYPE) linear_pts[i]; }, // shift g
            a, b, a-(IN_TYPE)linear_pts[i], b-(IN_TYPE)linear_pts[i], eps_tolerance<IN_TYPE>(), MAX_IT).first;
      }else{
        x = x-(g(x)-linear_pts[i])/gp(x);
      }
    }while(abs(x0-x) > tol);
    v[i]=x;
  }
  v[N-1] = b;
  return v;
}

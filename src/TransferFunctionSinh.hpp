/*
  builds a pair of functions g: [0,1] -> [0,1] and it's inverse g_inv from 
  a given function f and its domain [a,b].
    g(x) = \int_a^{a+x*(b-a)} \frac{dt}{sqrt(1+[f'(t)]^2)},
  and to keep evaluations of g^{-1} quick, it's being approximated using 
  some form of inverse polynomial interpolation.
  This class has been designed such that each method of approximating g^{-1} can be
  easily swapped.
  N defines the number of coefficients used to approximate g with inverse poly interp.

  Experimental atm.
  - TODO quantify the conditioning of any given implementation of g^{-1}
  - TODO check the monotonicity of g_inv and that g_inv(1)=1. Especially need
  0 <= g_inv(x) <= 1 for every x in [0,1] if we want to avoid strange seg faults
  - TODO if g'=0 where we need 1/g' then swap out that method as well
  - TODO decide on pros/cons of certain options.
  - TODO decide if helper functions should be in base class. Currently they're
  fairly general
 */
#pragma once
#include <cmath> // sqrt
#include <string>
#include <limits> // std::numeric_limits<T>::max
#include <functional> // std::function
#define ARMA_USE_CXX11
#include <armadillo>
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

  /* --- More member vars --- */
  std::string m_method_of_approx;
  std::function<IN_TYPE(IN_TYPE)> m_g;
  // using inverse polynomial interpolation and horners method for g_inv
  __attribute__((aligned)) std::unique_ptr<IN_TYPE[]> m_inv_coefs;

  /* --- Private Functions for approximating g^{-1} --- */
  /* Approximate g^{-1} using inverse polynomial interpolation */
  std::unique_ptr<IN_TYPE[]> inverse_poly_interp(
      std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp);

  /* Approximate g^{-1} using inverse polynomial interpolation
   * and specifying the slopes of inner points */
  std::unique_ptr<IN_TYPE[]> inverse_poly_interior_slopes_interp(
      std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp);

  /* Approximate g^{-1} using inverse polynomial interpolation
   * and specifying slopes at function endpoints */
  std::unique_ptr<IN_TYPE[]> inverse_hermite_interp(
      std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp);

  /* approximate g_inv by projecting g onto the corresponding polynomial space. */
  std::unique_ptr<IN_TYPE[]> inverse_polynomial_projection(
    std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp);

  /* Use Newton's method to build g's inverse function,
   * or use bisection if things go south */
  std::function<IN_TYPE(IN_TYPE)> newtons_inv(
      std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp);

  /* --- Private Helper functions --- */
  /* fill a vector with N linearly spaced points with respect to g in [0,1].
   * assumes g is monotone on [0,1], g(0)=0, and g(1)=1 */
  arma::vec gspace(unsigned int N, std::function<IN_TYPE(IN_TYPE)> g,
      std::function<IN_TYPE(IN_TYPE)> gp = NULL);

public:
  /* public constructor */
  template<typename OUT_TYPE>
  TransferFunctionSinh(FunctionContainer<IN_TYPE,OUT_TYPE> *fc, IN_TYPE minArg, IN_TYPE maxArg, IN_TYPE stepSize);

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
    out << "arcsinh transfer function approximating g_inv with ";
    out << m_method_of_approx;
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
  std::function<OUT_TYPE(IN_TYPE)> f_prime = [&fc](IN_TYPE x) -> OUT_TYPE {
    return (fc->autodiff1_func)(make_fvar<IN_TYPE,1>(x)).derivative(1);
  };

  // build our transfer function, starting with the integrand
  std::function<IN_TYPE(IN_TYPE)> temp_g = [f_prime](IN_TYPE x) -> IN_TYPE {
    return 1/((IN_TYPE) sqrt(1 + f_prime(x)*f_prime(x)));
  };

  // perform adaptive quadrature with a default tol of sqrt(epsilon).
  // Used to scale the actual g such that g(b)=b
  IN_TYPE c = gauss_kronrod<IN_TYPE, 15>::integrate(temp_g, m_minArg, m_maxArg);

  // build temp_g = a + \frac{b-a}{c}\int_a^x\frac{dt}{\sqrt{1+[f'(t)]^2}} which maps from [a,b] -> [a,b]
  temp_g = [this,f_prime,c](IN_TYPE x) -> IN_TYPE {
    if(x == m_minArg) return m_minArg; // boost gets upset if we don't do this
    return m_minArg + (m_maxArg - m_minArg)*gauss_kronrod<IN_TYPE, 15>::integrate(
        [f_prime](IN_TYPE t) -> IN_TYPE { return 1 / ((IN_TYPE)sqrt(1 + f_prime(t)*f_prime(t))); },
        m_minArg, x) / c;
  };

  // build temp_g_prime
  std::function<IN_TYPE(IN_TYPE)> temp_g_prime = [this, f_prime, c](IN_TYPE x) -> IN_TYPE
  {
    return (m_maxArg-m_minArg) / (IN_TYPE) sqrt(1 + f_prime(x)*f_prime(x)) / c;
  };

  /* We'll build g_inv by using some form of inverse polynomial interpolant.
     This is the experimental part, and so it has been made modular
     want a fast g_inv, but we also want an accurate g_inv in order to make
     good use of the nonuniform grid */
  // TODO make a vector of usable coefficients and then check which one is the best
  //std::unique_ptr<IN_TYPE[]> m_inv_coefs = inverse_poly_interior_slopes_interp(temp_g, temp_g_prime);
  m_inv_coefs = inverse_poly_interp(temp_g, temp_g_prime);

  // make a copy of inv_coefs to give to a functor
  std::array<IN_TYPE,NUM_COEFS> inv_coefs_copy;
  for(unsigned int j=0; j<NUM_COEFS; j++)
    inv_coefs_copy[j] = m_inv_coefs[j];

  // set temp_g_inv using horners method
  auto temp_g_inv = [inv_coefs_copy](IN_TYPE x) -> IN_TYPE
  {
    IN_TYPE sum = x*inv_coefs_copy[NUM_COEFS-1];
    for (int k=NUM_COEFS-2; k>0; k--)
      sum = x*(inv_coefs_copy[k] + sum);
    return sum + inv_coefs_copy[0];
  };

  // compute temp_g_inv's derivative
  std::array<IN_TYPE,NUM_COEFS-1> inv_prime_coefs;
  for(unsigned int j=1; j < NUM_COEFS; j++)
    inv_prime_coefs[j-1] = j*m_inv_coefs[j];

  // set g_inv_prime using horners method
  auto temp_g_inv_prime = [inv_prime_coefs](IN_TYPE x) -> IN_TYPE
  {
    IN_TYPE sum = x*inv_prime_coefs[NUM_COEFS-2];
    for (int k=NUM_COEFS-3; k>0; k--)
      sum = x*(inv_prime_coefs[k] + sum);
    return sum + inv_prime_coefs[0];
  };

  // check if this version of g_inv is any good and make a fuss if it's terrible
  //
  // Slow but accurate approximation of g_inv
  // std::function<IN_TYPE(IN_TYPE)> slow_g_inv = newtons_inv(temp_g, temp_g_prime);
  //{
  //  unsigned int N = 20;
  //  long double error_est;
  //  // base conditions
  //  if(abs(m_g_inv(0)) > tol || abs(m_g_inv(1) - 1) > tol){
  //    return std::numeric_limits<double>::max;
  //  }
  //  // check g_inv at N linearly spaced points
  //  for(unsigned int i=1; i<=N; i++){
  //    // check for monotonicity
  //    if((*m_g_inv)((i-1)/(IN_TYPE)N) > (*m_g_inv)(i/(IN_TYPE)N))
  //      return std::numeric_limits<double>::max;

  //    // estimate the error
  //    long double exact  = slow_g_inv(i/(IN_TYPE)N);
  //    long double approx = (*m_g_inv)(i/(IN_TYPE)N);
  //    error_est += 2*(abs(exact - approx)/(exact+approx));
  //  }
  //  // return the average relative error
  //  return error_est / N;
  //}

    
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
  m_g = [this,temp_g_inv,temp_g_inv_prime](IN_TYPE z) -> IN_TYPE {
    // FunC tables are often generated with points a bit past their right endpoint
    if(z >= m_maxArg)
      return z;

    // approx g with newton's method on g_inv
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
      if(temp_g_inv_prime(x) == 0.0 || NEWTON_IT > MAX_NEWTON_IT){
        // if g prime is undefined, do a hard switch to bisection
        boost::uintmax_t MAX_IT = MAX_BISECTION;
        x0 = x = toms748_solve(
            [this, temp_g_inv, z](IN_TYPE h) -> IN_TYPE { return temp_g_inv(h) - z; }, // shift g
            m_minArg, m_maxArg, m_minArg - z, m_maxArg - z, eps_tolerance<IN_TYPE>(), MAX_IT).first;
      }else{
        x = x-(temp_g_inv(x)-z)/temp_g_inv_prime(x);
      }
    }while(abs(x0-x) > tol);

    return x;
  };
}

template <typename IN_TYPE, unsigned int NUM_COEFS>
inline std::unique_ptr<IN_TYPE[]> TransferFunctionSinh<IN_TYPE,NUM_COEFS>::inverse_poly_interp(
    std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp)
{
  using arma::span;
  // check if this is possible
  if(NUM_COEFS < 2)
    throw std::invalid_argument("inverse_poly_coefs() must"
      " produce at least 2 polynomial coefficients");
  m_method_of_approx = std::string("inverse polynomial interpolation");

  // generate vandermonde system
  arma::mat A = arma::ones(NUM_COEFS,NUM_COEFS);
  A(span(0,NUM_COEFS-1), 1)=arma::linspace(m_minArg,m_maxArg,NUM_COEFS);
  for(unsigned int i=2; i<NUM_COEFS; i++)
    A(span(0,NUM_COEFS-1),i) = A(span(0,NUM_COEFS-1),i-1) % A(span(0,NUM_COEFS-1),1);

  // generate solution vector
  arma::vec y = arma::ones(NUM_COEFS);
  // note: gspace returns a vector of length NUM_COEFS whos points are all
  // between 0 and 1
  y.rows(0,NUM_COEFS-1) = gspace(NUM_COEFS, g, gp);

  y = arma::solve(A,y);

  // move from arma's vector to a std::unique_ptr<double[]>
  auto coefs = std::unique_ptr<IN_TYPE[]>(new IN_TYPE[NUM_COEFS]);
  for(unsigned int i = 0; i < NUM_COEFS; i++)
    coefs[i] = y[i];

  return coefs;
}

/* approximate g_inv with inverse polynomial interpolation and by specifying
 * slopes at interior points as 1/g_prime */
template <typename IN_TYPE, unsigned int NUM_COEFS>
inline std::unique_ptr<IN_TYPE[]> TransferFunctionSinh<IN_TYPE,NUM_COEFS>::inverse_poly_interior_slopes_interp(
    std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp)
{
  using arma::span;
  // generate vandermonde system
  // assert the user is asking for an even number of coefficients
  if(NUM_COEFS % 2 != 0)
    throw std::invalid_argument("NUM_COEFS is odd"
      " but inverse_poly_interior_slopes_coefs() can only produce an"
      " even number of polynomial coefficients");

  m_method_of_approx = std::string("inverse polynomial interpolation specifying "
      "slopes of inner points as 1/g'");
  unsigned int M = NUM_COEFS/2+1; // number of unique points being sampled from
  arma::mat A = arma::ones(NUM_COEFS,NUM_COEFS);

  A(span(0,M-1), 1)=arma::linspace(m_minArg,m_maxArg,M);
  for(unsigned int i=2; i<NUM_COEFS; i++)
    A(span(0,M-1),i) = A(span(0,M-1),i-1) % A(span(0,M-1),1);

  // set the bottom half to derivative values
  A(span(M,NUM_COEFS-1),0) = arma::zeros(M-2,1);
  for(unsigned int i=1; i<NUM_COEFS; i++)
    A(span(M,NUM_COEFS-1),i) = i*A(span(1,M-2),i-1);

  // generate solution vector
  arma::vec y = arma::ones(NUM_COEFS);
  y.rows(0,M-1) = gspace(M, g, gp);
  // using 1/g'(x_i) as an approximation for the slopes of the inverse functions
  for(int i=1; i<M-1; i++){
    y[M-1+i] = 1.0/gp(y[i]); // requires f(y[i])\neq\pm\infy
  }

  y = arma::solve(A,y);

  // move from arma's vector to a std::unique_ptr<double[]>
  auto coefs = std::unique_ptr<IN_TYPE[]>(new IN_TYPE[NUM_COEFS]);
  for(unsigned int i = 0; i < NUM_COEFS; i++)
    coefs[i] = y[i];

  return coefs;
}

template <typename IN_TYPE, unsigned int NUM_COEFS>
inline std::unique_ptr<IN_TYPE[]> TransferFunctionSinh<IN_TYPE,NUM_COEFS>::inverse_polynomial_projection(
    std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp)
{
  // TODO
  // move from arma's vector to a std::unique_ptr<double[]>
  auto coefs = std::unique_ptr<IN_TYPE[]>(new IN_TYPE[NUM_COEFS]);
  for(unsigned int i = 0; i < NUM_COEFS; i++)
    coefs[i] = 1;
  
  return coefs;
}

/* approximate g_inv with inverse hermite interpolation. ie specify
 * slopes at the endpoints as 1/g_prime */
template <typename IN_TYPE, unsigned int NUM_COEFS>
inline std::unique_ptr<IN_TYPE[]> TransferFunctionSinh<IN_TYPE,NUM_COEFS>::inverse_hermite_interp(
    std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp)
{
  // TODO this is outdated and I don't there are
  // many situations where this will be useful...
  using arma::span;
  // generate vandermonde system
  // assert the user is asking for at least 4 coefs
  if(NUM_COEFS < 4)
      throw std::invalid_argument("NUM_COEFS is too small."
      " inverse_hermite_interp() can only produce 4 or more coefs");
  m_method_of_approx = std::string("inverse polynomial interpolation specifying "
      "slopes of endpoints as 1/g'");
  unsigned int M = NUM_COEFS-2; // number of unique points being sampled from
  arma::mat A = arma::ones(NUM_COEFS,NUM_COEFS);

  A(span(0,M-1), 1)=arma::linspace(0,1,M);
  for(unsigned int i=2; i<NUM_COEFS; i++)
    A(span(0,M-1),i) = A(span(0,M-1),i-1) % A(span(0,M-1),1);

  // set the bottom half to derivative values.
  // TODO avoid the for loop
  A(span(M,NUM_COEFS-1),0) = arma::zeros(2,1);
  for(unsigned int i=1; i<NUM_COEFS; i++){
    A(M,i) = i*A(0,i-1);
    A(M+1,i) = i*A(M-1,i-1);
  }

  // generate solution vector
  arma::vec y = arma::ones(NUM_COEFS);
  y.rows(0,M-1) = gspace(M, g, gp);
  // using 1/g'(x_i) as an approximation for the slopes of the inverse functions
  y[M] = 1.0/gp(y[0]); // requires f(y[i])\neq\pm\infy
  y[M+1] = 1.0/gp(y[M-1]); // requires f(y[i])\neq\pm\infy

  y = arma::solve(A,y);

  if(abs(y[0]) > tol)
    throw std::runtime_error("inverse_hermite_interp() failed to"
      " solve for the inverse");

  // move from arma's vector to a std::unique_ptr<double[]>
  auto coefs = std::unique_ptr<IN_TYPE[]>(new IN_TYPE[NUM_COEFS]);
  for(unsigned int i = 0; i < NUM_COEFS; i++)
    coefs[i] = y[i];
  
  return coefs;
}

template <typename IN_TYPE, unsigned int NUM_COEFS>
inline std::function<IN_TYPE(IN_TYPE)> TransferFunctionSinh<IN_TYPE,NUM_COEFS>::newtons_inv(
    std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp)
{
  return [this, &g, &gp](IN_TYPE z) -> IN_TYPE {
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
            [this, &g, &z](IN_TYPE h) -> IN_TYPE { return g(h) - z; }, // shift g
            m_minArg, m_maxArg, m_minArg-z, m_maxArg - z, eps_tolerance<IN_TYPE>(), MAX_IT).first;
      }else{
        x = x-(g(x)-z)/gp(x);
      }
    }while(abs(x0-x) > tol);

    return x;
  };
}

/* fill a vector with N linearly spaced points with respect to g in [0,1].
 * assumes g is monotone on [0,1], g(0)=0, and g(1)=1 */
template <typename IN_TYPE, unsigned int NUM_COEFS>
inline arma::vec TransferFunctionSinh<IN_TYPE,NUM_COEFS>::gspace(unsigned int N,
    std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp)
{
  using boost::math::tools::eps_tolerance;
  using boost::math::tools::toms748_solve;

  arma::vec const linear_pts = arma::linspace(m_minArg,m_maxArg,N);
  arma::vec v = arma::ones(N,1);
  v[0] = m_minArg;
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
      if(gp == NULL || gp(x) == 0.0 || x < m_minArg || x > m_maxArg || NEWTON_IT > MAX_NEWTON_IT){
        // if g prime is undefined, do a hard switch to bisection
        boost::uintmax_t MAX_IT = MAX_BISECTION;
        x0 = x = toms748_solve(
            [&g, &linear_pts, &i](IN_TYPE z) -> IN_TYPE { return g(z) - linear_pts[i]; }, // shift g
            m_minArg, m_maxArg, m_minArg-linear_pts[i], m_maxArg-linear_pts[i], eps_tolerance<IN_TYPE>(), MAX_IT).first;
      }else{
        x = x-(g(x)-linear_pts[i])/gp(x);
      }
    }while(abs(x0-x) > tol);
    v[i]=x;
  }
  v[N-1] = m_maxArg;
  return v;
}

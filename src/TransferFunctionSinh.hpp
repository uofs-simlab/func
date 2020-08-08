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
#include "TransferFunction.hpp"
#include "FunctionContainer.hpp"

#define BOOST_MATH_GAUSS_NO_COMPUTE_ON_DEMAND
#include <boost/math/quadrature/gauss_kronrod.hpp> // gauss_kronrod::integrate
#include <boost/math/tools/toms748_solve.hpp> // toms748_solve

template <typename IN_TYPE, unsigned int N>
class HornersFunctor : public LightweightFunctor<IN_TYPE>
{
  __attribute__((aligned)) std::unique_ptr<IN_TYPE[]> m_coefs;

public:
  HornersFunctor(std::unique_ptr<IN_TYPE[]>& coefs) :
    m_coefs(std::move(coefs)) {}

  IN_TYPE operator()(IN_TYPE x) override
  {
      IN_TYPE sum = x*m_coefs[N-1];
      for (int k=N-2; k>0; k--)
        sum = x*(m_coefs[k] + sum);
      return sum + m_coefs[0];
  }
};

/* convenience function for wrapping lambdas */
template <typename IN_TYPE, class F>
class LambdaWrapper : public LightweightFunctor<IN_TYPE>
{
  F m_f;
public:
  LambdaWrapper(F f) : m_f(f) {}

  IN_TYPE operator()(IN_TYPE x){ return m_f(x); }
};

static double tol = 1e-8;

template <typename IN_TYPE>
class TransferFunctionSinh final : public TransferFunction<IN_TYPE>
{
  INHERIT_TRANSFER_FUNCTION(IN_TYPE);

  /* --- More member vars --- */
  std::unique_ptr<LightweightFunctor<IN_TYPE>> m_g_inv_prime;
  std::string m_method_of_approx;

  /* --- Private Functions for approximating g^{-1} --- */
  /* Approximate g^{-1} using inverse polynomial interpolation */
  std::unique_ptr<IN_TYPE[]> inverse_poly_interp(unsigned int num_coefs,
      std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp);

  /* Approximate g^{-1} using inverse polynomial interpolation
   * and specifying the slopes of inner points */
  std::unique_ptr<IN_TYPE[]> inverse_poly_interior_slopes_interp(
      unsigned int num_coefs, std::function<IN_TYPE(IN_TYPE)> g,
      std::function<IN_TYPE(IN_TYPE)> gp);

  /* Approximate g^{-1} using inverse polynomial interpolation
   * and specifying slopes at function endpoints */
  std::unique_ptr<IN_TYPE[]> inverse_hermite_interp(unsigned int num_coefs,
      std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp);

  /* approximate g_inv by projecting g onto the corresponding polynomial space. */
  std::unique_ptr<IN_TYPE[]> inverse_polynomial_projection(
    unsigned int num_coefs,
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
  TransferFunctionSinh(FunctionContainer<IN_TYPE,OUT_TYPE> *fc, IN_TYPE a, IN_TYPE b);

  void print_details(std::ostream& out) override
  {
    out << "arcsinh transfer function approximating g_inv with ";
    out << m_method_of_approx;
  }
};

template <typename IN_TYPE>
template <typename OUT_TYPE>
inline TransferFunctionSinh<IN_TYPE>::TransferFunctionSinh(FunctionContainer<IN_TYPE,OUT_TYPE> *fc, IN_TYPE a, IN_TYPE b) :
  TransferFunction<IN_TYPE>(a,b)
{
  using boost::math::quadrature::gauss_kronrod;
  using boost::math::differentiation::make_fvar;

  // build a function to return the first derivative of f
  std::function<OUT_TYPE(IN_TYPE)> f_prime = [&fc](IN_TYPE x) -> OUT_TYPE {
    return (fc->autodiff1_func)(make_fvar<IN_TYPE,1>(x)).derivative(1);
  };

  // build our transfer function
  // we'll be adjusting temp_g: [a,b] -> [a,b]
  std::function<IN_TYPE(IN_TYPE)> temp_g = [f_prime](IN_TYPE x) -> IN_TYPE {
    return 1/((IN_TYPE) sqrt(1 + f_prime(x)*f_prime(x)));
  };

  // perform adaptive quadrature with a default tol of sqrt(epsilon)
  IN_TYPE c = gauss_kronrod<IN_TYPE, 15>::integrate(temp_g, m_minArg, m_maxArg);

  // g:[0,1] -> [0,1] integrates temp_g over all of [a,b] and scales the answer
  // such that g(1) = 1. We'll re-adjust this later
  temp_g = [this,f_prime,c](IN_TYPE x) -> IN_TYPE {
    if(x == 0.0) return 0.0; // boost gets upset if we don't do this
    return gauss_kronrod<IN_TYPE, 15>::integrate(
        [f_prime,c](IN_TYPE t) -> IN_TYPE { return 1/((IN_TYPE)sqrt(1 + f_prime(t)*f_prime(t))) / c; },
        m_minArg, m_minArg+x*(m_maxArg-m_minArg));
  };

  // build g prime
  std::function<IN_TYPE(IN_TYPE)> g_prime = [this, f_prime, c](IN_TYPE x) -> IN_TYPE
  {
    IN_TYPE lerp = m_minArg + x*(m_maxArg-m_minArg);
    return (m_maxArg-m_minArg) / ((IN_TYPE) sqrt(1 + f_prime(lerp)*f_prime(lerp))) / c;
  };
  
  // use slow Newton's method as a slow but accurate initial approximation of g_inv
  std::function<IN_TYPE(IN_TYPE)> slow_g_inv = newtons_inv(temp_g, g_prime);

  /* We'll build g_inv by using some form of inverse polynomial interpolant.
     This is the experimental part, and so it has been made modular
     want a fast g_inv, but we also want an accurate g_inv in order to make
     good use of the nonuniform grid */
  constexpr unsigned int num_coefs = 2;
  for(unsigned int i=0; i < 1; i++){
    //auto inv_coefs = inverse_poly_interp(num_coefs[i], temp_g, g_prime);
    //std::unique_ptr<IN_TYPE[]> inv_coefs = inverse_poly_interior_slopes_interp(num_coefs, temp_g, g_prime);
    std::unique_ptr<IN_TYPE[]> inv_coefs = std::unique_ptr<IN_TYPE[]>(new IN_TYPE[]{0, 1});

    // compute this approximation's derivative while we still have access to inv_coefs
    std::unique_ptr<IN_TYPE[]> inv_p_coefs = std::unique_ptr<IN_TYPE[]>(new IN_TYPE[num_coefs-1]);
    for(unsigned int i=1; i < num_coefs + 1; i++)
      inv_p_coefs[i-1] = i*inv_coefs[i];
    
    // build g_inv. Not particularly efficient to create & delete all these
    // functors, but we currently lose access to inv_coefs upon creation
    mp_g_inv.reset(new HornersFunctor<IN_TYPE,num_coefs>(inv_coefs));

    // compute g_inv_prime
    m_g_inv_prime.reset(new HornersFunctor<IN_TYPE,num_coefs-1>(inv_p_coefs));

    // check if this version of g_inv is any good
    //{
    //  unsigned int N = 20;
    //  long double error_est;
    //  // base conditions
    //  if(abs(mp_g_inv(0)) > tol || abs(mp_g_inv(1) - 1) > tol){
    //    return std::numeric_limits<double>::max;
    //  }
    //  // check g_inv at N linearly spaced points
    //  for(unsigned int i=1; i<=N; i++){
    //    // check for monotonicity
    //    if((*mp_g_inv)((i-1)/(IN_TYPE)N) > (*mp_g_inv)(i/(IN_TYPE)N))
    //      return std::numeric_limits<double>::max;

    //    // estimate the error
    //    long double exact  = slow_g_inv(i/(IN_TYPE)N);
    //    long double approx = (*mp_g_inv)(i/(IN_TYPE)N);
    //    error_est += 2*(abs(exact - approx)/(exact+approx));
    //  }
    //  // return the average relative error
    //  return error_est / N;
    //}
  }

  /* Now that we have a fast approx to g_inv, we'll make it more "accurate" by
     setting our original g to g_inv_inv */
  temp_g = [this](IN_TYPE z) -> IN_TYPE {
    // FunC tables are often generated with points a bit past their right endpoint
    if(z >= 1.0)
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
      if((*m_g_inv_prime)(x) == 0.0 || NEWTON_IT > MAX_NEWTON_IT){
        // if g prime is undefined, do a hard switch to bisection
        boost::uintmax_t MAX_IT = MAX_BISECTION;
        x0 = x = toms748_solve(
            [this, &z](IN_TYPE h) -> IN_TYPE { return (*mp_g_inv)(h) - z; }, // shift g
            (IN_TYPE)0.0, (IN_TYPE)1.0, -z, (IN_TYPE)1.0 - z, eps_tolerance<IN_TYPE>(), MAX_IT).first;
      }else{
        x = x-((*mp_g_inv)(x)-z)/(*m_g_inv_prime)(x);
      }
    }while(abs(x0-x) > tol);

    return x;
  };
  mp_g.reset(new LambdaWrapper<IN_TYPE,decltype(temp_g)>(temp_g));
}

template <typename IN_TYPE>
inline std::unique_ptr<IN_TYPE[]> TransferFunctionSinh<IN_TYPE>::inverse_poly_interp(unsigned int num_coefs,
    std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp)
{
  using arma::span;
  // generate vandermonde system
  if(num_coefs < 2)
    throw std::invalid_argument("inverse_poly_coefs() must"
      " produce at least 2 polynomial coefficients");
  m_method_of_approx = std::string("inverse polynomial interpolation");
  arma::mat A = arma::ones(num_coefs,num_coefs);

  A(span(0,num_coefs-1), 1)=arma::linspace(0,1,num_coefs);
  for(unsigned int i=2; i<num_coefs; i++)
    A(span(0,num_coefs-1),i) = A(span(0,num_coefs-1),i-1) % A(span(0,num_coefs-1),1);

  // generate solution vector
  arma::vec y = arma::ones(num_coefs);
  // note: gspace returns a vector of length num_coefs whos points are all
  // between 0 and 1
  y.rows(0,num_coefs-1) = gspace(num_coefs, g, gp);

  y = arma::solve(A,y);

  if(abs(y[0]) > tol)
    throw std::runtime_error("inverse_poly_interp() failed to"
      " solve for the inverse");

  // move from arma's vector to a std::unique_ptr<double[]>
  auto coefs = std::unique_ptr<IN_TYPE[]>(new IN_TYPE[num_coefs]);
  for(unsigned int i = 0; i < num_coefs; i++)
    coefs[i] = y[i];

  return coefs;
}

/* approximate g_inv with inverse polynomial interpolation and by specifying
 * slopes at interior points as 1/g_prime */
template <typename IN_TYPE>
inline std::unique_ptr<IN_TYPE[]> TransferFunctionSinh<IN_TYPE>::inverse_poly_interior_slopes_interp(
    unsigned int num_coefs,
    std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp)
{
  using arma::span;
  // generate vandermonde system
  // assert the user is asking for an even number of coefficients
  if(num_coefs % 2 != 0)
    throw std::invalid_argument("num_coefs is odd"
      " but inverse_poly_interior_slopes_coefs() can only produce an"
      " even number of polynomial coefficients");

  m_method_of_approx = std::string("inverse polynomial interpolation specifying "
      "slopes of inner points as 1/g'");
  unsigned int M = num_coefs/2+1; // number of unique points being sampled from
  arma::mat A = arma::ones(num_coefs,num_coefs);

  A(span(0,M-1), 1)=arma::linspace(0,1,M);
  for(unsigned int i=2; i<num_coefs; i++)
    A(span(0,M-1),i) = A(span(0,M-1),i-1) % A(span(0,M-1),1);

  // set the bottom half to derivative values
  A(span(M,num_coefs-1),0) = arma::zeros(M-2,1);
  for(unsigned int i=1; i<num_coefs; i++)
    A(span(M,num_coefs-1),i) = i*A(span(1,M-2),i-1);

  // generate solution vector
  arma::vec y = arma::ones(num_coefs);
  y.rows(0,M-1) = gspace(M, g, gp);
  // using 1/g'(x_i) as an approximation for the slopes of the inverse functions
  for(int i=1; i<M-1; i++){
    y[M-1+i] = 1.0/gp(y[i]); // requires f(y[i])\neq\pm\infy
  }

  y = arma::solve(A,y);

  if(abs(y[0]) > tol)
    throw std::runtime_error("inverse_poly_interior_slopes_interp() failed to"
      " solve for the inverse");

  // move from arma's vector to a std::unique_ptr<double[]>
  auto coefs = std::unique_ptr<IN_TYPE[]>(new IN_TYPE[num_coefs]);
  for(unsigned int i = 0; i < num_coefs; i++)
    coefs[i] = y[i];

  return coefs;
}

template <typename IN_TYPE>
inline std::unique_ptr<IN_TYPE[]> TransferFunctionSinh<IN_TYPE>::inverse_polynomial_projection(
    unsigned int num_coefs,
    std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp)
{
  // TODO
  // move from arma's vector to a std::unique_ptr<double[]>
  auto coefs = std::unique_ptr<IN_TYPE[]>(new IN_TYPE[num_coefs]);
  for(unsigned int i = 0; i < num_coefs; i++)
    coefs[i] = 1;
  
  return coefs;
}

/* approximate g_inv with inverse hermite interpolation. ie specify
 * slopes at the endpoints as 1/g_prime */
template <typename IN_TYPE>
inline std::unique_ptr<IN_TYPE[]> TransferFunctionSinh<IN_TYPE>::inverse_hermite_interp(
    unsigned int num_coefs,
    std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp)
{
  using arma::span;
  // generate vandermonde system
  // assert the user is asking for at least 4 coefs
  if(num_coefs < 4)
      throw std::invalid_argument("num_coefs is too small."
      " inverse_hermite_interp() can only produce 4 or more coefs");
  m_method_of_approx = std::string("inverse polynomial interpolation specifying "
      "slopes of endpoints as 1/g'");
  unsigned int M = num_coefs-2; // number of unique points being sampled from
  arma::mat A = arma::ones(num_coefs,num_coefs);

  A(span(0,M-1), 1)=arma::linspace(0,1,M);
  for(unsigned int i=2; i<num_coefs; i++)
    A(span(0,M-1),i) = A(span(0,M-1),i-1) % A(span(0,M-1),1);

  // set the bottom half to derivative values.
  // TODO avoid the for loop
  A(span(M,num_coefs-1),0) = arma::zeros(2,1);
  for(unsigned int i=1; i<num_coefs; i++){
    A(M,i) = i*A(0,i-1);
    A(M+1,i) = i*A(M-1,i-1);
  }

  // generate solution vector
  arma::vec y = arma::ones(num_coefs);
  y.rows(0,M-1) = gspace(M, g, gp);
  // using 1/g'(x_i) as an approximation for the slopes of the inverse functions
  y[M] = 1.0/gp(y[0]); // requires f(y[i])\neq\pm\infy
  y[M+1] = 1.0/gp(y[M-1]); // requires f(y[i])\neq\pm\infy

  y = arma::solve(A,y);

  if(abs(y[0]) > tol)
    throw std::runtime_error("inverse_hermite_interp() failed to"
      " solve for the inverse");

  // move from arma's vector to a std::unique_ptr<double[]>
  auto coefs = std::unique_ptr<IN_TYPE[]>(new IN_TYPE[num_coefs]);
  for(unsigned int i = 0; i < num_coefs; i++)
    coefs[i] = y[i];
  
  return coefs;
}

template <typename IN_TYPE>
inline std::function<IN_TYPE(IN_TYPE)> TransferFunctionSinh<IN_TYPE>::newtons_inv(
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
            (IN_TYPE)0.0, (IN_TYPE)1.0, -z, (IN_TYPE)1.0 - z, eps_tolerance<IN_TYPE>(), MAX_IT).first;
      }else{
        x = x-(g(x)-z)/gp(x);
      }
    }while(abs(x0-x) > tol);

    return x;
  };
}

/* fill a vector with N linearly spaced points with respect to g in [0,1].
 * assumes g is monotone on [0,1], g(0)=0, and g(1)=1 */
template <typename IN_TYPE>
inline arma::vec TransferFunctionSinh<IN_TYPE>::gspace(unsigned int N,
    std::function<IN_TYPE(IN_TYPE)> g, std::function<IN_TYPE(IN_TYPE)> gp)
{
  using boost::math::tools::eps_tolerance;
  using boost::math::tools::toms748_solve;

  arma::vec const linear_pts = arma::linspace(0,1,N);
  arma::vec v = arma::zeros(N,1);
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
      if(gp == NULL || gp(x) == 0.0 || x < 0.0 || x > 1.0 || NEWTON_IT > MAX_NEWTON_IT){
        // if g prime is undefined, do a hard switch to bisection
        boost::uintmax_t MAX_IT = MAX_BISECTION;
        x0 = x = toms748_solve(
            [&g, &linear_pts, &i](IN_TYPE z) -> IN_TYPE { return g(z) - linear_pts[i]; }, // shift g
            (IN_TYPE)0.0, (IN_TYPE)1.0, (IN_TYPE)(-linear_pts[i]), (IN_TYPE)(1.0-linear_pts[i]), eps_tolerance<IN_TYPE>(), MAX_IT).first;
      }else{
        x = x-(g(x)-linear_pts[i])/gp(x);
      }
    }while(abs(x0-x) > tol);
    v[i]=x;
  }
  v[N-1] = 1.0;
  return v;
}

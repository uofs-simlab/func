#include <cmath>
#include <memory> // unique_ptr
#define BOOST_MATH_GAUSS_NO_COMPUTE_ON_DEMAND
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/differentiation/autodiff.hpp>
#include <boost/math/tools/toms748_solve.hpp>
#include "TransferFunctionSinh.hpp"

static double tol = 1e-8;

TransferFunctionSinh::TransferFunctionSinh(FunctionContainer *fc, double a, double b, unsigned int N) :
  TransferFunction(fc->double_func,a,b)
{
  using boost::math::quadrature::gauss_kronrod;
  using boost::math::differentiation::make_fvar;

  // load the first autodifferentiable function
  m_boost_func = fc->fvar1_func;

  // make a function for returning f's first derivative
  m_f_prime = [this](double x) -> double {
    return m_boost_func(make_fvar<double,1>(x)).derivative(1);
  };

  // build our transfer function
  // we'll be adjusting temp_g: [a,b] -> [a,b]
  std::function<double(double)> temp_g = [this](double x) -> double {
    return 1/sqrt(1 + pow(m_f_prime(x),2));
  };

  // perform adaptive quadrature with a default tol of sqrt(epsilon)
  m_c = gauss_kronrod<double, 15>::integrate(temp_g, m_minArg, m_maxArg);

  // g:[0,1] -> [0,1] integrates temp_g over all of [a,b]
  std::function<double(double)> g = [this](double x) -> double {
    if(x == 0.0) return 0;
    return gauss_kronrod<double, 15>::integrate(
        [this](double x) -> double { return 1/sqrt(1 + pow(m_f_prime(x),2)); },
        m_minArg, x*(m_maxArg-m_minArg)) / m_c;
  };

  // build g prime
  std::function<double(double)> gp = [this](double x) -> double { 
    return (m_maxArg-m_minArg)/sqrt(1 + pow(m_f_prime(m_minArg + x*(m_maxArg-m_minArg)),2)) / m_c;
  };

  // build g^{-1} by computing the inverse polynomial interpolant.
  // This is the experimental part: N (the number of sample points)
  // and the method of approximation (currently inverse poly interp
  // while specifying slopes of inner points)
  m_poly_coefs = inverse_poly_coefs(N, g, gp);
  m_numCoefs = 2*N-2;

  std::function<double(double)> g_inv = [this](double x) -> double {
      double sum = x*m_poly_coefs[m_numCoefs-1];
      for (int k=m_numCoefs-2; k>0; k--)
        sum = x*(m_poly_coefs[k] + sum);
      return sum;
  };

  mp_g_and_g_inv = std::make_pair(g, g_inv);
}

arma::vec TransferFunctionSinh::inverse_poly_coefs(unsigned int N, 
    std::function<double(double)> g, std::function<double(double)> gp)
{
  using arma::span;
  //generate vandermonde system
  unsigned int M = 2*N-2;
  arma::mat A = arma::ones(M,M);

  A(span(0,N-1), 1)=arma::linspace(0,1,N);
  for(unsigned int i=2; i<M; i++)
    A(span(0,N),i) = A(span(0,N),i-1) % A(span(0,N),1);

  // set the bottom half to derivative values
  A(span(N,M-1),0) = arma::zeros(N-2,1);
  for(unsigned int i=1; i<M; i++)
    A(span(N,M-1),i) = i*A(span(1,N-2),i-1);

  // generate solution vector
  arma::vec y = arma::ones(M);
  y.rows(0,N-1) = gspace(N, g, gp);
  for(int i=1; i<N-1; i++){
    y[N-1+i] = 1.0/gp(y[i]); // requires f(y[i])\neq\pm\infy
  }

  arma::vec poly_coefs = arma::solve(A,y);

  // evaluate the resulting polynomial using horners method
  return poly_coefs;
}

/* fill a vector with N linearly spaced points with respect to g in [0,1].
 * assumes g is monotone on [0,1], g(0)=0, and g(1)=1 */
arma::vec TransferFunctionSinh::gspace(unsigned int N, std::function<double(double)> g,
    std::function<double(double)> gp)
{
  using boost::math::tools::eps_tolerance;
  using boost::math::tools::toms748_solve;

  arma::vec const linear_pts = arma::linspace(0,1,N);
  arma::vec v = arma::zeros(N,1);
  for(int i=1; i<N-1; i++){
    // Find what g maps to linear_pts[i]
    double x0;
    double x = linear_pts[i];
    // Use Newton's method if we can get away with it. o.w. resort to bisection
    do{
      x0 = x;
      if(gp == NULL || gp(x) == 0.0 || x < 0.0 || x > 1.0){
        // if g prime is undefined, do a hard switch to bisection
        boost::uintmax_t const C_MAX_IT = 54;
        boost::uintmax_t MAX_IT = C_MAX_IT;
        x0 = x = toms748_solve(
            [&g, &linear_pts, &i](double z) -> double { return g(z) - linear_pts[i]; }, // shift g
            0.0, 1.0, -linear_pts[i], 1.0-linear_pts[i], eps_tolerance<double>(), MAX_IT).first;
      }else{
        x = x-(g(x)-linear_pts[i])/gp(x);
      }
    }while(abs(x0-x) > tol);
    v[i]=x;
  }
  v[N-1] = 1;
  return v;
}

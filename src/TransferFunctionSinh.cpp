#include <cmath>
#include <memory> // unique_ptr
#include <iostream>
#define BOOST_MATH_GAUSS_NO_COMPUTE_ON_DEMAND
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/differentiation/autodiff.hpp>
#include <boost/math/tools/toms748_solve.hpp>
#include "TransferFunctionSinh.hpp"
#include "FunctionContainer.hpp"

static double tol = 1e-8;

template <typename IN_TYPE, typename OUT_TYPE>
TransferFunctionSinh<IN_TYPE,OUT_TYPE>::TransferFunctionSinh(FunctionContainer<IN_TYPE,OUT_TYPE> *fc, IN_TYPE a, IN_TYPE b) :
  TransferFunction<IN_TYPE>(a,b)
{
  using boost::math::quadrature::gauss_kronrod;
  using boost::math::differentiation::make_fvar;

  // load the first autodifferentiable function
  m_boost_func = fc->autodiff1_func;

  // build a function to return the first derivative of f
  f_prime = [this](IN_TYPE x) -> IN_TYPE {
    return m_boost_func(make_fvar<IN_TYPE,1>(x)).derivative(1);
  };

  // build our transfer function
  // we'll be adjusting temp_g: [a,b] -> [a,b]
  std::function<IN_TYPE(IN_TYPE)> temp_g = [this](IN_TYPE x) -> IN_TYPE {
    return 1/sqrt(1 + pow(f_prime(x),2));
  };

  // perform adaptive quadrature with a default tol of sqrt(epsilon)
  m_scale_factor = gauss_kronrod<IN_TYPE, 15>::integrate(temp_g, this->m_minArg, this->m_maxArg);

  // g:[0,1] -> [0,1] integrates temp_g over all of [a,b] and scales the answer
  // such that g(1) = 1
  this->g = [this](IN_TYPE x) -> IN_TYPE {
    if(x == 0.0) return 0.0; // boost gets upset if we don't do this
    return gauss_kronrod<IN_TYPE, 15>::integrate(
        [this](IN_TYPE t) -> IN_TYPE { return 1/sqrt(1 + pow(f_prime(t),2)) / m_scale_factor; },
        this->m_minArg, this->m_minArg+x*(this->m_maxArg-this->m_minArg));
  };

  // build g prime
  g_prime = [this](IN_TYPE x) -> IN_TYPE { 
    return (this->m_maxArg-this->m_minArg) /
      sqrt(1 + pow(f_prime(this->m_minArg + x*(this->m_maxArg-this->m_minArg)),2)) / m_scale_factor;
  };

  /* build g^{-1} by computing the inverse polynomial interpolant.
   * This is the experimental part, and so it has been made modular.
   * Want g^{-1} evals to be quick so we currently approx it with
   * some type of inverse polynomial interpolation. */
  this->g_inv = newtons_g_inv();
  // g_inv = inverse_poly_interp(4, g, g_prime);
  // g_inv = inverse_poly_interior_slopes_interp(8, g, g_prime);
  // g_inv = inverse_hermite_interp(5, g, g_prime);
  //
  // TODO cycle through ways to approx. g_inv until we find the best fit
}

/* approximate g_inv with inverse polynomial interpolantion */
template <typename IN_TYPE, typename OUT_TYPE>
std::function<double(double)> TransferFunctionSinh<IN_TYPE,OUT_TYPE>::inverse_poly_interp(unsigned int num_coefs,
    std::function<double(double)> g, std::function<double(double)> gp)
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
  m_num_coefs = y.n_rows;
  m_polynomial_coefs.reset(new OUT_TYPE[m_num_coefs]);
  for(unsigned int i = 0; i < m_num_coefs; i++)
    m_polynomial_coefs[i] = y[i];

  // evaluate the resulting polynomial using horners
  return [this](double z) -> double {
      double sum = z*m_polynomial_coefs[m_num_coefs-1];
      for (int k=m_num_coefs-2; k>0; k--)
        sum = z*(m_polynomial_coefs[k] + sum);
      return sum;
  };
}

/* approximate g_inv with inverse polynomial interpolation and by specifying
 * slopes at interior points as 1/g_prime */
template <typename IN_TYPE, typename OUT_TYPE>
std::function<double(double)> TransferFunctionSinh<IN_TYPE,OUT_TYPE>::inverse_poly_interior_slopes_interp(
    unsigned int num_coefs,
    std::function<double(double)> g, std::function<double(double)> gp)
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
  m_num_coefs = y.n_rows;
  m_polynomial_coefs.reset(new OUT_TYPE[m_num_coefs]);
  for(unsigned int i = 0; i < m_num_coefs; i++)
    m_polynomial_coefs[i] = y[i];

  // evaluate the resulting polynomial using horners
  return [this](double z) -> double {
      double sum = z*m_polynomial_coefs[m_num_coefs-1];
      for (int k=m_num_coefs-2; k>0; k--)
        sum = z*(m_polynomial_coefs[k] + sum);
      return sum;
  };
}

/* approximate g_inv with inverse hermite interpolation. ie specify
 * slopes at the endpoints as 1/g_prime */
template <typename IN_TYPE, typename OUT_TYPE>
std::function<double(double)> TransferFunctionSinh<IN_TYPE,OUT_TYPE>::inverse_hermite_interp(
    unsigned int num_coefs,
    std::function<double(double)> g, std::function<double(double)> gp)
{
  using arma::span;
  // generate vandermonde system
  // assert the user is asking for at least 4 coefs
  if(num_coefs < 4)
      throw std::invalid_argument("num_coefs is too small"
      " or it is odd but inverse_poly_interior_slopes_coefs() can only"
      " produce an even number of polynomial coefficients");
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
  m_num_coefs = y.n_rows;
  m_polynomial_coefs.reset(new OUT_TYPE[m_num_coefs]);
  for(unsigned int i = 0; i < m_num_coefs; i++)
    m_polynomial_coefs[i] = y[i];

  // evaluate the resulting polynomial using horners
  return [this](double z) -> double {
      double sum = z*m_polynomial_coefs[m_num_coefs-1];
      for (int k=m_num_coefs-2; k>0; k--)
        sum = z*(m_polynomial_coefs[k] + sum);
      return sum;
  };
}

/* verrry expensive version of g^{-1} that is just Newton's method 
 * or bisection if things go south */
template <typename IN_TYPE, typename OUT_TYPE>
std::function<double(double)> TransferFunctionSinh<IN_TYPE,OUT_TYPE>::newtons_g_inv()
{
  m_method_of_approx = std::string("inverse polynomial interpolation");
  return [this](double z) -> double {
    using boost::math::tools::eps_tolerance;
    using boost::math::tools::toms748_solve;

    const unsigned int MAX_NEWTON_IT = 20;
    boost::uintmax_t const MAX_BISECTION = 54;

    // Find what g maps to x
    double x = z;
    double x0;
    unsigned int NEWTON_IT = 0;
    // Use Newton's method if we can get away with it. o.w. resort to bisection
    do{
      NEWTON_IT += 1;
      x0 = x;
      //if(g_prime == NULL || g_prime(x) == 0.0 || x < 0.0 || x > 1.0 || NEWTON_IT > MAX_NEWTON_IT){
      if(g_prime == NULL || g_prime(x) == 0.0 || NEWTON_IT > MAX_NEWTON_IT){
        // if g prime is undefined, do a hard switch to bisection
        boost::uintmax_t MAX_IT = MAX_BISECTION;
        x0 = x = toms748_solve(
            [this, &z](double h) -> double { return this->g(h) - z; }, // shift g
            0.0, 1.0, -z, 1.0 - z, eps_tolerance<double>(), MAX_IT).first;
      }else{
        x = x-(this->g(x)-z)/g_prime(x);
      }
    }while(abs(x0-x) > tol);

    return x;
  };
}

/* fill a vector with N linearly spaced points with respect to g in [0,1].
 * assumes g is monotone on [0,1], g(0)=0, and g(1)=1 */
template <typename IN_TYPE, typename OUT_TYPE>
arma::vec TransferFunctionSinh<IN_TYPE,OUT_TYPE>::gspace(unsigned int N,
    std::function<double(double)> g, std::function<double(double)> gp)
{
  using boost::math::tools::eps_tolerance;
  using boost::math::tools::toms748_solve;

  arma::vec const linear_pts = arma::linspace(0,1,N);
  arma::vec v = arma::zeros(N,1);
  const unsigned int MAX_NEWTON_IT = 20;
  boost::uintmax_t const MAX_BISECTION = 54;

  for(int i=1; i<N-1; i++){
    // Find what g maps to linear_pts[i]
    double x0;
    double x = linear_pts[i];
    unsigned int NEWTON_IT = 0;
    // Use Newton's method if we can get away with it. o.w. resort to bisection
    do{
      NEWTON_IT += 1;
      x0 = x;
      if(gp == NULL || gp(x) == 0.0 || x < 0.0 || x > 1.0 || NEWTON_IT > MAX_NEWTON_IT){
        // if g prime is undefined, do a hard switch to bisection
        boost::uintmax_t MAX_IT = MAX_BISECTION;
        x0 = x = toms748_solve(
            [&g, &linear_pts, &i](double z) -> double { return g(z) - linear_pts[i]; }, // shift g
            0.0, 1.0, -linear_pts[i], 1.0-linear_pts[i], eps_tolerance<double>(), MAX_IT).first;
      }else{
        x = x-(g(x)-linear_pts[i])/gp(x);
      }
    }while(abs(x0-x) > tol);
    v[i]=x;
  }
  v[N-1] = 1.0;
  return v;
}

template <typename IN_TYPE, typename OUT_TYPE>
void TransferFunctionSinh<IN_TYPE,OUT_TYPE>::print_details(std::ostream& out)
{
  out << "arcsinh transfer function approximating g_inv with ";
  out << m_method_of_approx;
}

template <typename IN_TYPE, typename OUT_TYPE>
void TransferFunctionSinh<IN_TYPE,OUT_TYPE>::print_debugging_details(std::ostream& out)
{
  print_details(out);
  out << std::endl;
  for(int i=0; i<m_num_coefs; i++)
    out << m_polynomial_coefs[i] << std::endl;
}

template class TransferFunctionSinh<double>;
template class TransferFunctionSinh<double,float>;
template class TransferFunctionSinh<float>;
template class TransferFunctionSinh<float,double>;

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
  // make a function for returning f's first derivative
  std::function<double(double)> fp = [&fc](double x) -> double {
    return fc->fvar1_func(make_fvar<double,1>(x)).derivative(1);
  };

  // build our transfer function
  std::function<double(double)> temp_g = [&fp](double x) -> double { 
    return 1/sqrt(1 + pow(fp(x),2));
  };

  // scale g such that g(1) = 1
  // perform adaptive quadrature with a default tol of sqrt(epsilon)
  double c = gauss_kronrod<double, 15>::integrate(temp_g, 0, 1);
  std::function<double(double)> g = [&c,&temp_g](double x) -> double {
    return gauss_kronrod<double, 15>::integrate(temp_g, 0, x)/c;
  };

  // build g prime
  std::function<double(double)> gp = [&a, &b, &c, &fp](double x) -> double { 
    return c*(b-a)/sqrt(1 + pow(fp(a + x*(b-a)),2));
  };

  // build g^{-1} by computing the inverse polynomial interpolant.
  // This is the experimental part: N (the number of sample points)
  // and the method of approximation (currently inverse poly interp
  // while specifying slopes of inner points)
  std::function<double(double)> g_inv = inverse_poly_coefs(5, g, gp);
  mp_g_and_g_inv = std::make_pair(g, g_inv);
}

std::function<double(double)> TransferFunctionSinh::inverse_poly_coefs(unsigned int N, 
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
  return [&poly_coefs, &N](double x) -> double {
      double sum = x*poly_coefs[N];
      for (int k=N-1; k>0; k--)
        sum = x*(poly_coefs[k] + sum);
      return poly_coefs[0]+sum;
  };
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
      if(gp == NULL || gp(x) == 0){
        // if g prime is undefined, do a hard switch to bisection
        boost::uintmax_t MAX_IT = 52;
        x = linear_pts[i];
        x0 = x = toms748_solve(
            [&g, &x](double z) -> double { return g(z) - x; }, // shift g
            0.0, 1.0, -x, 1.0-x, eps_tolerance<double>(), MAX_IT).first;
      }else{
        x = x-(g(x)-linear_pts[i])/gp(x);
      }
    }while(abs(x0-x) > tol);
    v[i]=x;
  }
  v[N-1] = 1;
  return v;
}

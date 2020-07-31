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
  - TODO need to be aware of the possibility of g' being so bad that it causes
  table evaluations to be out of bounds resulting in a seg fault
  - TODO If an implementation of g_inv is bad then swap it for something else.
  if they're all bad, alter the original g.
  - TODO decide on pros/cons of certain options.
  - TODO decide if helper functions should be in base class. Currently their
  fairly general
 */
#pragma once
#include <memory>
#include <string>
#include <functional>
#define ARMA_USE_CXX11
#include <armadillo>
#include "TransferFunction.hpp"
#include "FunctionContainer.hpp"

template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE>
class TransferFunctionSinh : public TransferFunction<IN_TYPE>
{
  /* --- More member vars --- */
  /* most of these functions will be lambdas that capture this because we
   * want to pass these functions as args to other functions, and we can't
   * make them static */
  std::function<fvar<OUT_TYPE,1>(fvar<IN_TYPE,1>)> m_boost_func;
  std::function<OUT_TYPE(IN_TYPE)> f_prime;
  IN_TYPE m_scale_factor; // scaling factor for g
  std::function<IN_TYPE(IN_TYPE)> g_prime;

  std::string m_method_of_approx;

  unsigned int m_num_coefs;
  __attribute__((aligned)) std::unique_ptr<OUT_TYPE[]> m_polynomial_coefs; 

  /* --- Private Functions for approximating g^{-1} --- */
  /* Approximate g^{-1} using inverse polynomial interpolation */
  std::function<double(double)> inverse_poly_interp(unsigned int num_coefs,
      std::function<double(double)> g, std::function<double(double)> gp);

  /* Approximate g^{-1} using inverse polynomial interpolation
   * and specifying the slopes of inner points */
  std::function<double(double)> inverse_poly_interior_slopes_interp(
      unsigned int num_coefs, std::function<double(double)> g,
      std::function<double(double)> gp);

  /* Approximate g^{-1} using inverse polynomial interpolation
   * and specifying slopes at function endpoints */
  std::function<double(double)> inverse_hermite_interp(unsigned int num_coefs,
      std::function<double(double)> g, std::function<double(double)> gp);

  /* Approx. g^{-1} using newton's method. v slow but accurate */
  std::function<double(double)> newtons_g_inv();

  /* --- Private Helper functions --- */
  /* fill a vector with N linearly spaced points with respect to g in [0,1].
   * assumes g is monotone on [0,1], g(0)=0, and g(1)=1 */
  arma::vec gspace(unsigned int N, std::function<double(double)> g,
      std::function<double(double)> gp = NULL);

public:
  /* public constructor */
  TransferFunctionSinh(FunctionContainer<IN_TYPE,OUT_TYPE> *fc, IN_TYPE a, IN_TYPE b);
  void print_details(std::ostream&) override;
  void print_debugging_details(std::ostream&);
};

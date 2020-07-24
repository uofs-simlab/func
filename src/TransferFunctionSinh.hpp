/*
  builds a pair of functions, g and g^{-1} from a function f and its domain [a,b].
  To keep evaluations of g^{-1} quick, it's being approximated using inverse polynomial
  interpolation.

  Somewhat experimental atm.
  TODO decide on pros/cons of certain options
  TODO decide if helper functions should be in base class
 */
#pragma once
#define ARMA_USE_CXX11
#include <armadillo>
#include <functional>
#include "TransferFunction.hpp"
#include "FunctionContainer.hpp"

class TransferFunctionSinh : public TransferFunction
{
  /* --- More member vars --- */
  std::function<double(double)> m_f_prime;
  std::function<fvar1(fvar1)> m_boost_func;
  unsigned int m_numCoefs;
  arma::vec m_poly_coefs; // TODO make this an __attribute__((aligned)) std::unique_ptr<double[]>
  double m_c; // scaling factor for g

  /* --- Private Functions for approximating g^{-1} --- */
  /* Approximate g^{-1} using inverse polynomial interpolation */
  arma::vec inverse_poly_coefs(unsigned int N, 
      std::function<double(double)> g, std::function<double(double)> gp);

  /* --- Private Helper functions --- */
  /* fill a vector with N linearly spaced points with respect to g in [0,1].
   * assumes g is monotone on [0,1], g(0)=0, and g(1)=1 */
  arma::vec gspace(unsigned int N, std::function<double(double)> g,
      std::function<double(double)> gp = NULL);

public:
  /* public constructor & destructor */
  TransferFunctionSinh(FunctionContainer *fc, double a, double b, unsigned int N);
};

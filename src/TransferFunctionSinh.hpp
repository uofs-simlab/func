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
  /* --- Private Functions for approximating g^{-1} --- */
  /* Approximate g^{-1} using inverse polynomial interpolation */
  std::function<double(double)> inverse_poly_coefs(unsigned int N, 
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

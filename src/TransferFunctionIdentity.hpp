/* 
   Simplest possible implementation of a TransferFunctionInterface.
   This is only used when the user does not have a version of Boost
   with automatic differentiation (any version older than 1.71.0)
   This will create a uniform grid and as such, can only slow down
   table evaluation.
*/
#pragma once
#include "TransferFunctionInterface.hpp"
#include "FunctionContainer.hpp"

template <typename IN_TYPE>
class TransferFunctionIdentity final : public TransferFunctionInterface<IN_TYPE>
{
  IMPLEMENT_TRANSFER_FUNCTION_INTERFACE(IN_TYPE);

public:
  /* public constructor */
  template<typename OUT_TYPE>
  TransferFunctionIdentity(FunctionContainer<IN_TYPE,OUT_TYPE> *fc,
      LookupTableParameters<IN_TYPE> par) :
    TransferFunctionInterface<IN_TYPE>(fc,par) {}

  IN_TYPE g(IN_TYPE x) override { return x; }
  // table hash is baked into g_inv
  IN_TYPE g_inv(IN_TYPE x) override { return (x-m_minArg)/m_stepSize; }

  void print_details(std::ostream& out) override
  {
    out << m_approx_method << " " << NUM_COEFS;
  }
};

/*
  Intermediate abstract class for LUTs with grid points spaced such that
  error is equally distributed throughout the interval [a,b]
*/
#pragma once
#include "FunctionContainer.hpp"
#include "UniformLookupTable.hpp"

struct NonUniformLookupTableParameters : public UniformLookupTableParameters
{
  // add a transfer function as a parameter
};

class NonUniformLookupTable : public UniformLookupTable
{
protected:
  TransferFunction *m_transferFunction;
  std::function<double(double)> m_g, m_g_inv;
public:
  // set the transfer function
  NonUniformLookupTable(FunctionContainer *func_container, UniformLookupTableParameters par) :
    UniformLookupTable(func_container, par), m_transferFunction(par.transferFunction),
    m_g(m_transferFunction->function_pair().first), m_g_inv(m_transferFunction->function_pair().second) {}
  virtual ~NonUniformLookupTable(){};
};

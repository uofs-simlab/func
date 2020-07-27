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
  std::shared_ptr<TransferFunction> m_transferFunction;
public:
  // set the transfer function
  NonUniformLookupTable(FunctionContainer *func_container, UniformLookupTableParameters par) :
    UniformLookupTable(func_container, par), m_transferFunction(par.transferFunction) {}
  virtual ~NonUniformLookupTable(){};
};

/*
  Intermediate abstract class for LUTs with grid points spaced such that
  error is equally distributed throughout the interval [a,b]
*/
#pragma once
#include "FunctionContainer.hpp"
#include "UniformLookupTable.hpp"
#include <exception>

class NonUniformLookupTable : public UniformLookupTable
{
protected:
  std::shared_ptr<TransferFunction> m_transferFunction;
public:
  // set the transfer function
  NonUniformLookupTable(FunctionContainer *func_container, UniformLookupTableParameters par) :
    UniformLookupTable(func_container, par), m_transferFunction(par.transferFunction)
  {
    if(m_transferFunction == nullptr)
      throw std::invalid_argument("NonUniformLookupTable needs a transfer function");
  }
  virtual ~NonUniformLookupTable(){};
};

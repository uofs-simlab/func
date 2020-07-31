/*
  Intermediate abstract class for LUTs with grid points spaced such that
  error is equally distributed throughout the interval [a,b]
*/
#pragma once
#include "FunctionContainer.hpp"
#include "UniformLookupTable.hpp"
#include <exception>

template <typename IN_TYPE, typename OUT_TYPE>
class NonUniformLookupTable : public UniformLookupTable<IN_TYPE,OUT_TYPE>
{
protected:
  std::shared_ptr<TransferFunction<IN_TYPE>> m_transferFunction;
public:
  // set the transfer function
  NonUniformLookupTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, UniformLookupTableParameters<IN_TYPE> par) :
    UniformLookupTable(func_container, par), m_transferFunction(par.transferFunction)
  {
    if(m_transferFunction == nullptr)
      throw std::invalid_argument("NonUniformLookupTable needs a transfer function");
  }
  virtual ~NonUniformLookupTable(){};
};

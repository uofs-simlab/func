/*
  An extension of the abstract class UniformLookupTable. Tables that
  inherit from this class will be equiped with a transfer function which
  will allow for nonuniform grid spacing.
  Ideally, the transfer function is constructed such that table error
  is equidistributed throughout the interval [a,b]
  Currently not an especially useful intermediate class.
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

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

#define INHERIT_NONUNIFORM_LUT(IN_TYPE,OUT_TYPE) \
  using NonUniformLookupTable<IN_TYPE,OUT_TYPE>::m_transferFunction

template <typename IN_TYPE, typename OUT_TYPE>
class NonUniformLookupTable : public UniformLookupTable<IN_TYPE,OUT_TYPE>
{
protected:
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE);
  std::shared_ptr<TransferFunction<IN_TYPE>> m_transferFunction;
public:
  // set the transfer function
  NonUniformLookupTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, UniformLookupTableParameters<IN_TYPE> par) :
    UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, par), m_transferFunction(par.transferFunction)
  {
    if(m_transferFunction == NULL)
      throw std::invalid_argument("NonUniformLookupTable needs a transfer function");
  }
  virtual ~NonUniformLookupTable(){};
};

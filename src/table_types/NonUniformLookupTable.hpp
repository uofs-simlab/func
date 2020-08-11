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
#include "TransferFunctionSinh.hpp"
#include <exception>

#define INHERIT_NONUNIFORM_LUT(IN_TYPE,OUT_TYPE) \
  using NonUniformLookupTable<IN_TYPE,OUT_TYPE,TRANSFER_FUNC_TYPE>::m_transferFunction

template <typename IN_TYPE, typename OUT_TYPE, class TRANSFER_FUNC_TYPE>
class NonUniformLookupTable : public UniformLookupTable<IN_TYPE,OUT_TYPE>
{
protected:
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE);
  std::unique_ptr<TRANSFER_FUNC_TYPE> m_transferFunction;
public:
  // set the transfer function
  NonUniformLookupTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, UniformLookupTableParameters<IN_TYPE> par) :
    UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, par)
  {
    // TODO check TRANSFER_FUNC_TYPE actually inherits from TransferFunctionInterface
    // TODO add more sophisticated transfer function generation here
    m_transferFunction = std::unique_ptr<TRANSFER_FUNC_TYPE>(new TRANSFER_FUNC_TYPE(func_container,par.minArg,par.maxArg,par.stepSize));
  }
  virtual ~NonUniformLookupTable(){};
};

// TODO register based on unsigned int?
#define REGISTER_NONUNIFORM_IMPL(classname,IN_TYPE,OUT_TYPE,TRANSFER_FUNC) \
  template<> const \
    UniformLookupTableRegistrar<classname<IN_TYPE,OUT_TYPE,TRANSFER_FUNC>,IN_TYPE,OUT_TYPE> \
    classname<IN_TYPE,OUT_TYPE,TRANSFER_FUNC>::registrar(STR(classname))

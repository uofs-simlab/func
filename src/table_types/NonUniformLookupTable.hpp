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
  TRANSFER_FUNC_TYPE m_transferFunction;
  virtual OUT_TYPE get_table_entry(unsigned int i, unsigned int j) override=0;
  virtual unsigned int get_num_coefs() override=0;

public:
  /* set the transfer function */
  NonUniformLookupTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, UniformLookupTableParameters<IN_TYPE> par) :
    UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, par),
    m_transferFunction(TRANSFER_FUNC_TYPE(func_container,par.minArg,par.maxArg,par.stepSize))
  {
    // TODO check if TRANSFER_FUNC_TYPE actually inherits from TransferFunctionInterface if the compiler doesn't already
  }

  /* Build from file */
  NonUniformLookupTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, std::string filename) :
    UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, filename),
    m_transferFunction(TRANSFER_FUNC_TYPE(func_container,m_minArg,m_maxArg,m_stepSize))
  {
    // recompute m_grid
    m_grid.reset(new IN_TYPE[m_numIntervals]);
    for (int ii=0; ii < m_numIntervals; ++ii)
      m_grid[ii] = m_transferFunction.g(m_minArg + ii*m_stepSize);
  }

  std::pair<IN_TYPE,IN_TYPE> arg_bounds_of_interval(unsigned intervalNumber) override
  {
    return std::make_pair(m_grid[intervalNumber], m_grid[intervalNumber+1]);
  }

  void print_details(std::ostream &out) override
  {
    out << m_name << "<";
    m_transferFunction.print_details(out);
    out << "> " << m_minArg << " " << m_maxArg << " "
        << m_stepSize << " " << m_numIntervals << " ";
  }

  virtual ~NonUniformLookupTable(){};
};

// We'll register TRANSFER_FUNC_TYPE as TransferFunctionSinh and each key will seemingly be
// templated on an unsigned int. That way we still have flexible NonUniform LUTs but it's not
// so stiff to write out
#define FUNC_REGISTER_NONUNIFORM_IMPL(classname,IN_TYPE,OUT_TYPE,N) \
  template<> const \
    UniformLookupTableRegistrar<classname<IN_TYPE,OUT_TYPE,TransferFunctionSinh<IN_TYPE,N>>,IN_TYPE,OUT_TYPE> \
    classname<IN_TYPE,OUT_TYPE,TransferFunctionSinh<IN_TYPE,N>>::registrar(FUNC_STR(classname<N>))

#define FUNC_REGISTER_TEMPLATE_NONUNIFORM_IMPL(classname,IN_TYPE,OUT_TYPE,N,other...) \
  template<> const \
    UniformLookupTableRegistrar<classname<IN_TYPE,OUT_TYPE,TransferFunctionSinh<IN_TYPE,N>>,IN_TYPE,OUT_TYPE> \
    classname<IN_TYPE,OUT_TYPE,TransferFunctionSinh<IN_TYPE,N>,other>::registrar(FUNC_STR(classname<N,other>))

#define FUNC_REGISTER_DOUBLE_AND_FLOAT_NONUNIFORM_IMPL(classname,N) \
  FUNC_REGISTER_NONUNIFORM_IMPL(classname,double,double,N); \
  FUNC_REGISTER_NONUNIFORM_IMPL(classname,double,float,N);  \
  FUNC_REGISTER_NONUNIFORM_IMPL(classname,float,double,N);  \
  FUNC_REGISTER_NONUNIFORM_IMPL(classname,float,float,N);

#define FUNC_REGISTER_TEMPLATED_NONUNIFORM_IMPL(classname,N,other...) \
  FUNC_REGISTER_NONUNIFORM_IMPL(classname,double,double,N); \
  FUNC_REGISTER_NONUNIFORM_IMPL(classname,double,float,N);  \
  FUNC_REGISTER_NONUNIFORM_IMPL(classname,float,double,N);  \
  FUNC_REGISTER_NONUNIFORM_IMPL(classname,float,float,N);

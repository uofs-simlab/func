/*
  Abstract base class for implementations.

  Derived classes generally do the following:
  - set up anything needed for 'evaluation'
  - determine size of data needed for 'evaluation'
  - override the brackets operator to perform 'evaluation'
  - cleanup in destructor
*/
#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <functional>

/* macro to get the EvaluationImplementation's member variables
   without having to sprinkle "this->" throughout our code.
   Remember to at least make these "using" statments protected */
#define INHERIT_EVALUATION_IMPL(TIN,TOUT) \
  using EvaluationImplementation<TIN,TOUT>::m_func; \
  using EvaluationImplementation<TIN,TOUT>::m_order; \
  using EvaluationImplementation<TIN,TOUT>::m_name; \
  using EvaluationImplementation<TIN,TOUT>::m_dataSize; \
  using EvaluationImplementation<TIN,TOUT>::m_minArg; \
  using EvaluationImplementation<TIN,TOUT>::m_maxArg

namespace func {

template <typename TIN, typename TOUT = TIN>
class EvaluationImplementation
{
protected:

  std::string  m_name;     // name of implementation type

  std::function<TOUT(TIN)>   m_func; // mathematical function to evaluate

  TIN      m_minArg, m_maxArg; // bounds of evaluation

  unsigned     m_order;    // order of accuracy of implementation
  unsigned     m_dataSize; // size of relevant data for impl evaluation

public:

  // Every class inheriting from this one use a FunctionContainer as
  // their first arg (aside from UniformFailureProofTable).
  EvaluationImplementation(std::function<TOUT(TIN)> func = nullptr, std::string name = "") :
    m_name(name), m_func(func), m_minArg(0), m_maxArg(0) {}

  virtual ~EvaluationImplementation(){};
  EvaluationImplementation(EvaluationImplementation&&) = default;

  virtual TOUT operator()(TIN x) = 0;
  virtual void print_details(std::ostream& out)
  {
    out << m_minArg << " " << m_maxArg << " ";
  };
  virtual void print_details_json(std::ostream& out)=0;

  /* public access of protected data */
  virtual TIN min_arg() { return m_minArg; }
  virtual TIN max_arg() { return m_maxArg; }
  unsigned int order() const { return m_order; }
  unsigned int size() const { return m_dataSize; }
  std::string name() const { return m_name; }
  std::function<TOUT(TIN)> function() const { return m_func; }
};
}

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

template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE>
class EvaluationImplementation
{
protected:

  std::function<OUT_TYPE(IN_TYPE)>   mp_func; // mathematical function to evaluate

  IN_TYPE            m_minArg, m_maxArg; // bounds of evaluation
  unsigned           m_order;    // order of accuracy of implementation
  std::string        m_name;     // name of implementation type
  unsigned           m_dataSize; // size of relevant data for impl evaluation

public:

  // Every class inheriting from this one use a FunctionContainer as 
  // their first arg (aside from UniformFailureProofTable).
  EvaluationImplementation(std::function<OUT_TYPE(IN_TYPE)> func = NULL, std::string name = "") :
    m_name(name), mp_func(func), m_minArg(0), m_maxArg(0) {}

  virtual ~EvaluationImplementation(){};

  virtual OUT_TYPE operator()(IN_TYPE x) = 0;
  virtual void print_details(std::ostream& out)
  {
    out << m_minArg << " " << m_maxArg << " ";
  };

  /* public access of protected data */
  IN_TYPE min_arg(){ return m_minArg; };
  IN_TYPE max_arg(){ return m_maxArg; };
  unsigned order(){ return m_order; };
  unsigned size(){ return m_dataSize; };
  std::string name(){ return m_name; };
  std::function<OUT_TYPE(IN_TYPE)> function(){ return mp_func; };
};

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
#include <iostream>
#include <functional>

class EvaluationImplementation
{
protected:

  std::function<double(double)> mp_func;   // mathematical function to evaluate
  double             m_minArg, m_maxArg; // bounds of evaluation

  unsigned           m_order;    // order of accuracy of implementation
  std::string        m_name;     // name of implementation type

  unsigned           m_dataSize; // size of relevant data for impl evaluation

public:

  // Every class inheriting from this one use a FunctionContainer as their first arg (aside from UniformFailureProofTable)
  EvaluationImplementation(std::function<double(double)> func = NULL, std::string name = "");
  virtual ~EvaluationImplementation(){};

  virtual double operator()(double x) = 0;
  virtual void print_details(std::ostream& out)
  {
    out << m_minArg << " " << m_maxArg << " ";
  };

  /* public access of protected data */
  double min_arg();
  double max_arg();
  unsigned order();
  unsigned size();
  std::string name();
  std::function<double(double)> function();
};

/*
  Abstract base class for implementations.

  Derived classes generally do the following:
  - set up anything needed for 'evaluation'
  - determine size of data needed for 'evaluation'
  - override the brackets operator to perform 'evaluation'
  - cleanup in destructor
  - Can be constructed with an optional vector of special points 
  to specify discontinuities
*/
#pragma once
#include "SpecialPoint.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <functional>

class EvaluationImplementation
{
protected:

  std::function<double(double)>   mp_func; // mathematical function to evaluate
  std::vector<SpecialPoint>       m_special_points; // specify any discontinuities

  double             m_minArg, m_maxArg; // bounds of evaluation
  unsigned           m_order;    // order of accuracy of implementation
  std::string        m_name;     // name of implementation type
  unsigned           m_dataSize; // size of relevant data for impl evaluation

public:

  // Every class inheriting from this one use a FunctionContainer as 
  // their first arg (aside from UniformFailureProofTable and CompoundLookupTable).
  EvaluationImplementation(std::function<double(double)> func = NULL, std::string name = "", 
      std::vector<SpecialPoint> points = std::vector<SpecialPoint>());

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
  std::vector<SpecialPoint> special_points();
};

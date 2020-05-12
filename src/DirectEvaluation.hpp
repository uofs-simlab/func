/*
  Direct evaluation of a function

  Usage example:
    DirectEvaluation de(&function,0,10);
    double f = de(0.87354);

  Notes:
  - basically just a wrapper around a function
  - requires a max and min value (for consistency with EvaluationImplementation)
*/
#pragma once

/* A macro for recording where the given function is being evaluated.  */
#ifdef FUNC_RECORD
  #define FUNC_RECORD_ARG(X) (record_arg(X))
#else
  #define FUNC_RECORD_ARG(X)
#endif

#include "EvaluationImplementation.hpp"
#include <cmath>

class DirectEvaluation final : public EvaluationImplementation
{

private:
  // the fields methods needed to record where the function is being evaluated
  unsigned int *m_histogram;
  unsigned int m_histSize; 
  double m_peak;
  unsigned int calcIndex(double x);
  void record_arg(double x);
  std::string histogramToString();
 
public:
  DirectEvaluation(EvaluationFunctor<double,double> *func, double min = 0, double max = 1, unsigned int histSize = 10);
  double operator()(double x) override;
  void print_details(std::ostream& out) override;
  void print_details_json(std::ostream& out);
  ~DirectEvaluation();
};

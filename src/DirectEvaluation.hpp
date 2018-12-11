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
#include "EvaluationImplementation.hpp"
#include <cmath>

class DirectEvaluation final : public EvaluationImplementation
{
public:
  DirectEvaluation(EvaluationFunctor<double,double> *func, double min = 0, double max = 1);
  double operator()(double x) override;
  ~DirectEvaluation();
  void print_details(std::ostream& out) override;
};

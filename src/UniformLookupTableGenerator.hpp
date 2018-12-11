/*
  Generate an UniformLookupTable of a table type T (template
  parameter) that is accurate to a tolerance TOL

  NOTES:
*/
#pragma once
#include "UniformLookupTable.hpp"

#include <memory>

class UniformLookupTableGenerator
{
private:

  EvaluationFunctor<double,double> *mp_func;
  double m_min;
  double m_max;

  /* Nested functor for error evaluation */
  struct LookupTableErrorFunctor;

  /* Nested functor for optimal grid spacing determination */
  struct OptimalStepSizeFunctor;

public:

  UniformLookupTableGenerator(EvaluationFunctor<double,double> *func, double minArg, double maxArg);
  ~UniformLookupTableGenerator();
  std::unique_ptr<UniformLookupTable> generate_by_step(std::string tableKey, double stepSize);
  std::unique_ptr<UniformLookupTable> generate_by_tol(std::string tableKey, double desiredTolerance);
  std::unique_ptr<UniformLookupTable> generate_by_impl_size(std::string tableKey, unsigned long desiredSize);
  double error_at_step_size(std::string tableKey, double stepSize);
  void plot_implementation_at_step_size(std::string tableKey, double stepSize);

};

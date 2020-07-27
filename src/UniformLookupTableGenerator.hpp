/*
  Generate an UniformLookupTable of a table type T (template
  parameter) that is accurate to a tolerance TOL

  NOTES:
*/
#pragma once
#include "UniformLookupTable.hpp"
#include "TransferFunction.hpp"
#include <memory>

class UniformLookupTableGenerator
{
private:
  FunctionContainer *mp_func_container;

  std::shared_ptr<TransferFunction> mp_transfer_function;
  double m_min;
  double m_max;

  /* Nested functor for error evaluation */
  struct LookupTableErrorFunctor;

  /* Nested functor for optimal grid spacing determination */
  struct OptimalStepSizeFunctor;

public:
  UniformLookupTableGenerator(FunctionContainer *func_container,
      double minArg, double maxArg, std::shared_ptr<TransferFunction> transferFunction = nullptr);

  ~UniformLookupTableGenerator();
  std::unique_ptr<UniformLookupTable> generate_by_step(std::string tableKey, double stepSize);
  std::unique_ptr<UniformLookupTable> generate_by_tol(std::string tableKey, double desiredTolerance);
  std::unique_ptr<UniformLookupTable> generate_by_impl_size(std::string tableKey, unsigned long desiredSize);
  double error_at_step_size(std::string tableKey, double stepSize);
  void plot_implementation_at_step_size(std::string tableKey, double stepSize);

};

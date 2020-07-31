/*
  Generate an UniformLookupTable of a table type T (template
  parameter) that is accurate to a tolerance TOL

  NOTES:
*/
#pragma once
#include "UniformLookupTable.hpp"
#include "TransferFunction.hpp"
#include <memory>

template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE>
class UniformLookupTableGenerator
{
private:
  FunctionContainer<IN_TYPE,OUT_TYPE> *mp_func_container;

  std::shared_ptr<TransferFunction<IN_TYPE>> mp_transfer_function;
  IN_TYPE m_min;
  IN_TYPE m_max;

  /* Nested functor for error evaluation */
  struct LookupTableErrorFunctor;

  /* Nested functor for optimal grid spacing determination */
  struct OptimalStepSizeFunctor;

public:
  UniformLookupTableGenerator(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container,
      IN_TYPE minArg, IN_TYPE maxArg,
      std::shared_ptr<TransferFunction<IN_TYPE>> transferFunction = nullptr);

  ~UniformLookupTableGenerator(){}
  std::unique_ptr<UniformLookupTable<IN_TYPE,OUT_TYPE>> generate_by_step(std::string tableKey, IN_TYPE stepSize);
  std::unique_ptr<UniformLookupTable<IN_TYPE,OUT_TYPE>> generate_by_tol(std::string tableKey, double desiredTolerance);
  std::unique_ptr<UniformLookupTable<IN_TYPE,OUT_TYPE>> generate_by_impl_size(std::string tableKey, unsigned long desiredSize);
  double error_at_step_size(std::string tableKey, IN_TYPE stepSize);
  void plot_implementation_at_step_size(std::string tableKey, IN_TYPE stepSize);

};

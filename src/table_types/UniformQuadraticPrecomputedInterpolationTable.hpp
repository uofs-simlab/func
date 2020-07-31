/*
  Quadratic Interpolation LUT with uniform sampling (precomputed coefficients)

  Usage example:
    UniformQuadraticPrecomputedInterpolationTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "UniformLookupTable.hpp"

template <typename IN_TYPE, typename OUT_TYPE>
class UniformQuadraticPrecomputedInterpolationTable final : public UniformLookupTable<IN_TYPE,OUT_TYPE>
{
  REGISTER_ULUT(UniformQuadraticPrecomputedInterpolationTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,3,32>[]> m_table;
public:
  UniformQuadraticPrecomputedInterpolationTable(FunctionContainer *func_container, UniformLookupTableParameters par);
  OUT_TYPE operator()(IN_TYPE x) override;
};

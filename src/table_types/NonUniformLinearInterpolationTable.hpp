/*
  Linear Interpolation LUT with nonuniform sampling

  Usage example:
    NonUniformLinearInterpolationTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "NonUniformLookupTable.hpp"

template <typename IN_TYPE, typename OUT_TYPE>
class NonUniformLinearInterpolationTable final : public NonUniformLookupTable<IN_TYPE,OUT_TYPE>
{
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE);

  REGISTER_LUT(NonUniformLinearInterpolationTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,1>[]> m_table;
public:
  NonUniformLinearInterpolationTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container,
      UniformLookupTableParameters<IN_TYPE> par);
  OUT_TYPE operator()(IN_TYPE x) override;
};

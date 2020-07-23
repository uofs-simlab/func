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

class NonUniformLinearInterpolationTable final : public NonUniformLookupTable
{
  REGISTER_ULUT(NonUniformLinearInterpolationTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<1,8>[]> m_table;
public:
  NonUniformLinearInterpolationTable(FunctionContainer *func_container, UniformLookupTableParameters par);
  double operator()(double x) override;
};

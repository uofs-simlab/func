/*
  Linear Interpolation LUT with uniform sampling

  Usage example:
    UniformLinearInterpolationTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "UniformLookupTable.hpp"

class UniformLinearInterpolationTable final : public UniformLookupTable
{
  REGISTER_ULUT(UniformLinearInterpolationTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<1,8>[]> m_table;
public:
  UniformLinearInterpolationTable(EvaluationFunctor<double,double> *func, UniformLookupTableParameters par);
  double operator()(double x) override;
};

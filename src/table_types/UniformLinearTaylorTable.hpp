/*
  Linear Taylor LUT with uniform sampling

  Usage example:
    UniformLinearTaylorTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "UniformAutoDiffTable.hpp"

class UniformLinearTaylorTable final : public UniformAutoDiffTable<1>
{
  REGISTER_ULUT_DIFF(1,UniformLinearTaylorTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<2,16>[]> m_table;
public:
  UniformLinearTaylorTable(EvaluationFunctor<autodiff_fvar<double,1>,autodiff_fvar<double,1>> *func, UniformLookupTableParameters par);
  double operator()(double x) override;
};

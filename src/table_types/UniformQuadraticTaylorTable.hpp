/*
  Quadratic Taylor LUT with uniform sampling

  Usage example:
    UniformQuadraticTaylorTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "UniformAutoDiffTable.hpp"

class UniformQuadraticTaylorTable final : public UniformAutoDiffTable<2>
{
  REGISTER_ULUT_DIFF(2,UniformQuadraticTaylorTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<3,32>[]> m_table;
public:
  UniformQuadraticTaylorTable(EvaluationFunctor<autodiff_fvar<double,2>,autodiff_fvar<double,2>> *func, UniformLookupTableParameters par);
  double operator()(double x) override;
};

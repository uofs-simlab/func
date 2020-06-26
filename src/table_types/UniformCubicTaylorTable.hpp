/*
  Cubic Taylor LUT with uniform sampling

  Usage example:
    UniformCubicTaylorTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "UniformAutoDiffTable.hpp"
#include <boost/math/differentiation/autodiff.hpp>

class UniformCubicTaylorTable final : public UniformAutoDiffTable<3>
{
  REGISTER_ULUT_DIFF(3, UniformCubicTaylorTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<4,32>[]> m_table;
public:
  UniformCubicTaylorTable(EvaluationFunctor<autodiff_fvar<double,3>,autodiff_fvar<double,3>> *func, UniformLookupTableParameters par);
  double operator()(double x) override;
};

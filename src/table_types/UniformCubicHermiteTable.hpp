/*
  Cubic Interpolation LUT with uniform sampling (precomputed coefficients)

  Usage example:
    UniformCubicHermiteTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#include "UniformAutoDiffTable.hpp"

class UniformCubicHermiteTable final : public UniformAutoDiffTable<1>
{
  REGISTER_ULUT_DIFF(1,UniformCubicHermiteTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<4,32>[]> m_table;
public:
  UniformCubicHermiteTable(EvaluationFunctor<autodiff_fvar<double,1>,autodiff_fvar<double,1>> *func, UniformLookupTableParameters par);
  double operator()(double x) override;
};

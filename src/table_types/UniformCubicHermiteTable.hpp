/*
  Linear Interpolation LUT with uniform sampling (precomputed coefficients)

  Usage example:
    UniformCubicHermiteTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#include "UniformLookupTable.hpp"

class UniformCubicHermiteTable final : public UniformLookupTable
{
  REGISTER_ULUT(UniformCubicHermiteTable);

public:
  UniformCubicHermiteTable(EvaluationFunctor<double,double> *func, UniformLookupTableParameters par);
  double operator()(double x) override;
};

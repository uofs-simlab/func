/*
  LUT using [M/N] pade approximants with uniform sampling. Polynomial coefficients are calculated using
  Armadillo. 

  Usage example using [4/3] approximants:
    UniformPadeTable<4,3> look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
  - Available template values are M such that 0 < N <= M and M+N<=7 and the EvaluationFunctor's M+Nth
  derivative must be defined.
*/
#include "UniformLookupTable.hpp"

template <unsigned int M, unsigned int N> 
class UniformPadeTable final : public UniformLookupTable
{
  REGISTER_ULUT(UniformPadeTable);

public:
  UniformPadeTable(EvaluationFunctor<double,double> *func, UniformLookupTableParameters par);
  double operator()(double x) override;
};

/*
  4th to 7th degree polynomial interpolation LUT with uniform sampling (precomputed coefficients using an
  Armadillo matrix)

  Usage example for a 4th degree interpolant:
    UniformArmadilloPrecomputedInterpolationTable<4> look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
  - the template implementation is only for N=4,5,6,7 (ie, available polynomial interpolation is 
  of degrees 4 up to degree 7)
*/
#include "UniformLookupTable.hpp"

template <unsigned int N>
class UniformArmadilloPrecomputedInterpolationTable final : public UniformLookupTable
{
  REGISTER_ULUT(UniformArmadilloPrecomputedInterpolationTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<N+1,64>[]> m_table;
public:
  UniformArmadilloPrecomputedInterpolationTable(EvaluationFunctor<double,double> *func, UniformLookupTableParameters par);
  double operator()(double x) override;
};

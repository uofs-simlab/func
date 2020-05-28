/*
  Linear Interpolation LUT with uniform sampling (precomputed coefficients)

  Usage example:
    UniformFailureProofLinearPITable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
  - identical to the standard cubic table but samples from outside the
  table's range are done with the original function
  - if the NDEBUG flag is not specified then arguments outside the table's range
  are recorded and printed by the table's destructor
*/
#include "UniformLookupTable.hpp"

#ifndef NDEBUG
#include <vector>
#endif

class UniformFailureProofLinearPITable final : public UniformLookupTable
{
  REGISTER_ULUT(UniformFailureProofLinearPITable);

  #ifndef NDEBUG
  std::vector<double> m_args;
  #endif
  __attribute__((aligned)) std::unique_ptr<polynomial<2,16>[]> m_table;
public:
  UniformFailureProofLinearPITable(EvaluationFunctor<double,double> *func, UniformLookupTableParameters par);
  double operator()(double x) override;
  ~UniformFailureProofLinearPITable();
};

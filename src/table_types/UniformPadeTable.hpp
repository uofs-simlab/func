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
  - Available template values are all M,N such that 0 < N <= M and M+N<=7
*/
#include "UniformLookupTable.hpp"

template <typename IN_TYPE, typename OUT_TYPE, unsigned int M, unsigned int N> 
class UniformPadeTable final : public UniformLookupTable<IN_TYPE,OUT_TYPE>
{
  REGISTER_ULUT(UniformPadeTable);

  // mind the lesser than sign                                         -> |
  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE, M+N+1, M+N<4 ? 32:64>[]> m_table;
  std::function<fvar<OUT_TYPE,M+N>(fvar<OUT_TYPE,M+N>)> mp_boost_func;
public:
  UniformPadeTable(FunctionContainer *func_container, UniformLookupTableParameters par);
  OUT_TYPE operator()(IN_TYPE x) override;
  std::function<fvar<OUT_TYPE,M+N>(fvar<OUT_TYPE,M+N>)> boost_function(){ return mp_boost_func; }
};

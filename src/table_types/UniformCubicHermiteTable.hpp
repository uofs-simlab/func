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
#include "UniformLookupTable.hpp"

template <typename IN_TYPE, typename OUT_TYPE>
class UniformCubicHermiteTable final : public UniformLookupTable<IN_TYPE,OUT_TYPE>
{
  REGISTER_ULUT(UniformCubicHermiteTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,4,32>[]> m_table;
  std::function<fvar<OUT_TYPE,1>(fvar<IN_TYPE,1>)> mp_boost_func;
public:
  UniformCubicHermiteTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, UniformLookupTableParameters<IN_TYPE> par);
  OUT_TYPE operator()(IN_TYPE x) override;
  std::function<fvar<OUT_TYPE,1>(fvar<OUT_TYPE,1>)> boost_function(){ return mp_boost_func; }
};

/*
  Linear Interpolation LUT with uniform sampling

  Usage example:
    UniformLinearInterpolationTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "UniformLookupTable.hpp"

template <typename IN_TYPE, typename OUT_TYPE>
class UniformLinearInterpolationTable final : public UniformLookupTable<IN_TYPE,OUT_TYPE>
{
  REGISTER_ULUT(UniformLinearInterpolationTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,1,8>[]> m_table;
public:
  UniformLinearInterpolationTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, UniformLookupTableParameters<IN_TYPE> par);
  OUT_TYPE operator()(IN_TYPE x) override;
};

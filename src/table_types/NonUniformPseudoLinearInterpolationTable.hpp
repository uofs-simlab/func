/*
  Linear Interpolation LUT with nonuniform sampling

  Usage example:
    NonUniformPseudoLinearInterpolationTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "NonUniformLookupTable.hpp"

template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE, class TRANSFER_FUNC_TYPE = TransferFunctionSinh<IN_TYPE>>
class NonUniformPseudoLinearInterpolationTable final : public NonUniformLookupTable<IN_TYPE,OUT_TYPE,TRANSFER_FUNC_TYPE>
{
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE);
  INHERIT_NONUNIFORM_LUT(IN_TYPE,OUT_TYPE);

  REGISTER_LUT(NonUniformPseudoLinearInterpolationTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,1>[]> m_table;
public:
  NonUniformPseudoLinearInterpolationTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container,
      UniformLookupTableParameters<IN_TYPE> par) :
    NonUniformLookupTable<IN_TYPE,OUT_TYPE,TRANSFER_FUNC_TYPE>(func_container, par)
  {
    /* Base class variables */
    m_name = STR(NonUniformPseudoLinearInterpolationTable);
    m_order = 2;
    m_numTableEntries = m_numIntervals;
    m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;
    
    /* Allocate and set table */
    m_table.reset(new polynomial<OUT_TYPE,1>[m_numTableEntries]);
    for (int ii=0; ii<m_numIntervals; ++ii) {
      // transform the previously used uniform grid to a nonuniform grid
      const IN_TYPE x = m_transferFunction.g(m_minArg + ii*m_stepSize);
      m_grid[ii]  = x;
      m_table[ii].coefs[0] = mp_func(x);
    }
  }

  OUT_TYPE operator()(IN_TYPE x) override
  {
    // set x0 = floor((g_inv(x)-m_minArg)/m_stepSize)
    // and dx = fractional((g_inv(x)-m_minArg)/m_stepSize)
    // where each of the above member vars are encoded into g_inv
    // The source of the pseudolinearity is from the way we compute dx
    OUT_TYPE dx = m_transferFunction.g_inv(x);
    unsigned x0 = (unsigned) dx;
    dx -= x0;

    // value of table entries around x position
    OUT_TYPE y1  = m_table[x0].coefs[0];
    OUT_TYPE y2  = m_table[x0+1].coefs[0];
    // linear interpolation
    return y1+dx*(y2-y1);
  }
};

REGISTER_NONUNIFORM_IMPL(NonUniformPseudoLinearInterpolationTable,double,double,TransferFunctionSinh<double>);

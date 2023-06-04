/*
  Linear Interpolation LUT. Coefficients are computed at lookup time.
  Approx 50% less memory usage compared to EqSpaceInterpTable<1> but the hash involves an additional subtraction.

  Usage example:
    LinearRawInterpTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
  - Does not have a nonuniform variant and it's not obvious how to make one unless we make the operator() far slower (lookup data from m_grid?)
*/
#pragma once
#include "MetaTable.hpp"

namespace func {

template <typename TIN, typename TOUT=TIN, GridTypes GT=GridTypes::UNIFORM>
class LinearRawInterpTable final : public MetaTable<1,TIN,TOUT,GT>
{
  INHERIT_META(1,TIN,TOUT,GT);
public:
  // build the LUT from scratch or look in jsonStats for an existing LUT
  LinearRawInterpTable(const FunctionContainer<TIN,TOUT>& func_container, const LookupTableParameters<TIN>& par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<1,TIN,TOUT,GT>(jsonStats.empty() ? // use the default move constructor for MetaTable (probably not elided...)
      std::move(MetaTable<1,TIN,TOUT,GT>(func_container, par)) :
      std::move(MetaTable<1,TIN,TOUT,GT>(jsonStats)))
  {
    if(!jsonStats.empty())
      return; // all our work is already done

    /* Base class variables */
    m_name = grid_type_to_string<GT>() + "LinearRawInterpTable";
    m_order = 2;
    m_numTableEntries = m_numIntervals+2;
    m_dataSize = sizeof(m_table[0]) * m_numTableEntries;

    auto fun = func_container.standard_fun;

    /* Allocate and set table */
    m_grid.reset(new TIN[m_numTableEntries]);
    m_table.reset(new polynomial<TOUT,1>[m_numTableEntries]);

    FUNC_BUILDPAR
    for (unsigned int ii=0; ii<m_numTableEntries-1; ++ii) {
      TIN x = m_minArg + ii*m_stepSize;
      m_grid[ii] = x;
      m_table[ii].coefs[0] = fun(x);
    }
    // special case to make lut(tableMaxArg) work
    m_table[m_numTableEntries-1].coefs[0] = m_table[m_numTableEntries-2].coefs[0];
  }

  /* this operator() is slightly different from MetaTable's provided Horner's method
   * TODO is there a good way to make this work with nonuniform grids? */
  TOUT operator()(TIN x) const override
  {
    unsigned int x0; TOUT dx;
    std::tie(x0,dx) = MetaTable<1,TIN,TOUT,GT>::template hash<GT>(x);

    // linear interpolation
    TOUT y1 = m_table[x0].coefs[0];
    TOUT y2 = m_table[x0+1].coefs[0];
    return y1+dx*(y2-y1);
  }
};

template <typename TIN, typename TOUT=TIN>
using UniformLinearRawInterpTable = LinearRawInterpTable<TIN,TOUT,GridTypes::UNIFORM>;
} // namespace func

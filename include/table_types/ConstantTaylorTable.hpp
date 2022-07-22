/*
  Constant Taylor LUT with uniform sampling

  Usage example:
    ConstantTaylorTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "MetaTable.hpp"

namespace func {

template <typename TIN, typename TOUT=TIN, GridTypes GT=UNIFORM>
class ConstantTaylorTable final : public MetaTable<TIN,TOUT,1,GT>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,1,GT);

  static const std::string classname;
public:
  // build the LUT from scratch or look in filename for an existing LUT
  ConstantTaylorTable(FunctionContainer<TIN,TOUT> *func_container, LookupTableParameters<TIN> par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<TIN,TOUT,1,GT>(jsonStats.empty() ? // use the default move constructor for MetaTable (probably not elided...)
      std::move(MetaTable<TIN,TOUT,1,GT>(func_container, par)) :
      std::move(MetaTable<TIN,TOUT,1,GT>(jsonStats, classname, func_container)))
  {
    if(!jsonStats.empty())
      return; // all our work is already done

    /* Base class default variables */
    m_name = classname;
    m_order = 1;
    m_numTableEntries = m_numIntervals;
    m_dataSize = (unsigned) sizeof(m_table[0]) * (m_numTableEntries);

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,1>[m_numTableEntries]);
    for (unsigned int ii=0; ii<m_numTableEntries; ++ii) {
      TIN x;
      // (possibly) transform the uniform grid into a nonuniform grid
      if (GT == UNIFORM)
        x = m_minArg + ii*m_stepSize;
      else
        x = m_transferFunction.g(m_minArg + ii*m_stepSize);

      m_grid[ii] = x;
      m_table[ii].coefs[0] = m_func(x);
    }
  }

  /* Constant interpolation from table point immediately below x */
  TOUT operator()(TIN x) override
  {
    return m_table[(unsigned)((x-m_minArg)/m_stepSize+0.5)].coefs[0];
  }
};

template <typename TIN, typename TOUT, GridTypes GT>
const std::string ConstantTaylorTable<TIN,TOUT,GT>::classname = grid_type_to_string<GT>() + "ConstantTaylorTable";

// define friendlier names
template <typename TIN, typename TOUT=TIN>
using UniformConstantTaylorTable = ConstantTaylorTable<TIN,TOUT,UNIFORM>;
} // namespace func

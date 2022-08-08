/*
  Linear Interpolation LUT with uniform sampling (precomputed coefficients)

  Usage example:
    LinearPrecomputedInterpolationTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "MetaTable.hpp"

namespace func {

template <typename TIN, typename TOUT=TIN, GridTypes GT=GridTypes::UNIFORM>
class LinearPrecomputedInterpolationTable final : public MetaTable<TIN,TOUT,2,GT>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,2,GT);

  static const std::string classname;
public:
  // build the LUT from scratch or look in filename for an existing LUT
  LinearPrecomputedInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, LookupTableParameters<TIN> par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<TIN,TOUT,2,GT>(jsonStats.empty() ? // use the default move constructor for MetaTable (probably not elided...)
      std::move(MetaTable<TIN,TOUT,2,GT>(func_container, par)) :
      std::move(MetaTable<TIN,TOUT,2,GT>(jsonStats, classname, func_container)))
  {
    if(!jsonStats.empty())
      return; // all our work is already done

    /* Base class default variables */
    m_name = classname;
    m_order = 2;
    m_numTableEntries = m_numIntervals;
    m_dataSize =static_cast<unsigned>(sizeof(m_table[0]) * (m_numTableEntries));

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,2>[m_numTableEntries]);
    for (unsigned int ii=0; ii < m_numTableEntries; ++ii) {
      TIN x;
      TIN h = m_stepSize;
      // (possibly) transform the uniform grid into a nonuniform grid
      if (GT == GridTypes::UNIFORM)
        x = m_minArg + ii*m_stepSize;
      else{
        x = m_transferFunction.g(m_minArg + ii*m_stepSize);
        h = m_transferFunction.g(m_minArg + (ii+1)*m_stepSize) - x;
      }

      m_grid[ii] = x;
      m_table[ii].coefs[0] = m_func(x);
      m_table[ii].coefs[1] = m_func(x+h) - m_table[ii].coefs[0];
    }
  }

  /* build this table from a file. Everything other than m_table is built by MetaTable */
  LinearPrecomputedInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    MetaTable<TIN,TOUT,2,GT>(func_container, filename,
        grid_type_to_string<GT>() + "LinearPrecomputedInterpolationTable") {}
  // operator() comes from MetaTable
};

template <typename TIN, typename TOUT, GridTypes GT>
const std::string LinearPrecomputedInterpolationTable<TIN,TOUT,GT>::classname = grid_type_to_string<GT>() + "LinearPrecomputedInterpolationTable";

// define friendlier names
template <typename TIN, typename TOUT=TIN>
using UniformLinearPrecomputedInterpolationTable = LinearPrecomputedInterpolationTable<TIN,TOUT,GridTypes::UNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformLinearPrecomputedInterpolationTable = LinearPrecomputedInterpolationTable<TIN,TOUT,GridTypes::NONUNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformPseudoLinearPrecomputedInterpolationTable = LinearPrecomputedInterpolationTable<TIN,TOUT,GridTypes::NONUNIFORM_PSEUDO>;
} // namespace func

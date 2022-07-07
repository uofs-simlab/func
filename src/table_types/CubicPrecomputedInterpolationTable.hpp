/*
  Cubic Interpolation LUT with precomputed coefficients

  Usage example:
    CubicPrecomputedInterpolationTable look(&function,{0,10,0.0001});
    double val = look(0.87354);

  Notes:
  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - points sampled in each subintervale will be uniform b/c that can
    reuse uniform code and seems to be more well behaved in the end.
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "MetaTable.hpp"

template <typename TIN, typename TOUT=TIN, GridTypes GT=UNIFORM>
class CubicPrecomputedInterpolationTable final : public MetaTable<TIN,TOUT,4,GT>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,4,GT);

  static const std::string classname;
public:
  // build the LUT from scratch or look in filename for an existing LUT
  CubicPrecomputedInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, LookupTableParameters<TIN> par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<TIN,TOUT,4,GT>(jsonStats.empty() ? // use the default move constructor for MetaTable (probably not elided...)
      std::move(MetaTable<TIN,TOUT,4,GT>(func_container, par)) :
      std::move(MetaTable<TIN,TOUT,4,GT>(jsonStats, classname, func_container)))
  {
    if(!jsonStats.empty())
      return; // all our work is already done

    /* Base class default variables */
    m_name = classname;
    m_order = 4;
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = (unsigned) sizeof(m_table[0]) * (m_numTableEntries);

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,4>[m_numTableEntries]);
    for (unsigned int ii=0;ii<m_numIntervals;++ii) {
      TIN x;
      TIN h = m_stepSize;
      // (possibly) transform the uniform grid into a nonuniform grid
      if (GT == UNIFORM)
        x = m_minArg + ii*m_stepSize;
      else{
        x = m_transferFunction.g(m_minArg + ii*m_stepSize);
        h = m_transferFunction.g(m_minArg + (ii+1)*m_stepSize) - x;
      }

      // grid points
      m_grid[ii] = x;
      // polynomial coefficients:
      const TOUT y0 = m_func(x);
      const TOUT y1 = m_func(x+h/3);
      const TOUT y2 = m_func(x+2*h/3);
      const TOUT y3 = m_func(x+h);
      m_table[ii].coefs[0] = y0;
      m_table[ii].coefs[1] = -11*y0/2+9*y1-9*y2/2+y3;
      m_table[ii].coefs[2] = 9*y0-45*y1/2+18*y2-9*y3/2;
      m_table[ii].coefs[3] = -9*y0/2+27*y1/2-27*y2/2+9*y3/2;
    }
  }

  /* build this table from a file. Everything other than m_table is built by MetaTable */
  CubicPrecomputedInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    MetaTable<TIN,TOUT,4,GT>(func_container, filename,
        grid_type_to_string<GT>() + "CubicPrecomputedInterpolationTable") {}
  // operator() comes straight from the MetaTable
};

template <typename TIN, typename TOUT, GridTypes GT>
const std::string CubicPrecomputedInterpolationTable<TIN,TOUT,GT>::classname = grid_type_to_string<GT>() + "CubicPrecomputedInterpolationTable";

// define friendlier names
template <typename TIN, typename TOUT=TIN>
using UniformCubicPrecomputedInterpolationTable = CubicPrecomputedInterpolationTable<TIN,TOUT,UNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformCubicPrecomputedInterpolationTable = CubicPrecomputedInterpolationTable<TIN,TOUT,NONUNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformPseudoCubicPrecomputedInterpolationTable = CubicPrecomputedInterpolationTable<TIN,TOUT,NONUNIFORM_PSEUDO>;

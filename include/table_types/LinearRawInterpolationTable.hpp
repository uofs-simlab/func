/*
  Linear Interpolation LUT. Coefficients are computed at lookup time.
  Approx 50% less memory usage compared to LinearPrecomputedInterpolationTable
  but the hash is slower.

  Usage example:
    LinearRawInterpolationTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "MetaTable.hpp"

namespace func {

template <typename TIN, typename TOUT=TIN, GridTypes GT=GridTypes::UNIFORM>
class LinearRawInterpolationTable final : public MetaTable<TIN,TOUT,1,GT>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,1,GT);

  static const std::string classname;
public:
  // build the LUT from scratch or look in filename for an existing LUT
  //#pragma omp declare simd
  LinearRawInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, LookupTableParameters<TIN> par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<TIN,TOUT,1,GT>(jsonStats.empty() ? // use the default move constructor for MetaTable (probably not elided...)
      std::move(MetaTable<TIN,TOUT,1,GT>(func_container, par)) :
      std::move(MetaTable<TIN,TOUT,1,GT>(jsonStats, classname, func_container)))
  {
    if(!jsonStats.empty())
      return; // all our work is already done

    /* Base class variables */
    m_name  = classname;
    m_order = 2;
    m_numTableEntries = m_numIntervals+2; // +2 because hash uses two array entries
    m_dataSize = static_cast<unsigned>(sizeof(m_table[0]) * (m_numTableEntries));

    /* Allocate and set table */
    m_grid.reset(new TIN[m_numTableEntries]);
    m_table.reset(new polynomial<TOUT,1>[m_numTableEntries]);
    for (unsigned int ii=0; ii<m_numTableEntries-1; ++ii) {
      TIN x;
      // (possibly) transform the uniform grid into a nonuniform grid
      if (GT == GridTypes::UNIFORM)
        x = m_minArg + ii*m_stepSize;
      else
        x = m_transferFunction.g(m_minArg + ii*m_stepSize);

      m_grid[ii] = x;
      m_table[ii].coefs[0] = m_func(x);
    }
    // special case to make lut(tableMaxArg) work
    m_table[m_numTableEntries-1].coefs[0] = m_table[m_numTableEntries-2].coefs[0];
  }

  /* build this table from a file. Everything other than m_table is built by MetaTable */
  LinearRawInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    MetaTable<TIN,TOUT,1,GT>(func_container, filename,
        grid_type_to_string<GT>() + "LinearRawInterpolationTable") {}

  // operator() is slightly different from MetaTable's provided Horner's method
  TOUT operator()(TIN x) override
  {
    //enum class GridTypes {UNIFORM, NONUNIFORM, NONUNIFORM_PSEUDO};
    // hash is copied from MetaTable
    TOUT dx;
    unsigned int x0;
    switch(GT){
    case GridTypes::UNIFORM:
      {
      // nondimensionalized x position, scaled by step size
      dx = static_cast<TOUT>(m_stepSize_inv*(x-m_minArg));
      // index of previous table entry
      x0 = static_cast<unsigned>(dx);
      // value of table entries around x position
      dx -= x0;
      break;
      }
    case GridTypes::NONUNIFORM:
      {
      x0 = static_cast<unsigned>(m_transferFunction.g_inv(x));
      TIN h   = m_grid[x0+1] - m_grid[x0];
      dx = (x - m_grid[x0])/h;
      break;
      }
    case GridTypes::NONUNIFORM_PSEUDO:
      {
      dx = m_transferFunction.g_inv(x);
      x0 = static_cast<unsigned>(dx);
      dx -= x0;
      break;
      }
    }

    // linear interpolation
    TOUT y1  = m_table[x0].coefs[0];
    TOUT y2  = m_table[x0+1].coefs[0];
    return y1+dx*(y2-y1);
  }
};

template <typename TIN, typename TOUT, GridTypes GT>
const std::string LinearRawInterpolationTable<TIN,TOUT,GT>::classname = grid_type_to_string<GT>() + "LinearRawInterpolationTable";

template <typename TIN, typename TOUT=TIN>
using UniformLinearRawInterpolationTable = LinearRawInterpolationTable<TIN,TOUT,GridTypes::UNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformLinearRawInterpolationTable = LinearRawInterpolationTable<TIN,TOUT,GridTypes::NONUNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformPseudoLinearRawInterpolationTable = LinearRawInterpolationTable<TIN,TOUT,GridTypes::NONUNIFORM_PSEUDO>;
} // namespace func

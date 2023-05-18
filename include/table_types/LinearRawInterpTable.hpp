/*
  Linear Interpolation LUT. Coefficients are computed at lookup time.
  Approx 50% less memory usage compared to LinearInterpTable
  but the hash is slower.

  Usage example:
    LinearRawInterpTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "MetaTable.hpp"

namespace func {

template <typename TIN, typename TOUT=TIN, GridTypes GT=GridTypes::UNIFORM>
class LinearRawInterpTable final : public MetaTable<TIN,TOUT,1,GT>
{
  //INHERIT_META(TIN,TOUT,1,GT);

  //static const std::string classname;
  static constexpr const char * classname = "LinearRawInterpTable";
public:
  // build the LUT from scratch or look in filename for an existing LUT
  //#pragma omp declare simd
  LinearRawInterpTable(const FunctionContainer<TIN,TOUT>& func_container, const LookupTableParameters<TIN>& par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<TIN,TOUT,1,GT>(jsonStats.empty() ? // use the default move constructor for MetaTable (probably not elided...)
      std::move(MetaTable<TIN,TOUT,1,GT>(func_container, par)) :
      std::move(MetaTable<TIN,TOUT,1,GT>(jsonStats))),
    m_name(classname), m_order(2), m_numTableEntries(m_numIntervals+2)
  {
    if(!jsonStats.empty())
      return; // all our work is already done

    /* Base class variables */
    m_dataSize = static_cast<unsigned>(sizeof(m_table[0]) * m_numTableEntries);

    /* Allocate and set table */
    m_grid.reset(new TIN[m_numTableEntries]);
    m_table.reset(new polynomial<TOUT,1>[m_numTableEntries]);

    FUNC_BUILDPAR
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

  // this operator() is slightly different from MetaTable's provided Horner's method
  TOUT operator()(TIN x) override
  {
    unsigned int x0; TOUT dx;
    std::tie(x0,dx) = MetaTable<TIN,TOUT,1,GT>::hash(x);

    // linear interpolation
    TOUT y1 = m_table[x0].coefs[0];
    TOUT y2 = m_table[x0+1].coefs[0];
    return y1+dx*(y2-y1);
  }
};

//template <typename TIN, typename TOUT, GridTypes GT>
//const std::string LinearRawInterpTable<TIN,TOUT,GT>::classname = grid_type_to_string<GT>() + "LinearRawInterpTable";

template <typename TIN, typename TOUT=TIN>
using UniformLinearRawInterpTable = LinearRawInterpTable<TIN,TOUT,GridTypes::UNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformLinearRawInterpTable = LinearRawInterpTable<TIN,TOUT,GridTypes::NONUNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformPseudoLinearRawInterpTable = LinearRawInterpTable<TIN,TOUT,GridTypes::NONUNIFORM_PSEUDO>;
} // namespace func

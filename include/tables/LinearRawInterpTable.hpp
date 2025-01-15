#pragma once
#include "MetaTable.hpp"

namespace func {

/**
  \brief Linear Interpolation LUT where coefficients are computed when calling operator().
  Uses approx 50% less memory than an equivalent UniformExactInterpTable<1>
  but the hash involves an additional subtraction.
  \ingroup MetaTable

  \code{.cpp}
  // LinearRawInterpTable does not benefit from templated functions
  double foo(double x){ return x; }
 
  int main(){
    double min = 0.0, max = 10.0, step = 0.0001;
    UniformLinearRawInterpTable<double> L({foo}, {min, max, step});
    auto val = L(0.87354);
  }
  \endcode


  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
  - Does not have a nonuniform variant and it's not obvious how to make
    this LookupTable implementation nonuniform unless we make the operator()
    far slower (basically defeating the purpose of this LUT type e.g. lookup
    breakpoints from m_grid?)
*/
template <typename TIN, typename TOUT=TIN, GridTypes GT=GridTypes::UNIFORM>
class LinearRawInterpTable final : public MetaTable<1,TIN,TOUT,GT>
{
  INHERIT_META(1,TIN,TOUT,GT);
public:
  LinearRawInterpTable() = default;
  // This LUT _must_ use a different operator() than the base class MetaTable
  //LinearRawInterpTable(const MetaTable<N+1,TIN,TOUT,GT>& L): MetaTable<N+1,TIN,TOUT,GT>(L) {}

  // Either build the LUT from scratch or read data from json
  LinearRawInterpTable(const FunctionContainer<TIN,TOUT>& func_container, const LookupTableParameters<TIN>& par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<1,TIN,TOUT,GT>(func_container, par, jsonStats)
  {
    if(!jsonStats.empty())
      return; // all our work is already done

    /* Base class variables */
    m_name = grid_type_to_string<GT>() + "LinearRawInterpTable";
    m_order = 2;
    m_numTableEntries = m_numIntervals+2;
    m_dataSize = sizeof(m_table[0]) * m_numTableEntries;

    auto fun = func_container.standard_fun;
    if(fun == nullptr)
      throw std::invalid_argument("Error in func::LinearRawInterpTable: Given an invalid FunctionContainer");

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,1>[m_numTableEntries]);

    FUNC_BUILDPAR
    for (unsigned int ii=0; ii<m_numTableEntries-1; ++ii) {
      TIN x = m_minArg + ii*m_stepSize;
      m_table[ii].coefs[0] = fun(x);
    }

    /* special case to make lut(tableMaxArg) work */
    m_table[m_numTableEntries-1] = taylor_shift(m_table[m_numTableEntries-2], static_cast<TIN>(1), static_cast<TIN>(2), static_cast<TIN>(0), static_cast<TIN>(1));
  }

  /* this operator() is slightly different from MetaTable's provided Horner's method
   * TODO is there a way to make this work with nonuniform grids in a way that works with our model? */
  TOUT operator()(TIN x) const override
  {
    unsigned int x0; TIN dx;
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

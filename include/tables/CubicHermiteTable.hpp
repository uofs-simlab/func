#pragma once
#include "MetaTable.hpp"
#include "config.hpp" // FUNC_USE_BOOST
#include <stdexcept>

namespace func {

/** \brief A LUT using cubic splines on each subinterval
 * \ingroup MetaTable
 *
 * \code{.cpp}
 * // CubicHermiteTable requires the user's mathematical function is templated
 * template <typename T>
 * T foo(T x){ return x; }
 *
 * int main(){
 *   double min = 0.0, max = 10.0, step = 0.0001;
 *   UniformCubicHermiteTable<double>    L({FUNC_SET_F(foo,double)}, {min, max, step}); // uniform partition
 *   NonUniformCubicHermiteTable<double> L({FUNC_SET_F(foo,double)}, {min, max, step}); // nonuniform partition
 *   auto val = L(0.87354);
 * }
 * \endcode */
template <typename TIN, typename TOUT=TIN, GridTypes GT=GridTypes::UNIFORM>
class CubicHermiteTable final : public MetaTable<4,TIN,TOUT,GT>
{
  INHERIT_META(4,TIN,TOUT,GT);
public:
  CubicHermiteTable() = default;
  CubicHermiteTable(const MetaTable<4,TIN,TOUT,GT>& L): MetaTable<4,TIN,TOUT,GT>(L) {}

  // Either build the LUT from scratch or read data from json
  CubicHermiteTable(const FunctionContainer<TIN,TOUT>& func_container, const LookupTableParameters<TIN>& par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<4,TIN,TOUT,GT>(func_container, par, jsonStats)
  {
#ifndef FUNC_USE_BOOST
    if(jsonStats.empty())
      throw std::invalid_argument("Error in func::CubicHermiteTable: CubicHermite LUTs need Armadillo to be generated");
#else
    if(!jsonStats.empty())
      return; // all our work is already done

    using boost::math::differentiation::make_fvar;
    /* Base class default variables */
    m_name = grid_type_to_string<GT>() + "CubicHermiteTable";
    m_order = 4;
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = sizeof(m_table[0]) * m_numTableEntries;

    auto fun = func_container.standard_fun;
    auto boost_fun = func_container.autodiff1_fun;
    if(boost_fun == nullptr)
      throw std::invalid_argument("Error in func::CubicHermiteTable: 1st derivative of given function is not provided");

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,4>[m_numTableEntries]);
    FUNC_BUILDPAR
    for (unsigned int ii=0; ii<m_numTableEntries-1; ++ii) {
      TIN x;
      TIN h = m_stepSize;
      // (possibly) transform the uniform grid into a nonuniform grid
      FUNC_IF_CONSTEXPR(GT == GridTypes::UNIFORM)
        x = m_minArg + ii*m_stepSize;
      else{
        x = m_transferFunction(m_minArg + ii*m_stepSize);
        h = m_transferFunction(m_minArg + (ii+1)*m_stepSize) - x;
      }

      const auto derivs0 = boost_fun(make_fvar<TIN,1>(x));
      const TOUT y0 = derivs0.derivative(0);
      const TOUT m0 = derivs0.derivative(1);
      const auto derivs1 = boost_fun(make_fvar<TIN,1>(x+h));
      const TOUT y1 = derivs1.derivative(0);
      const TOUT m1 = derivs1.derivative(1);

      m_table[ii].coefs[0] = y0;
      m_table[ii].coefs[1] = h*m0;
      m_table[ii].coefs[2] = -3*y0+3*y1-(2*m0+m1)*h;
      m_table[ii].coefs[3] = 2*y0-2*y1+(m0+m1)*h;

      FUNC_IF_CONSTEXPR(GT == GridTypes::NONUNIFORM){
        auto p = m_table[ii];
        for(unsigned int k=0; k<4; k++)
          m_table[ii].coefs[k] = polynomial_diff(p,-x/h,k)/static_cast<TIN>(pow(h,k))/static_cast<TIN>(factorial(k));
      }
    }
    // special case to make lut(tableMaxArg) work
    m_table[m_numTableEntries-1].coefs[0] = func_container.standard_fun(m_tableMaxArg);
    for (unsigned int k=1; k<4; k++)
      m_table[m_numTableEntries-1].coefs[k] = static_cast<TIN>(0)*m_table[m_numTableEntries-1].coefs[0];
#endif
  }

  // operator() is in MetaTable
};

// define friendlier names
template <typename TIN, typename TOUT=TIN>
using UniformCubicHermiteTable = CubicHermiteTable<TIN,TOUT,GridTypes::UNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformCubicHermiteTable = CubicHermiteTable<TIN,TOUT,GridTypes::NONUNIFORM>;
} // namespace func

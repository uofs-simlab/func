/*
  Truncated Taylor LUT of degree N

  Usage example:
    TaylorTable look(&function,0,10,0.0001);
    double val = look(0.87354);
  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "MetaTable.hpp"
#include "FunctionContainer.hpp"
#include "config.hpp" // FUNC_USE_BOOST
#include <stdexcept>

namespace func {

template <unsigned int N, typename TIN, typename TOUT=TIN, GridTypes GT=GridTypes::UNIFORM>
class TaylorTable final : public MetaTable<N+1,TIN,TOUT,GT>
{
  INHERIT_META(N+1,TIN,TOUT,GT);
public:
  // build the LUT from scratch or look in filename for an existing LUT
  TaylorTable(const FunctionContainer<TIN,TOUT>& func_container, const LookupTableParameters<TIN>& par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<N+1,TIN,TOUT,GT>(jsonStats.empty() ? // use the default move constructor for MetaTable (probably not elided...)
      std::move(MetaTable<N+1,TIN,TOUT,GT>(func_container, par)) :
      std::move(MetaTable<N+1,TIN,TOUT,GT>(jsonStats)))
  {
#ifndef FUNC_USE_BOOST
    /* This could theoretically be a compile time error; however, that will only stop us from registering this table (which is not useful!) */
    if(jsonStats.empty())
      throw std::invalid_argument("Error in func::TaylorTable: Boost version 1.71.0 or newer is not available but FunC must use Boost's automatic differentiation to compute Taylor sums");
#else
    if(!jsonStats.empty())
      return; // all our work is already done

    using boost::math::differentiation::make_fvar;

    /* Base class default variables */
    m_name = grid_type_to_string<GT>() + "TaylorTable<" + std::to_string(N) + ">";
    m_order = N+1;
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = static_cast<unsigned>(sizeof(m_table[0]) * m_numTableEntries);

    auto boost_fun = func_container.template get_nth_func<N>();
    if(boost_fun == nullptr)
      throw std::invalid_argument(m_name+" needs the " + std::to_string(N) + "th derivative but this is not defined");

    /* Allocate and set table */
    m_grid.reset(new TIN[m_numTableEntries]);
    m_table.reset(new polynomial<TOUT,N+1>[m_numTableEntries]);
    FUNC_BUILDPAR
    for (unsigned int ii=0; ii<m_numTableEntries-1; ++ii) {
      auto x = m_minArg + ii*m_stepSize;
      auto h = m_stepSize;
      // (possibly) transform the uniform grid into a nonuniform grid
      if (GT != GridTypes::UNIFORM){
        x = m_transferFunction(x);
        h = m_transferFunction(m_minArg + (ii+1)*m_stepSize) - x;
      }
      m_grid[ii] = x;

      /* Taylor expansion of f over the basis: (xh+0.5h)^k for k=0,1,...,N */
      auto const derivs = boost_fun(make_fvar<TIN,N>(x + 0.5*h));
      for(unsigned int k=0; k<N+1; k++)
        m_table[ii].coefs[k] = derivs.derivative(k)*pow(h,k)/boost::math::factorial<double>(k);

      /* Taylor expansion of the above polynomial over the basis: x^k */
      auto p = m_table[ii];
      for(unsigned int k=0; k<N+1; k++)
        m_table[ii].coefs[k] = polynomial_diff(p,-0.5,k)/boost::math::factorial<double>(k);

      if(GT == GridTypes::NONUNIFORM){
        auto p = m_table[ii];
        for(unsigned int k=0; k<N+1; k++)
          m_table[ii].coefs[k] = polynomial_diff(p,-x/h,k)/pow(h,k)/boost::math::factorial<double>(k);
      }
    }
    // special case to make lut(tableMaxArg) work
    m_grid[m_numTableEntries-1] = m_tableMaxArg;
    m_table[m_numTableEntries-1].coefs[0] = func_container.standard_fun(m_tableMaxArg);
    for (unsigned int k=1; k<N+1; k++)
      m_table[m_numTableEntries-1].coefs[k] = 0;
#endif
  }
};

/* define friendlier names */
template <unsigned int N, typename TIN, typename TOUT=TIN>
using UniformTaylorTable = TaylorTable<N,TIN,TOUT,GridTypes::UNIFORM>;
template <unsigned int N, typename TIN, typename TOUT=TIN>
using NonUniformTaylorTable = TaylorTable<N,TIN,TOUT,GridTypes::NONUNIFORM>;
} // namespace func

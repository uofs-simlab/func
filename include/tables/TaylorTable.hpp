#pragma once
#include "MetaTable.hpp"
#include "FunctionContainer.hpp"
#include "config.hpp" // FUNC_USE_BOOST
#include <stdexcept>

namespace func {

/** \brief LUT using degree 1 to 7 truncated Taylor series 
    \ingroup MetaTable

  Usage example:
    TaylorTable look(&function,0,10,0.0001);
    double val = look(0.87354);
  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
template <unsigned int N, typename TIN, typename TOUT=TIN, GridTypes GT=GridTypes::UNIFORM>
class TaylorTable final : public MetaTable<N+1,TIN,TOUT,GT>
{
  INHERIT_META(N+1,TIN,TOUT,GT);
public:
  TaylorTable() = default;
  TaylorTable(const MetaTable<N+1,TIN,TOUT,GT>& L): MetaTable<N+1,TIN,TOUT,GT>(L) {}

  // build the LUT from scratch or look in filename for an existing LUT
  TaylorTable(const FunctionContainer<TIN,TOUT>& func_container, const LookupTableParameters<TIN>& par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<N+1,TIN,TOUT,GT>(func_container, par, jsonStats)
  {
#ifndef FUNC_USE_BOOST
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
    m_dataSize = sizeof(m_table[0]) * m_numTableEntries;

    auto boost_fun = func_container.template get_nth_func<N>();
    if(boost_fun == nullptr)
      throw std::invalid_argument(m_name+" needs the " + std::to_string(N) + "th derivative but this is not defined");

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,N+1>[m_numTableEntries]);
    FUNC_BUILDPAR
    for (unsigned int ii=0; ii<m_numTableEntries-1; ++ii) {
      TIN x = m_minArg + ii*m_stepSize;
      TIN h = m_stepSize;
      // (possibly) transform the uniform grid into a nonuniform grid
      FUNC_IF_CONSTEXPR(GT != GridTypes::UNIFORM){
        x = m_transferFunction(x);
        h = m_transferFunction(m_minArg + (ii+1)*m_stepSize) - x;
      }

      /* Taylor expansion of f over the basis: (x+0.5h)^k for k=0,1,...,N. Resulting polynomial maps [-0.5,0.5]->\R */
      auto const derivs = boost_fun(make_fvar<TIN,N>(x + 0.5*h));
      for(unsigned int k=0; k<N+1; k++)
        m_table[ii].coefs[k] = derivs.derivative(k)/static_cast<TIN>(factorial(k));
        //m_table[ii].coefs[k] = derivs.derivative(k)*static_cast<TIN>(pow(h,k))/static_cast<TIN>(factorial(k));

      /* Taylor expansion of the above polynomial over the basis: x^k */
      FUNC_IF_CONSTEXPR(GT == GridTypes::UNIFORM)
        m_table[ii] = taylor_shift(m_table[ii], static_cast<TIN>(-0.5)*h, static_cast<TIN>(0.5)*h, static_cast<TIN>(0.0), static_cast<TIN>(1.0));
      else
        m_table[ii] = taylor_shift(m_table[ii], static_cast<TIN>(-0.5)*h, static_cast<TIN>(0.5)*h, x, x+h);
    }
    /* special case to make lut(tableMaxArg) work. Move the second last polynomial into the last interval (shifting args) */
    FUNC_IF_CONSTEXPR(GT == GridTypes::UNIFORM)
      m_table[m_numTableEntries-1] = taylor_shift(m_table[m_numTableEntries-2], static_cast<TIN>(1), static_cast<TIN>(2), static_cast<TIN>(0), static_cast<TIN>(1));
    else
      m_table[m_numTableEntries-1] = m_table[m_numTableEntries-2];
#endif
  }
};

/* define friendlier names */
template <unsigned int N, typename TIN, typename TOUT=TIN>
using UniformTaylorTable = TaylorTable<N,TIN,TOUT,GridTypes::UNIFORM>;
template <unsigned int N, typename TIN, typename TOUT=TIN>
using NonUniformTaylorTable = TaylorTable<N,TIN,TOUT,GridTypes::NONUNIFORM>;
} // namespace func

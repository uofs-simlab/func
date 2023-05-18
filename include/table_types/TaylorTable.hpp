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

template <typename TIN, typename TOUT, unsigned int N, GridTypes GT=GridTypes::UNIFORM>
class TaylorTable final : public MetaTable<TIN,TOUT,N+1,GT>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,N+1,GT);

  static const std::string classname;
#ifdef FUNC_USE_BOOST
  std::function<adVar<TOUT,N>(adVar<TOUT,N>)> mp_boost_func;
#endif

public:
  // build the LUT from scratch or look in filename for an existing LUT
  TaylorTable(const FunctionContainer<TIN,TOUT>& func_container, const LookupTableParameters<TIN>& par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<TIN,TOUT,N+1,GT>(jsonStats.empty() ? // use the default move constructor for MetaTable (probably not elided...)
      std::move(MetaTable<TIN,TOUT,N+1,GT>(func_container, par)) :
      std::move(MetaTable<TIN,TOUT,N+1,GT>(jsonStats, classname, func_container)))
  {
#ifndef FUNC_USE_BOOST
    /* This could theoretically be a compile time error; however, that will only stop us from registering this table (not useful!) */
    if(jsonStats.empty())
      throw std::invalid_argument("Error in func::TaylorTable: Boost version 1.71.0 or newer is not available but FunC must use Boost's automatic differentiation to compute Taylor sums");
#else
    if(!jsonStats.empty())
      return; // all our work is already done

    using boost::math::differentiation::make_fvar;

    /* Base class default variables */
    m_name = classname;
    m_order = N+1;
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = static_cast<unsigned>(sizeof(m_table[0]) * m_numTableEntries);

    if(func_container->template get_nth_func<N>() == nullptr)
      throw std::invalid_argument(m_name+" needs the " + std::to_string(N) + "th derivative but this is not defined");
    mp_boost_func = func_container->template get_nth_func<N>();

    /* Allocate and set table */
    m_grid.reset(new TIN[m_numTableEntries]);
    m_table.reset(new polynomial<TOUT,N+1>[m_numTableEntries]);
    FUNC_BUILDPAR
    for (unsigned int ii=0; ii<m_numTableEntries-1; ++ii) {
      auto xgrid = m_minArg + ii*m_stepSize;
      auto h = m_stepSize;
      // (possibly) transform the uniform grid into a nonuniform grid
      if (GT != GridTypes::UNIFORM){
        xgrid = m_transferFunction.g(xgrid);
        h = m_transferFunction.g(m_minArg + (ii+1)*m_stepSize) - xgrid;
      }
      auto xcenter = xgrid + 0.5*h;
      m_grid[ii] = xgrid;

      auto const derivs = (mp_boost_func)(make_fvar<TIN,N>(xcenter));
      // convenient rename:
      std::array<TOUT,N+1> d;
      for(unsigned int k=0; k<N+1; k++)
        d[k] = derivs.derivative(k);

      /* TODO there should be a nice enough closed form solution for the coefs of 
       * a Taylor approx. of arbitrary degree */
      switch(N){
        case 1:
          {
          m_table[ii].coefs[1] = h*d[1];
          m_table[ii].coefs[0] = d[0] - 0.5*h*d[1];
          break;
          }
        case 2:
          {
          m_table[ii].coefs[2] = 0.5*h*h*d[2];
          m_table[ii].coefs[1] = h*(d[1] - 0.5*h*d[2]);
          m_table[ii].coefs[0] = d[0] - 0.5*h*d[1] + 0.125*h*h*d[2];
          break;
          }
        case 3:
          {
          m_table[ii].coefs[3] = h*h*h*d[3]/6;
          m_table[ii].coefs[2] = h*h*(0.5*d[2] - 0.25*h*d[3]);
          m_table[ii].coefs[1] = h*(d[1] - 0.5*h*d[2] + 0.125*h*h*d[3]);
          m_table[ii].coefs[0] = d[0] - 0.5*h*d[1] + 0.125*h*h*d[2] - h*h*h*d[3]/48;
          break;
          }
        default: { throw std::invalid_argument("Broken switch case in func::TaylorTable"); }
      }
    }
    // special case to make lut(tableMaxArg) work
    m_grid[m_numTableEntries-1] = m_tableMaxArg;
    m_table[m_numTableEntries-1].coefs[0] = m_func(m_tableMaxArg);
    for (unsigned int k=1; k<N+1; k++)
      m_table[m_numTableEntries-1].coefs[k] = 0;
#endif
  }

#ifdef FUNC_USE_BOOST
  std::function<adVar<TOUT,N>(adVar<TOUT,N>)> boost_function(){ return mp_boost_func; }
#endif
};

template <std::size_t N, typename TIN, typename TOUT, GridTypes GT>
const std::string TaylorTable<N,TIN,TOUT,GT>::classname = grid_type_to_string<GT>() + "TaylorTable<" + std::to_string(N) + ">";

/* define friendlier names */
template <std::size_t N, typename TIN, typename TOUT=TIN>
using UniformTaylorTable = TaylorTable<N,TIN,TOUT,GridTypes::UNIFORM>;
template <std::size_t N, typename TIN, typename TOUT=TIN>
using NonUniformTaylorTable = TaylorTable<N,TIN,TOUT,GridTypes::NONUNIFORM>;
template <std::size_t N, typename TIN, typename TOUT=TIN>
using NonUniformPseudoTaylorTable = TaylorTable<N,TIN,TOUT,GridTypes::NONUNIFORM_PSEUDO>;
} // namespace func

/*
  Quadratic Taylor LUT with uniform sampling

  Usage example:
    QuadraticTaylorTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "MetaTable.hpp"
#include "FunctionContainer.hpp"
#include "config.hpp" // FUNC_USE_BOOST

namespace func {

template <typename TIN, typename TOUT=TIN, GridTypes GT=GridTypes::UNIFORM>
class QuadraticTaylorTable final : public MetaTable<TIN,TOUT,3,GT>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,3,GT);

  static const std::string classname;
#ifdef FUNC_USE_BOOST
  std::function<adVar<TOUT,2>(adVar<TIN,2>)> mp_boost_func;
#endif

public:
  QuadraticTaylorTable(FunctionContainer<TIN,TOUT> *func_container, LookupTableParameters<TIN> par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<TIN,TOUT,3,GT>(jsonStats.empty() ? // use the default move constructor for MetaTable (probably not elided...)
      std::move(MetaTable<TIN,TOUT,3,GT>(func_container, par)) :
      std::move(MetaTable<TIN,TOUT,3,GT>(jsonStats, classname, func_container)))
  {
#ifndef FUNC_USE_BOOST
    static_assert(sizeof(TIN)!=sizeof(TIN), "Cannot generate a QuadraticTaylorTable without Boost version 1.71.0 or newer");
#else
    if(!jsonStats.empty())
      return; // all our work is already done
    using boost::math::differentiation::make_fvar;

    /* Base class default variables */
    m_name = classname;
    m_order = 3;
    m_numTableEntries = m_numIntervals;
    m_dataSize = static_cast<unsigned>(sizeof(m_table[0]) * (m_numTableEntries));

    if(func_container->autodiff1_func == NULL)
      throw std::invalid_argument("QuadraticTaylorTable needs the 2nd derivative but this is not defined");

    mp_boost_func = func_container->autodiff2_func;

    /* Allocate and set table */
    m_grid.reset(new TIN[m_numTableEntries]);
    m_table.reset(new polynomial<TOUT,3>[m_numTableEntries]);
    for (unsigned int ii=0;ii<m_numTableEntries;++ii) {
      // nonuniform grids are not supported for Taylor tables
      TIN xgrid = m_minArg + ii*m_stepSize;
      TIN xcenter = xgrid + 0.5*m_stepSize;
      m_grid[ii] = xgrid;
      auto const derivs = (mp_boost_func)(make_fvar<TIN,2>(xcenter));
      auto const d2 = derivs.derivative(2);
      auto const d1 = derivs.derivative(1);
      auto const d0 = derivs.derivative(0);
      auto const h  = m_stepSize;
      m_table[ii].coefs[2] = 0.5*h*h*d2;
      m_table[ii].coefs[1] = h*(d1 - 0.5*h*d2);
      m_table[ii].coefs[0] = d0 - 0.5*h*d1 + 0.125*h*h*d2;
    }
#endif
  }

#ifdef FUNC_USE_BOOST
  std::function<adVar<TOUT,2>(adVar<TIN,2>)> boost_function(){ return mp_boost_func; }
#endif
};

template <typename TIN, typename TOUT, GridTypes GT>
const std::string QuadraticTaylorTable<TIN,TOUT,GT>::classname = grid_type_to_string<GT>() + "QuadraticTaylorTable";

template <typename TIN, typename TOUT=TIN>
using UniformQuadraticTaylorTable = QuadraticTaylorTable<TIN,TOUT,GridTypes::UNIFORM>;
} // namespace func

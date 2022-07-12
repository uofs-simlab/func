/*
  Linear Taylor LUT with uniform sampling

  Usage example:
    LinearTaylorTable look(&function,0,10,0.0001);
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

template <typename TIN, typename TOUT=TIN, GridTypes GT=UNIFORM>
class LinearTaylorTable final : public MetaTable<TIN,TOUT,2,GT>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,2,GT);

  static const std::string classname;
#ifdef FUNC_USE_BOOST
  std::function<adVar<TOUT,1>(adVar<TIN,1>)> mp_boost_func;
#endif

public:
  LinearTaylorTable(FunctionContainer<TIN,TOUT> *func_container, LookupTableParameters<TIN> par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<TIN,TOUT,2,GT>(jsonStats.empty() ? // use the default move constructor for MetaTable (probably not elided...)
      std::move(MetaTable<TIN,TOUT,2,GT>(func_container, par)) :
      std::move(MetaTable<TIN,TOUT,2,GT>(jsonStats, classname, func_container)))
  {
#ifndef FUNC_USE_BOOST
    static_assert(sizeof(TIN)!=sizeof(TIN), "Cannot generate a LinearTaylorTable without Boost version 1.71.0 or newer");
#else
    using boost::math::differentiation::make_fvar;

    /* Base class variables */
    m_name = classname;
    m_order = 2;
    m_numTableEntries = m_numIntervals;
    m_dataSize = (unsigned) sizeof(m_table[0]) * (m_numTableEntries);

    if(func_container->autodiff1_func == NULL)
      throw std::invalid_argument("LinearTaylorTable needs the 1st derivative but this is not defined");

    mp_boost_func = func_container->autodiff1_func;

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,2>[m_numTableEntries]);
    for (unsigned int ii=0; ii<m_numIntervals; ++ii) {
      // nonuniform grids are not supported for Taylor tables
      TIN xgrid = m_minArg + ii*m_stepSize;
      TIN xcenter = xgrid + 0.5*m_stepSize;
      m_grid[ii] = xgrid;
      // get every derivative up to the first
      auto const derivs = (mp_boost_func)(make_fvar<TIN,1>(xcenter));
      auto const d1 = derivs.derivative(1);
      auto const d0 = derivs.derivative(0);
      auto const h  = m_stepSize;
      m_table[ii].coefs[1] = h*d1;
      m_table[ii].coefs[0] = d0 - 0.5*h*d1;
    }
#endif
  }

  // operator() comes from MetaTable
#ifdef FUNC_USE_BOOST
  std::function<adVar<TOUT,1>(adVar<TIN,1>)> boost_function(){ return mp_boost_func; }
#endif
};

template <typename TIN, typename TOUT, GridTypes GT>
const std::string LinearTaylorTable<TIN,TOUT,GT>::classname = grid_type_to_string<GT>() + "LinearTaylorTable";

// define friendlier names
template <typename TIN, typename TOUT=TIN>
using UniformLinearTaylorTable = LinearTaylorTable<TIN,TOUT,UNIFORM>;
} // namespace func

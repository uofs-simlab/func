/*
  Cubic Taylor LUT

  Usage example:
    CubicTaylorTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "MetaTable.hpp"
#include "FunctionContainer.hpp"
#include "config.hpp" // FUNC_USE_BOOST


template <typename TIN, typename TOUT=TIN, GridTypes GT=UNIFORM>
class CubicTaylorTable final : public MetaTable<TIN,TOUT,4,TAYLOR,GT>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,4,TAYLOR,GT);

#ifdef FUNC_USE_BOOST
  std::function<adVar<TOUT,3>(adVar<TOUT,3>)> mp_boost_func;
#endif

public:
  CubicTaylorTable(FunctionContainer<TIN,TOUT> *func_container, LookupTableParameters<TIN> par) :
    MetaTable<TIN,TOUT,4,TAYLOR,GT>(func_container, par)
  {
#ifndef FUNC_USE_BOOST
    static_assert(sizeof(TIN)!=sizeof(TIN), "Cannot generate a CubicTaylorTable without Boost version 1.71.0 or newer");
#else
    using boost::math::differentiation::make_fvar;

    /* Base class default variables */
    m_name = grid_type_to_string<GT>() + "CubicTaylorTable";
    m_order = 4;
    m_numTableEntries = m_numIntervals;
    m_dataSize = (unsigned) sizeof(m_table[0]) * (m_numTableEntries);

    if(func_container->autodiff3_func == NULL)
      throw std::invalid_argument("CubicTaylorTable needs the 3rd derivative but this is not defined");

    mp_boost_func = func_container->autodiff3_func;

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,4>[m_numTableEntries]);
    for (unsigned int ii=0;ii<m_numIntervals;++ii) {
      // nonuniform grids are not supported for Taylor tables
      TIN xgrid = m_minArg + ii*m_stepSize;
      TIN xcenter = xgrid + 0.5*m_stepSize;
      m_grid[ii] = xgrid;
      auto const derivs = (mp_boost_func)(make_fvar<TIN,3>(xcenter));
      auto const d3 = derivs.derivative(3);
      auto const d2 = derivs.derivative(2);
      auto const d1 = derivs.derivative(1);
      auto const d0 = derivs.derivative(0);
      auto const h  = m_stepSize;
      m_table[ii].coefs[3] = h*h*h*d3/6;
      m_table[ii].coefs[2] = h*h*(0.5*d2 - 0.25*h*d3);
      m_table[ii].coefs[1] = h*(d1 - 0.5*h*d2 + 0.125*h*h*d3);
      m_table[ii].coefs[0] = d0 - 0.5*h*d1 + 0.125*h*h*d2 - h*h*h*d3/48;
    }
#endif
  }

  /* build this table from a file. Everything other than m_table is built by MetaTable */
  CubicTaylorTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    MetaTable<TIN,TOUT,4,TAYLOR,GT>(func_container, filename,
        grid_type_to_string<GT>() + "CubicTaylorTable") {}
  // operator() comes from MetaTable

#ifdef FUNC_USE_BOOST
  std::function<adVar<TOUT,3>(adVar<TOUT,3>)> boost_function(){ return mp_boost_func; }
#endif
};

// define friendlier names
template <typename TIN, typename TOUT=TIN>
using UniformCubicTaylorTable = CubicTaylorTable<TIN,TOUT,UNIFORM>;

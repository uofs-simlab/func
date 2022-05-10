/*
  Cubic Interpolation LUT with precomputed polynomial coefficients

  Usage example:
    CubicHermiteTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "MetaTable.hpp"
#include "config.hpp" // FUNC_USE_BOOST

template <typename TIN, typename TOUT=TIN, GridTypes GT=UNIFORM>
class CubicHermiteTable final : public MetaTable<TIN,TOUT,4,HORNER,GT>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,4,HORNER,GT);

  FUNC_REGISTER_LUT(CubicHermiteTable);

#ifdef FUNC_USE_BOOST
  std::function<adVar<TOUT,1>(adVar<TIN,1>)> mp_boost_func;
#endif

public:
  CubicHermiteTable(FunctionContainer<TIN,TOUT> *func_container, LookupTableParameters<TIN> par) :
      MetaTable<TIN,TOUT,4,HORNER,GT>(func_container, par)
  {
#ifndef FUNC_USE_BOOST
    static_assert(sizeof(TIN)!=sizeof(TIN), "Cannot generate a CubicHermiteTable without Boost version 1.71.0 or newer");
#else
    using boost::math::differentiation::make_fvar;
    /* Base class default variables */
    m_name = grid_type_to_string<GT>() + "CubicHermiteTable";
    m_order = 4;
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = (unsigned) sizeof(m_table[0]) * (m_numTableEntries);

    if(func_container->autodiff1_func == NULL)
      throw std::invalid_argument("CubicHermiteTable needs the 1st derivative but this is not defined");
    mp_boost_func = func_container->autodiff1_func;

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,4>[m_numTableEntries]);
    for (int ii=0; ii<m_numIntervals; ++ii) {
      TIN x;
      TIN h = m_stepSize;
      // (possibly) transform the uniform grid into a nonuniform grid
      if (GT == UNIFORM)
        x = m_minArg + ii*m_stepSize;
      else{
        x = m_transferFunction.g(m_minArg + ii*m_stepSize);
        h = m_transferFunction.g(m_minArg + (ii+1)*m_stepSize) - x;
      }
      m_grid[ii] = x;

      const auto derivs0 = (mp_boost_func)(make_fvar<TIN,1>(x));
      const TOUT y0    = derivs0.derivative(0);
      const TOUT m0    = derivs0.derivative(1);
      const auto derivs1 = (mp_boost_func)(make_fvar<TIN,1>(x+h));
      const TOUT y1    = derivs1.derivative(0);
      const TOUT m1    = derivs1.derivative(1);

      m_table[ii].coefs[0] = y0;
      m_table[ii].coefs[1] = h*m0;
      m_table[ii].coefs[2] = -3*y0+3*y1-(2*m0+m1)*h;
      m_table[ii].coefs[3] = 2*y0-2*y1+(m0+m1)*h;
    }
#endif
  }

  /* build this table from a file. Everything other than m_table is built by MetaTable */
  CubicHermiteTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    MetaTable<TIN,TOUT,4,HORNER,GT>(func_container, filename,
        grid_type_to_string<GT>() + "CubicHermiteTable") {}
  // operator() comes straight from the MetaTable

#ifdef FUNC_USE_BOOST
  std::function<adVar<TOUT,1>(adVar<TOUT,1>)> boost_function(){ return mp_boost_func; }
#endif
};

// define friendlier names
template <typename TIN, typename TOUT=TIN>
using UniformCubicHermiteTable = CubicHermiteTable<TIN,TOUT,UNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformCubicHermiteTable = CubicHermiteTable<TIN,TOUT,NONUNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformPseudoCubicHermiteTable = CubicHermiteTable<TIN,TOUT,NONUNIFORM_PSEUDO>;

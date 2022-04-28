/*
  Cubic Interpolation LUT with uniform sampling (precomputed coefficients)

  Usage example:
    UniformCubicHermiteTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "MetaTable.hpp"
#include "config.hpp"

#ifndef FUNC_USE_BOOST_AUTODIFF
#error "UniformCubicHermiteTable needs boost version >= 1.71"
#endif

template <typename TIN, typename TOUT = TIN>
class UniformCubicHermiteTable final : public MetaTable<TIN,TOUT,4,HORNER>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_UNIFORM_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,4,HORNER);

  FUNC_REGISTER_LUT(UniformCubicHermiteTable);

  std::function<adVar<TOUT,1>(adVar<TIN,1>)> mp_boost_func;

public:
  UniformCubicHermiteTable(FunctionContainer<TIN,TOUT> *func_container, UniformLookupTableParameters<TIN> par) :
      MetaTable<TIN,TOUT,4,HORNER>(func_container, par)
  {
    using boost::math::differentiation::make_fvar;
    /* Base class default variables */
    m_name = "UniformCubicHermiteTable";
    m_order = 4;
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = (unsigned) sizeof(m_table[0]) * (m_numTableEntries);

    __IS_NULL(func_container->autodiff1_func);
    mp_boost_func = func_container->autodiff1_func;

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,4>[m_numTableEntries]);
    for (int ii=0; ii<m_numIntervals; ++ii) {
      const TIN x = m_minArg + ii*m_stepSize;
      m_grid[ii] = x;

      const auto derivs0 = (mp_boost_func)(make_fvar<TIN,1>(x));
      const TOUT y0    = derivs0.derivative(0);
      const TOUT m0    = derivs0.derivative(1);
      const auto derivs1 = (mp_boost_func)(make_fvar<TIN,1>(x+m_stepSize));
      const TOUT y1    = derivs1.derivative(0);
      const TOUT m1    = derivs1.derivative(1);

      m_table[ii].coefs[0] = y0;
      m_table[ii].coefs[1] = m_stepSize*m0;
      m_table[ii].coefs[2] = -3*y0+3*y1-(2*m0+m1)*m_stepSize;
      m_table[ii].coefs[3] = 2*y0-2*y1+(m0+m1)*m_stepSize;
    }
  }

  /* build this table from a file. Everything other than m_table is built by MetaTable */
  UniformCubicHermiteTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    MetaTable<TIN,TOUT,4,HORNER>(func_container, filename, "UniformCubicHermiteTable") {}
  // operator() comes straight from the MetaTable

  std::function<adVar<TOUT,1>(adVar<TOUT,1>)> boost_function(){ return mp_boost_func; }
};

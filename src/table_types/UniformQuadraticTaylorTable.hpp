/*
  Quadratic Taylor LUT with uniform sampling

  Usage example:
    UniformQuadraticTaylorTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "MetaTable.hpp"
#include "config.hpp"

#ifndef FUNC_USE_BOOST_AUTODIFF
#error "UniformQuadraticTaylorTable needs boost version >= 1.71"
#endif

template <typename TIN, typename TOUT = TIN>
class UniformQuadraticTaylorTable final : public MetaTable<TIN,TOUT,3,TAYLOR>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_UNIFORM_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,3,TAYLOR);
  FUNC_REGISTER_LUT(UniformQuadraticTaylorTable);

  std::function<adVar<TOUT,2>(adVar<TIN,2>)> mp_boost_func;
public:
  UniformQuadraticTaylorTable(FunctionContainer<TIN,TOUT> *func_container, UniformLookupTableParameters<TIN> par) :
    MetaTable<TIN,TOUT,3,TAYLOR>(func_container, par)
  {
    using boost::math::differentiation::make_fvar;

    /* Base class default variables */
    m_name = "UniformQuadraticTaylorTable";
    m_order = 3;
    m_numTableEntries = m_numIntervals;
    m_dataSize = (unsigned) sizeof(m_table[0]) * (m_numTableEntries);

    __IS_NULL(func_container->autodiff2_func);
    mp_boost_func = func_container->autodiff2_func;

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,3>[m_numTableEntries]);
    for (int ii=0;ii<m_numIntervals;++ii) {
      TIN x = (m_minArg + ii*m_stepSize);
      m_grid[ii] = x;
      auto const derivs = (mp_boost_func)(make_fvar<TIN,2>(x));
      m_table[ii].coefs[0] = derivs.derivative(0);
      m_table[ii].coefs[1] = derivs.derivative(1);
      m_table[ii].coefs[2] = derivs.derivative(2)/2;
    }
  }

  UniformQuadraticTaylorTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    MetaTable<TIN,TOUT,3,TAYLOR>(func_container, filename, "UniformQuadraticTaylorTable") {}

  std::function<adVar<TOUT,2>(adVar<TIN,2>)> boost_function(){ return mp_boost_func; }
};

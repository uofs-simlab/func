/*
  Linear Taylor LUT with uniform sampling

  Usage example:
    UniformLinearTaylorTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "MetaTable.hpp"
#include "config.hpp"

#ifndef FUNC_USE_BOOST_AUTODIFF
#error "UniformLinearTaylorTable needs boost version >= 1.71"
#endif

template <typename TIN, typename TOUT = TIN>
class UniformLinearTaylorTable final : public MetaTable<TIN,TOUT,2,TAYLOR>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_UNIFORM_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,2,TAYLOR);
  FUNC_REGISTER_LUT(UniformLinearTaylorTable);

  std::function<adVar<TOUT,1>(adVar<TIN,1>)> mp_boost_func;
public:
  UniformLinearTaylorTable(FunctionContainer<TIN,TOUT> *func_container, UniformLookupTableParameters<TIN> par) :
    MetaTable<TIN,TOUT,2,TAYLOR>(func_container, par)
  {
    using boost::math::differentiation::make_fvar;

    /* Base class variables */
    m_name = "UniformLinearTaylorTable";
    m_order = 2;
    m_numTableEntries = m_numIntervals;
    m_dataSize = (unsigned) sizeof(m_table[0]) * (m_numTableEntries);

    __IS_NULL(func_container->autodiff1_func);
    mp_boost_func = func_container->autodiff1_func;

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,2>[m_numTableEntries]);
    for (int ii=0; ii<m_numIntervals; ++ii) {
      const TIN x = m_minArg + ii*m_stepSize;
      m_grid[ii]     = x;
      // get every derivative up to the first
      auto const derivs = (mp_boost_func)(make_fvar<TIN,1>(x));
      m_table[ii].coefs[0] = derivs.derivative(0);
      m_table[ii].coefs[1] = derivs.derivative(1);
    }
  }

  /* build this table from a file. Everything other than m_table is built by MetaTable */
  UniformLinearTaylorTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    MetaTable<TIN,TOUT,2,TAYLOR>(func_container, filename, "UniformLinearTaylorTable") {}
  // operator() comes from MetaTable

  std::function<adVar<TOUT,1>(adVar<TIN,1>)> boost_function(){ return mp_boost_func; }
};

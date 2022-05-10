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
#include "config.hpp" // FUNC_USE_BOOST


template <typename TIN, typename TOUT=TIN, GridTypes GT=UNIFORM>
class LinearTaylorTable final : public MetaTable<TIN,TOUT,2,TAYLOR,GT>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,2,TAYLOR,GT);
  FUNC_REGISTER_LUT(LinearTaylorTable);

  std::function<adVar<TOUT,1>(adVar<TIN,1>)> mp_boost_func;
public:
  LinearTaylorTable(FunctionContainer<TIN,TOUT> *func_container, LookupTableParameters<TIN> par) :
    MetaTable<TIN,TOUT,2,TAYLOR,GT>(func_container, par)
  {
#ifndef FUNC_USE_BOOST
    static_assert(sizeof(TIN)!=sizeof(TIN), "Cannot generate a LinearTaylorTable without Boost version 1.71.0 or newer");
#else
    using boost::math::differentiation::make_fvar;

    /* Base class variables */
    m_name = grid_type_to_string<GT>() + "LinearTaylorTable";
    m_order = 2;
    m_numTableEntries = m_numIntervals;
    m_dataSize = (unsigned) sizeof(m_table[0]) * (m_numTableEntries);

    if(func_container->autodiff1_func == NULL)
      throw std::invalid_argument("LinearTaylorTable needs the 1st derivative but this is not defined");

    mp_boost_func = func_container->autodiff1_func;

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,2>[m_numTableEntries]);
    for (int ii=0; ii<m_numIntervals; ++ii) {
      // nonuniform grids are not supported for Taylor tables
      TIN x = m_minArg + ii*m_stepSize;

      m_grid[ii] = x;
      // get every derivative up to the first
      auto const derivs = (mp_boost_func)(make_fvar<TIN,1>(x));
      m_table[ii].coefs[0] = derivs.derivative(0);
      m_table[ii].coefs[1] = derivs.derivative(1);
    }
#endif
  }

  /* build this table from a file. Everything other than m_table is built by MetaTable */
  LinearTaylorTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    MetaTable<TIN,TOUT,2,TAYLOR,GT>(func_container, filename,
        grid_type_to_string<GT>() + "LinearTaylorTable") {}
  // operator() comes from MetaTable

  std::function<adVar<TOUT,1>(adVar<TIN,1>)> boost_function(){ return mp_boost_func; }
};

// define friendlier names
template <typename TIN, typename TOUT=TIN>
using UniformLinearTaylorTable = LinearTaylorTable<TIN,TOUT,UNIFORM>;

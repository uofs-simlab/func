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

template <typename TIN, typename TOUT, unsigned int N, GridTypes GT=UNIFORM>
class TaylorTable final : public MetaTable<TIN,TOUT,N+1,GT>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,N+1,GT);

  static TOUT constexpr fact[] = {1,1,2,6,24,120,720,5040};

#ifdef FUNC_USE_BOOST
  std::function<adVar<TOUT,N>(adVar<TOUT,N>)> mp_boost_func;
#endif

  std::string degree_to_string(){
    switch(N){
      case 1:
        return "Linear";
      case 2:
        return "Quadratic";
      case 3:
        return "Cubic";
      default:
        return "Degree" + std::to_string(N);
    }
  }

public:
  TaylorTable(FunctionContainer<TIN,TOUT> *func_container, LookupTableParameters<TIN> par) :
    MetaTable<TIN,TOUT,N+1,GT>(func_container, par)
  {
#ifndef FUNC_USE_BOOST
    static_assert(sizeof(TIN)!=sizeof(TIN), "Cannot generate a TaylorTable without Boost version 1.71.0 or newer");
#else
    using boost::math::differentiation::make_fvar;

    /* Base class default variables */
    m_name = grid_type_to_string<GT>()+degree_to_string()+"TaylorTable";
    m_order = N+1;
    m_numTableEntries = m_numIntervals;
    m_dataSize = (unsigned) sizeof(m_table[0]) * (m_numTableEntries);

    if(func_container->template get_nth_func<N>() == nullptr)
      throw std::invalid_argument(m_name+" needs the "+std::to_string(N)+"th derivative but this is not defined");
    mp_boost_func = func_container->template get_nth_func<N>();

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,N+1>[m_numTableEntries]);
    for (unsigned int ii=0;ii<m_numIntervals;++ii) {
      // nonuniform grids are not supported for Taylor tables
      TIN x = m_minArg + ii*m_stepSize;

      m_grid[ii] = x;
      auto const derivs = (mp_boost_func)(make_fvar<TIN,N>(x));
      for(unsigned int k=0; k<N+1; k++)
        m_table[ii].coefs[k] = derivs.derivative(k)/(fact[k]);     
    }
#endif
  }

  /* build this table from a file. Everything other than m_table is built by MetaTable */
  TaylorTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    MetaTable<TIN,TOUT,N+1,GT>(func_container, filename,
        grid_type_to_string<GT>() + degree_to_string() + "TaylorTable") {}

  TOUT operator()(TIN x) override
  {
    TOUT dx;
    unsigned int x0;
    // nondimensionalized x position
    dx = (x-m_minArg);
    TOUT x0r = dx/m_stepSize+0.5;
    // index of previous table entry
    x0 = (unsigned) x0r;
    dx -= x0*m_stepSize; 

    // general degree horners method, evaluated from the inside out.
    // TODO time to see if the added generality leads to any slowdown
    TOUT sum = 0;
    for (unsigned int k=N; k>0; k--)
      sum = dx*(m_table[x0].coefs[k] + sum);
    return m_table[x0].coefs[0]+sum;
  }

#ifdef FUNC_USE_BOOST
  std::function<adVar<TOUT,N>(adVar<TOUT,N>)> boost_function(){ return mp_boost_func; }
#endif
};

// define friendlier names
template <typename TIN, typename TOUT=TIN>
using UniformLinearTaylorTable = TaylorTable<TIN,TOUT,1,UNIFORM>;
template <typename TIN, typename TOUT=TIN>
using UniformQuadraticTaylorTable = TaylorTable<TIN,TOUT,2,UNIFORM>;
template <typename TIN, typename TOUT=TIN>
using UniformCubicTaylorTable = TaylorTable<TIN,TOUT,3,UNIFORM>;

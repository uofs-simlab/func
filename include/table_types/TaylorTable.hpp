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

  static const std::string degree_to_string() {
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
  // build the LUT from scratch or look in filename for an existing LUT
  TaylorTable(FunctionContainer<TIN,TOUT> *func_container, LookupTableParameters<TIN> par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<TIN,TOUT,N+1,GT>(jsonStats.empty() ? // use the default move constructor for MetaTable (probably not elided...)
      std::move(MetaTable<TIN,TOUT,N+1,GT>(func_container, par)) :
      std::move(MetaTable<TIN,TOUT,N+1,GT>(jsonStats, classname, func_container)))
  {
#ifndef FUNC_USE_BOOST
    static_assert(sizeof(TIN)!=sizeof(TIN), "Cannot generate a TaylorTable without Boost version 1.71.0 or newer");
#else
    if(!jsonStats.empty())
      return; // all our work is already done

    using boost::math::differentiation::make_fvar;

    /* Base class default variables */
    m_name = classname;
    m_order = N+1;
    m_numTableEntries = m_numIntervals;
    m_dataSize = static_cast<unsigned>(sizeof(m_table[0]) * (m_numTableEntries));

    if(func_container->template get_nth_func<N>() == nullptr)
      throw std::invalid_argument(m_name+" needs the "+std::to_string(N)+"th derivative but this is not defined");
    mp_boost_func = func_container->template get_nth_func<N>();

    /* Allocate and set table */
    m_grid.reset(new TIN[m_numTableEntries]);
    m_table.reset(new polynomial<TOUT,N+1>[m_numTableEntries]);
    for (unsigned int ii=0; ii<m_numTableEntries; ++ii) {
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

      /* TODO there is certainly a nice enough closed form solution for the coefs of 
       * an arbitrary degree shifted Taylor approx. but we can figure that out another day */
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
#endif
  }

#ifdef FUNC_USE_BOOST
  std::function<adVar<TOUT,N>(adVar<TOUT,N>)> boost_function(){ return mp_boost_func; }
#endif
};

template <typename TIN, typename TOUT, unsigned int N, GridTypes GT>
const std::string TaylorTable<TIN,TOUT,N,GT>::classname = grid_type_to_string<GT>()+TaylorTable<TIN,TOUT,N,GT>::degree_to_string()+"TaylorTable";

/* define friendlier names */
template <typename TIN, typename TOUT=TIN>
using UniformLinearTaylorTable = TaylorTable<TIN,TOUT,1,GridTypes::UNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformLinearTaylorTable = TaylorTable<TIN,TOUT,1,GridTypes::NONUNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformPseudoLinearTaylorTable = TaylorTable<TIN,TOUT,1,GridTypes::NONUNIFORM_PSEUDO>;

template <typename TIN, typename TOUT=TIN>
using UniformQuadraticTaylorTable = TaylorTable<TIN,TOUT,2,GridTypes::UNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformQuadraticTaylorTable = TaylorTable<TIN,TOUT,2,GridTypes::NONUNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformPseudoQuadraticTaylorTable = TaylorTable<TIN,TOUT,2,GridTypes::NONUNIFORM_PSEUDO>;

template <typename TIN, typename TOUT=TIN>
using UniformCubicTaylorTable = TaylorTable<TIN,TOUT,3,GridTypes::UNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformCubicTaylorTable = TaylorTable<TIN,TOUT,3,GridTypes::NONUNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformPseudoCubicTaylorTable = TaylorTable<TIN,TOUT,3,GridTypes::NONUNIFORM_PSEUDO>;
} // namespace func

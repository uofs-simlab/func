/*
  Linear Interpolation LUT with uniform sampling (precomputed coefficients)

  Usage example:
    UniformEqSpaceInterpTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "MetaTable.hpp"

namespace func {

template <unsigned int N, typename TIN, typename TOUT, GridTypes GT=GridTypes::UNIFORM>
class EqSpaceInterpTable final : public MetaTable<N+1,TIN,TOUT,GT>
{
  INHERIT_META(N+1,TIN,TOUT,GT);
public:
  // build the LUT from scratch or look in filename for an existing LUT
  EqSpaceInterpTable(const FunctionContainer<TIN,TOUT>& fun_container, const LookupTableParameters<TIN>& par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<N+1,TIN,TOUT,GT>(jsonStats.empty() ? // use the default move constructor for MetaTable (probably not elided...)
      std::move(MetaTable<N+1,TIN,TOUT,GT>(fun_container, par)) :
      std::move(MetaTable<N+1,TIN,TOUT,GT>(jsonStats)))
  {
    if(!jsonStats.empty())
      return; // all our work is already done

    /* Base class default variables */
    m_name = grid_type_to_string<GT>() + "EqSpaceInterpTable<" + std::to_string(N) + ">";
    m_order = N+1;
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = static_cast<unsigned>(sizeof(m_table[0]) * (m_numTableEntries));

    auto fun = fun_container.standard_fun;

    /* Allocate and set table */
    m_grid.reset(new TIN[m_numTableEntries]);
    m_table.reset(new polynomial<TOUT,N+1>[m_numTableEntries]);
    FUNC_BUILDPAR
    for (unsigned int ii=0; ii<m_numTableEntries-1; ++ii) {
      TIN x;
      TIN h = m_stepSize;
      // (possibly) transform the uniform grid into a nonuniform grid
      if (GT == GridTypes::UNIFORM)
        x = m_minArg + ii*m_stepSize;
      else{
        x = m_transferFunction(m_minArg + ii*m_stepSize);
        h = m_transferFunction(m_minArg + (ii+1)*m_stepSize) - x;
      }

      m_grid[ii] = x;
      switch(N){
        case 0:
          m_table[ii].coefs[0] = fun(x + 0.5*h);
          break;
        case 1:
          {
          m_table[ii].coefs[0] = fun(x);
          m_table[ii].coefs[1] = fun(x+h) - m_table[ii].coefs[0];
          break;
          }
        case 2:
          {
          const TOUT y0 = fun(x);
          const TOUT y1 = fun(x+h/2);
          const TOUT y2 = fun(x+h);
          m_table[ii].coefs[0] = y0;
          m_table[ii].coefs[1] = -3*y0+4*y1-y2;
          m_table[ii].coefs[2] = 2*y0+-4*y1+2*y2;
          break;
          }
        case 3:
          {
          const TOUT y0 = fun(x);
          const TOUT y1 = fun(x+h/3);
          const TOUT y2 = fun(x+2*h/3);
          const TOUT y3 = fun(x+h);
          m_table[ii].coefs[0] = y0;
          m_table[ii].coefs[1] = -11*y0/2+9*y1-9*y2/2+y3;
          m_table[ii].coefs[2] = 9*y0-45*y1/2+18*y2-9*y3/2;
          m_table[ii].coefs[3] = -9*y0/2+27*y1/2-27*y2/2+9*y3/2;
          break;
          }
        default: { throw std::invalid_argument("Broken switch case in func::TaylorTable"); }
      }
      if(GT == GridTypes::NONUNIFORM){
        auto p = m_table[ii];
        for(unsigned int k=0; k<N+1; k++)
          m_table[ii].coefs[k] = polynomial_diff(p,-x/h,k)/std::pow(h,k)/boost::math::factorial<double>(k);
      }
    }
    // special case to make lut(tableMaxArg) work
    m_grid[m_numTableEntries-1] = m_tableMaxArg;
    m_table[m_numTableEntries-1].coefs[0] = fun(m_tableMaxArg);
    for (unsigned int k=1; k<N+1; k++)
      m_table[m_numTableEntries-1].coefs[k] = 0;
  }

  // operator() is in MetaTable
};

// define friendlier names
template <unsigned int N, typename TIN, typename TOUT=TIN>
using UniformEqSpaceInterpTable = EqSpaceInterpTable<N,TIN,TOUT,GridTypes::UNIFORM>;
template <unsigned int N, typename TIN, typename TOUT=TIN>
using NonUniformEqSpaceInterpTable = EqSpaceInterpTable<N,TIN,TOUT,GridTypes::NONUNIFORM>;

} // namespace func

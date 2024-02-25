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
  EqSpaceInterpTable() = default;
  EqSpaceInterpTable(const MetaTable<N+1,TIN,TOUT,GT>& L): MetaTable<N+1,TIN,TOUT,GT>(L) {}

  // build the LUT from scratch or look in filename for an existing LUT
  EqSpaceInterpTable(const FunctionContainer<TIN,TOUT>& fun_container, const LookupTableParameters<TIN>& par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<N+1,TIN,TOUT,GT>(fun_container, par, jsonStats)
  {
    if(!jsonStats.empty())
      return; // all our work is already done

    /* Base class default variables */
    m_name = grid_type_to_string<GT>() + "EqSpaceInterpTable<" + std::to_string(N) + ">";
    m_order = N+1;
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = m_numTableEntries * sizeof(m_table[0]);

    auto fun = fun_container.standard_fun;
    if(fun == nullptr)
      throw std::invalid_argument("Error in func::EqSpaceInterpTable: Given an invalid FunctionContainer");

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,N+1>[m_numTableEntries]);
    FUNC_BUILDPAR
    for(unsigned int ii=0; ii<m_numTableEntries-1; ++ii) {
      TIN x;
      TIN h = m_stepSize;
      // (possibly) transform the uniform grid into a nonuniform grid
      FUNC_IF_CONSTEXPR(GT == GridTypes::UNIFORM)
        x = m_minArg + ii*m_stepSize;
      else{
        x = m_transferFunction(m_minArg + ii*m_stepSize);
        h = m_transferFunction(m_minArg + (ii+1)*m_stepSize) - x;
      }

      /* how many equally spaced nodes do we use? */
      switch(N){
        case 0:
          m_table[ii].coefs[0] = fun(x + h/2);
          break;
        case 1:
          {
          m_table[ii].coefs[0] = fun(x);
          m_table[ii].coefs[1] = fun(x+h) - m_table[ii].coefs[0];
          break;
          }
        case 2:
          {
          TIN a,b,c;
          const TOUT y0 = fun(x);
          const TOUT y1 = fun(x+h/2);
          const TOUT y2 = fun(x+h);
          m_table[ii].coefs[0] = y0;
          a = -3, b = 4,  c = -1; m_table[ii].coefs[1] = a*y0+b*y1+c*y2;
          a = 2,  b = -4, c = 2;  m_table[ii].coefs[2] = a*y0+b*y1+c*y2;
          break;
          }
        case 3:
          {
          TIN a,b,c,d;
          const TOUT y0 = fun(x);
          const TOUT y1 = fun(x+h/3);
          const TOUT y2 = fun(x+2*h/3);
          const TOUT y3 = fun(x+h);
          m_table[ii].coefs[0] = y0;
          a = -11/2.0, b = 9,       c = -9/2.0,  d = 1.0;    m_table[ii].coefs[1] = a*y0+b*y1+c*y2+  y3;
          a = 9,       b = -45/2.0, c = 18,      d = -9/2.0; m_table[ii].coefs[2] = a*y0+b*y1+c*y2+d*y3;
          a = -9/2.0,  b = 27/2.0,  c = -27/2.0, d = 9/2.0;  m_table[ii].coefs[3] = a*y0+b*y1+c*y2+d*y3;
          break;
          }
        default: { throw std::invalid_argument(std::string("EqSpaceInterpTables<N> only support N=0,1,2,3 but given N=") + std::to_string(N)); }
      }

      /* TODO This formula is too unstable for this table type as given in this form when N>2 and h is small. */
      FUNC_IF_CONSTEXPR(GT == GridTypes::NONUNIFORM){
        auto p = m_table[ii];
        for(unsigned int k=0; k<N+1; k++)
          m_table[ii].coefs[k] = polynomial_diff(p,-x/h,k)/static_cast<TIN>(pow(h,k))/static_cast<TIN>(factorial(k));
      }
    }
    /* special case to make lut(tableMaxArg) work. All other coefs are copied
     * from the previous subinterval to get diff() to work predictably */
    m_table[m_numTableEntries-1].coefs[0] = fun(m_tableMaxArg);
    for(unsigned int k=1; k<N+1; k++){
      m_table[m_numTableEntries-1].coefs[k] = m_table[m_numTableEntries-2].coefs[k];
      //m_table[m_numTableEntries-1].coefs[k] = static_cast<TIN>(0)*m_table[m_numTableEntries-2].coefs[k];
    }
  }

  // operator() is in MetaTable
};

// define friendlier names
template <unsigned int N, typename TIN, typename TOUT=TIN>
using UniformEqSpaceInterpTable = EqSpaceInterpTable<N,TIN,TOUT,GridTypes::UNIFORM>;
template <unsigned int N, typename TIN, typename TOUT=TIN>
using NonUniformEqSpaceInterpTable = EqSpaceInterpTable<N,TIN,TOUT,GridTypes::NONUNIFORM>;

} // namespace func

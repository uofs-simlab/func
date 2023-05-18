/*
  Linear Interpolation LUT with uniform sampling (precomputed coefficients)

  Usage example:
    LinearInterpolationTable look(&function,0,10,0.0001);
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

template <typename TIN, typename TOUT, unsigned int N, GridTypes GT=GridTypes::UNIFORM>
class LinearInterpolationTable final : public MetaTable<TIN,TOUT,N+1,GT>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,N+1,GT);

  static const std::string classname;
public:
  // build the LUT from scratch or look in filename for an existing LUT
  LinearInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, LookupTableParameters<TIN> par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<TIN,TOUT,N+1,GT>(jsonStats.empty() ? // use the default move constructor for MetaTable (probably not elided...)
      std::move(MetaTable<TIN,TOUT,N+1,GT>(func_container, par)) :
      std::move(MetaTable<TIN,TOUT,N+1,GT>(jsonStats, classname, func_container)))
  {
    if(!jsonStats.empty())
      return; // all our work is already done

    /* Base class default variables */
    m_name = classname;
    m_order = N+1;
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = static_cast<unsigned>(sizeof(m_table[0]) * (m_numTableEntries));

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
        x = m_transferFunction.g(m_minArg + ii*m_stepSize);
        h = m_transferFunction.g(m_minArg + (ii+1)*m_stepSize) - x;
      }

      m_grid[ii] = x;
      switch(N){
        case 0:
          m_table[ii].coefs[0] = m_func(x + 0.5*h);
          break;
        case 1:
          {
          m_table[ii].coefs[0] = m_func(x);
          m_table[ii].coefs[1] = m_func(x+h) - m_table[ii].coefs[0];
          break;
          }
        case 2:
          {
          const TOUT y0  = m_func(x);
          const TOUT y1  = m_func(x+h/2);
          const TOUT y2  = m_func(x+h);
          m_table[ii].coefs[0] = y0;
          m_table[ii].coefs[1] = -3*y0+4*y1-y2;
          m_table[ii].coefs[2] = 2*y0+-4*y1+2*y2;
          break;
          }
        case 3:
          {
          const TOUT y0 = m_func(x);
          const TOUT y1 = m_func(x+h/3);
          const TOUT y2 = m_func(x+2*h/3);
          const TOUT y3 = m_func(x+h);
          m_table[ii].coefs[0] = y0;
          m_table[ii].coefs[1] = -11*y0/2+9*y1-9*y2/2+y3;
          m_table[ii].coefs[2] = 9*y0-45*y1/2+18*y2-9*y3/2;
          m_table[ii].coefs[3] = -9*y0/2+27*y1/2-27*y2/2+9*y3/2;
          break;
          }
        default: { throw std::invalid_argument("Broken switch case in func::TaylorTable"); }
      }
    }
    // special case to make lut(tableMaxArg) work
    m_grid[m_numTableEntries-1] = m_tableMaxArg;
    m_table[m_numTableEntries-1].coefs[0] = m_func(m_tableMaxArg);
    for (unsigned int k=1; k<N+1; k++)
      m_table[m_numTableEntries-1].coefs[k] = 0;
  }

  /* build this table from a file. Everything other than m_table is built by MetaTable */
  LinearInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    MetaTable<TIN,TOUT,N+1,GT>(func_container, filename, grid_type_to_string<GT>() + degree_to_string<N>() + "InterpolationTable") {}
  // operator() is in MetaTable
};

template <typename TIN, typename TOUT, unsigned int N, GridTypes GT>
const std::string LinearInterpolationTable<TIN,TOUT,N,GT>::classname = grid_type_to_string<GT>() + degree_to_string<N>() + "InterpolationTable";

// define friendlier names
template <typename TIN, typename TOUT=TIN>
using UniformLinearInterpolationTable = LinearInterpolationTable<TIN,TOUT,1,GridTypes::UNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformLinearInterpolationTable = LinearInterpolationTable<TIN,TOUT,1,GridTypes::NONUNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformPseudoLinearInterpolationTable = LinearInterpolationTable<TIN,TOUT,1,GridTypes::NONUNIFORM_PSEUDO>;

template <typename TIN, typename TOUT=TIN>
using UniformQuadraticInterpolationTable = InterpolationTable<TIN,TOUT,2,GridTypes::UNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformQuadraticInterpolationTable = InterpolationTable<TIN,TOUT,2,GridTypes::NONUNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformPseudoQuadraticInterpolationTable = InterpolationTable<TIN,TOUT,2,GridTypes::NONUNIFORM_PSEUDO>;

template <typename TIN, typename TOUT=TIN>
using UniformCubicInterpolationTable = InterpolationTable<TIN,TOUT,3,GridTypes::UNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformCubicInterpolationTable = InterpolationTable<TIN,TOUT,3,GridTypes::NONUNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformPseudoCubicInterpolationTable = InterpolationTable<TIN,TOUT,3,GridTypes::NONUNIFORM_PSEUDO>;

} // namespace func

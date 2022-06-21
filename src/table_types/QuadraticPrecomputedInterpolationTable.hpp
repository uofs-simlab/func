/*
  Quadratic Interpolation LUT with precomputed coefficients

  Usage example:
    QuadraticPrecomputedInterpolationTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "MetaTable.hpp"

template <typename TIN, typename TOUT=TIN, GridTypes GT=UNIFORM>
class QuadraticPrecomputedInterpolationTable final : public MetaTable<TIN,TOUT,3,HORNER,GT>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,3,HORNER,GT);

public:
  QuadraticPrecomputedInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, LookupTableParameters<TIN> par) :
    MetaTable<TIN,TOUT,3,HORNER,GT>(func_container, par)
  {
    /* Base class default variables */
    m_name = grid_type_to_string<GT>() + "QuadraticPrecomputedInterpolationTable";
    m_order = 3;
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = (unsigned) sizeof(m_table[0]) * (m_numTableEntries);

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,3>[m_numTableEntries]);
    for (int ii=0;ii<m_numIntervals;++ii) {
      TIN x;
      TIN h = m_stepSize;
      // (possibly) transform the uniform grid into a nonuniform grid
      if (GT == UNIFORM)
        x = m_minArg + ii*m_stepSize;
      else{
        x = m_transferFunction.g(m_minArg + ii*m_stepSize);
        h = m_transferFunction.g(m_minArg + (ii+1)*m_stepSize) - x;
      }

      // grid points
      m_grid[ii] = x;
      // polynomial coefficients
      const TOUT y0  = m_func(x);
      const TOUT y1  = m_func(x+h/2);
      const TOUT y2  = m_func(x+h);
      m_table[ii].coefs[0] = y0;
      m_table[ii].coefs[1] = -3*y0+4*y1-y2;
      m_table[ii].coefs[2] = 2*y0+-4*y1+2*y2;
    }
  }

  QuadraticPrecomputedInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    MetaTable<TIN,TOUT,3,HORNER,GT>(func_container, filename,
        grid_type_to_string<GT>() + "QuadraticPrecomputedInterpolationTable") {}
};

template <typename TIN, typename TOUT=TIN>
using UniformQuadraticPrecomputedInterpolationTable = QuadraticPrecomputedInterpolationTable<TIN,TOUT,UNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformQuadraticPrecomputedInterpolationTable = QuadraticPrecomputedInterpolationTable<TIN,TOUT,NONUNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformPseudoQuadraticPrecomputedInterpolationTable = QuadraticPrecomputedInterpolationTable<TIN,TOUT,NONUNIFORM_PSEUDO>;

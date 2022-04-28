/*
  Quadratic Interpolation LUT with uniform sampling (precomputed coefficients)

  Usage example:
    UniformQuadraticPrecomputedInterpolationTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "MetaTable.hpp"

template <typename TIN, typename TOUT = TIN>
class UniformQuadraticPrecomputedInterpolationTable final : public MetaTable<TIN,TOUT,3,HORNER>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_UNIFORM_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,3,HORNER);
  FUNC_REGISTER_LUT(UniformQuadraticPrecomputedInterpolationTable);

public:
  UniformQuadraticPrecomputedInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, UniformLookupTableParameters<TIN> par) :
    MetaTable<TIN,TOUT,3,HORNER>(func_container, par)
  {
    /* Base class default variables */
    m_name = "UniformQuadraticPrecomputedInterpolationTable";
    m_order = 3;
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = (unsigned) sizeof(m_table[0]) * (m_numTableEntries);

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,3>[m_numTableEntries]);
    for (int ii=0;ii<m_numIntervals;++ii) {
      const TIN x = m_minArg + ii*m_stepSize;
      // grid points
      m_grid[ii] = x;
      // polynomial coefficients
      const TOUT y0  = m_func(x);
      const TOUT y1  = m_func(x+m_stepSize/2);
      const TOUT y2  = m_func(x+m_stepSize);
      m_table[ii].coefs[0] = y0;
      m_table[ii].coefs[1] = -3*y0+4*y1-y2;
      m_table[ii].coefs[2] = 2*y0+-4*y1+2*y2;
    }
  }

  UniformQuadraticPrecomputedInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    MetaTable<TIN,TOUT,3,HORNER>(func_container, filename, "UniformQuadraticPrecomputedInterpolationTable") {}
};

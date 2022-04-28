/*
  Linear Interpolation LUT with uniform sampling (precomputed coefficients)

  Usage example:
    UniformLinearPrecomputedInterpolationTable look(&function,0,10,0.0001);
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
class UniformLinearPrecomputedInterpolationTable final : public MetaTable<TIN,TOUT,2,HORNER>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_UNIFORM_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,2,HORNER);
  FUNC_REGISTER_LUT(UniformLinearPrecomputedInterpolationTable);

public:
  UniformLinearPrecomputedInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, UniformLookupTableParameters<TIN> par) :
    MetaTable<TIN,TOUT,2,HORNER>(func_container, par)
  {
    /* Base class default variables */
    m_name = "UniformLinearPrecomputedInterpolationTable";
    m_order = 2;
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = (unsigned) sizeof(m_table[0]) * (m_numTableEntries);

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,2>[m_numTableEntries]);
    for (int ii=0; ii < m_numIntervals; ++ii) {
      TIN x = m_minArg + ii*m_stepSize;
      m_grid[ii] = x;
      m_table[ii].coefs[0] = m_func(x);
      x = m_minArg + (ii+1)*(m_stepSize);
      m_table[ii].coefs[1] = m_func(x) - m_table[ii].coefs[0];
    }
  }

  /* build this table from a file. Everything other than m_table is built by MetaTable */
  UniformLinearPrecomputedInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    MetaTable<TIN,TOUT,2,HORNER>(func_container, filename, "UniformLinearPrecomputedInterpolationTable") {}
  // operator() comes from MetaTable
};

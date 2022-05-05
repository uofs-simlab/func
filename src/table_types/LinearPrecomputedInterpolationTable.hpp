/*
  Linear Interpolation LUT with uniform sampling (precomputed coefficients)

  Usage example:
    LinearPrecomputedInterpolationTable look(&function,0,10,0.0001);
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
class LinearPrecomputedInterpolationTable final : public MetaTable<TIN,TOUT,2,HORNER,GT>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,2,HORNER,GT);
  FUNC_REGISTER_LUT(LinearPrecomputedInterpolationTable);

public:
  LinearPrecomputedInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, LookupTableParameters<TIN> par) :
    MetaTable<TIN,TOUT,2,HORNER,GT>(func_container, par)
  {
    /* Base class default variables */
    m_name = grid_type_to_string<GT>() + "LinearPrecomputedInterpolationTable";
    m_order = 2;
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = (unsigned) sizeof(m_table[0]) * (m_numTableEntries);

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,2>[m_numTableEntries]);
    for (int ii=0; ii < m_numIntervals; ++ii) {
      TIN x;
      // (possibly) transform the uniform grid into a nonuniform grid
      if (GT == UNIFORM)
        x = m_minArg + ii*m_stepSize;
      else
        x = m_transferFunction.g(m_minArg + ii*m_stepSize);

      m_grid[ii] = x;
      m_table[ii].coefs[0] = m_func(x);
      x = m_minArg + (ii+1)*(m_stepSize);
      m_table[ii].coefs[1] = m_func(x) - m_table[ii].coefs[0];
    }
  }

  /* build this table from a file. Everything other than m_table is built by MetaTable */
  LinearPrecomputedInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    MetaTable<TIN,TOUT,2,HORNER,GT>(func_container, filename,
        grid_type_to_string<GT>() + "LinearPrecomputedInterpolationTable") {}
  // operator() comes from MetaTable
};

// define friendlier names
template <typename TIN, typename TOUT=TIN>
using UniformLinearPrecomputedInterpolationTable = LinearPrecomputedInterpolationTable<TIN,TOUT,UNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformLinearPrecomputedInterpolationTable = LinearPrecomputedInterpolationTable<TIN,TOUT,NONUNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformPseudoLinearPrecomputedInterpolationTable = LinearPrecomputedInterpolationTable<TIN,TOUT,NONUNIFORM_PSEUDO>;

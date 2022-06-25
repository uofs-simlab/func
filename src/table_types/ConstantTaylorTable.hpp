/*
  Constant Taylor LUT with uniform sampling

  Usage example:
    ConstantTaylorTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "MetaTable.hpp"

template <typename TIN, typename TOUT=TIN, GridTypes GT=UNIFORM>
class ConstantTaylorTable final : public MetaTable<TIN,TOUT,1,GT>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,1,GT);

public:
  ConstantTaylorTable(FunctionContainer<TIN,TOUT> *func_container, LookupTableParameters<TIN> par) :
    MetaTable<TIN,TOUT,1,GT>(func_container, par)
  {
    /* Base class default variables */
    m_name = grid_type_to_string<GT>() + "ConstantTaylorTable";
    m_order = 1;
    m_numTableEntries = m_numIntervals;
    m_dataSize = (unsigned) sizeof(m_table[0]) * (m_numTableEntries);

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,1>[m_numTableEntries]);
    for (unsigned int ii=0; ii<m_numTableEntries; ++ii) {
      TIN x;
      // (possibly) transform the uniform grid into a nonuniform grid
      if (GT == UNIFORM)
        x = m_minArg + ii*m_stepSize;
      else
        x = m_transferFunction.g(m_minArg + ii*m_stepSize);

      m_grid[ii] = x;
      m_table[ii].coefs[0] = m_func(x);
    }
  }

  /* build this table from a file */
  ConstantTaylorTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    MetaTable<TIN,TOUT,1,GT>(func_container, filename,
        grid_type_to_string<GT>() + "ConstantTaylorTable") {}

  /* Constant interpolation from table point immediately below x */
  TOUT operator()(TIN x) override
  {
    return m_table[(unsigned)((x-m_minArg)/m_stepSize+0.5)].coefs[0];
  }
};

// define friendlier names
template <typename TIN, typename TOUT=TIN>
using UniformConstantTaylorTable = ConstantTaylorTable<TIN,TOUT,UNIFORM>;

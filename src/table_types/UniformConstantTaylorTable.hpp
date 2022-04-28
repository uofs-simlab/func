/*
  Constant Taylor LUT with uniform sampling

  Usage example:
    UniformConstantTaylorTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "MetaTable.hpp"

template <typename TIN, typename TOUT = TIN>
class UniformConstantTaylorTable final : public MetaTable<TIN,TOUT,1,TAYLOR>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_UNIFORM_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,1,TAYLOR);

  FUNC_REGISTER_LUT(UniformConstantTaylorTable);

  //TOUT get_table_entry(unsigned int i, unsigned int j) override { return m_table[i].coefs[0]; (void) j; }

public:
  UniformConstantTaylorTable(FunctionContainer<TIN,TOUT> *func_container, UniformLookupTableParameters<TIN> par) :
    MetaTable<TIN,TOUT,1,TAYLOR>(func_container, par)
  {
    /* Base class default variables */
    m_name = "UniformConstantTaylorTable";
    m_order = 1;
    m_numTableEntries = m_numIntervals;
    m_dataSize = (unsigned) sizeof(m_table[0]) * (m_numTableEntries);

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,1>[m_numTableEntries]);
    for (int ii=0; ii<m_numTableEntries; ++ii) {
      TIN x = m_minArg + ii*m_stepSize;
      m_grid[ii] = x;
      m_table[ii].coefs[0] = m_func(x);
    }
  }

  /* build this table from a file */
  UniformConstantTaylorTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    MetaTable<TIN,TOUT,1,TAYLOR>(func_container, filename, "UniformConstantTaylorTable") {}

  /* Constant interpolation from table point immediately below x */
  TOUT operator()(TIN x) override
  {
    return m_table[(unsigned)((x-m_minArg)/m_stepSize+0.5)].coefs[0];
  }
};

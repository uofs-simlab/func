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
#include "UniformLookupTable.hpp"

template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE>
class UniformLinearPrecomputedInterpolationTable final : public UniformLookupTable<IN_TYPE,OUT_TYPE>
{
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE);
  FUNC_REGISTER_LUT(UniformLinearPrecomputedInterpolationTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,2>[]> m_table;
  OUT_TYPE get_table_entry(unsigned int i, unsigned int j) override { return m_table[i].coefs[j]; }
  unsigned int get_num_coefs() override { return m_table[0].num_coefs; }

public:
  UniformLinearPrecomputedInterpolationTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, UniformLookupTableParameters<IN_TYPE> par) :
    UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, par)
  {
    /* Base class default variables */
    m_name = FUNC_STR(UniformLinearPrecomputedInterpolationTable);
    m_order = 2;
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;

    /* Allocate and set table */
    m_table.reset(new polynomial<OUT_TYPE,2>[m_numTableEntries]);
    for (int ii=0; ii < m_numIntervals; ++ii) {
      IN_TYPE x = m_minArg + ii*m_stepSize;
      m_grid[ii] = x;
      m_table[ii].coefs[0] = mp_func(x);
      x = m_minArg + (ii+1)*(m_stepSize);
      m_table[ii].coefs[1] = mp_func(x) - m_table[ii].coefs[0];
    }
  }

  /* build this table from a file. Everything other than m_table is built by UniformLookupTable */
  UniformLinearPrecomputedInterpolationTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, std::string filename) :
    UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, filename)
  {
    std::ifstream file_reader(filename);
    using nlohmann::json;
    json jsonStats;
    file_reader >> jsonStats;

    // double check the names match
    std::string temp_name = jsonStats["name"].get<std::string>();
    if(temp_name != "UniformLinearPrecomputedInterpolationTable")
      throw std::invalid_argument("Error while reading " + filename + ": "
          "Cannot build a " + temp_name + " from a UniformLinearPrecomputedInterpolationTable");

    m_table.reset(new polynomial<OUT_TYPE,2>[m_numTableEntries]);
    for(unsigned int i=0; i<m_numTableEntries; i++)
      for(unsigned int j=0; j<m_table[i].num_coefs; j++)
        m_table[i].coefs[j] = jsonStats["table"][std::to_string(i)]["coefs"][std::to_string(j)].get<OUT_TYPE>();
  }


  OUT_TYPE operator()(IN_TYPE x) override
  {
    // nondimensionalized x position, scaled by step size
    OUT_TYPE dx = m_stepSize_inv*(x-m_minArg);
    // index of previous table entry
    unsigned x0  = (unsigned) dx;
    dx -= x0;
    // linear interpolation
    return m_table[x0].coefs[0]+dx*m_table[x0].coefs[1];
  }
};

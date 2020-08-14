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
#include "UniformLookupTable.hpp"

template <typename IN_TYPE, typename OUT_TYPE>
class UniformConstantTaylorTable final : public UniformLookupTable<IN_TYPE,OUT_TYPE>
{
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE);

  REGISTER_LUT(UniformConstantTaylorTable);

  __attribute__((aligned)) std::unique_ptr<OUT_TYPE[]> m_table;
  OUT_TYPE get_table_entry(unsigned int i, unsigned int j) override { return m_table[i]; (void) j; }
  unsigned int get_num_coefs() override { return 1; }

public:
  UniformConstantTaylorTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, UniformLookupTableParameters<IN_TYPE> par) :
    UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, par)
  {
    /* Base class default variables */
    m_name = STR(UniformConstantTaylorTable);
    m_order = 1;
    m_numTableEntries = m_numIntervals;
    m_dataSize = (unsigned) sizeof(OUT_TYPE) * m_numTableEntries;

    /* Allocate and set table */
    m_table.reset(new OUT_TYPE[m_numTableEntries]);
    for (int ii=0; ii<m_numTableEntries; ++ii) {
      IN_TYPE x = m_minArg + ii*m_stepSize;
      m_grid[ii] = x;
      m_table[ii] = mp_func(x);
    }
  }

  /* build this table from a file. Everything other than m_table is built by UniformLookupTable */
  UniformConstantTaylorTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, std::string filename) :
    UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, filename)
  {
    std::ifstream file_reader(filename);
    using nlohmann::json;
    json jsonStats;
    file_reader >> jsonStats;

    // double check the names match
    std::string temp_name = jsonStats["name"].get<std::string>();
    if(temp_name != "UniformConstantTaylorTable")
      throw std::invalid_argument("Error while reading " + filename + ": "
          "Cannot build a " + temp_name + " from a UniformConstantTaylorTable");

    m_table.reset(new polynomial<OUT_TYPE,1>[m_numTableEntries]);
    for(unsigned int i=0; i<m_numTableEntries; i++)
      for(unsigned int j=0; j<m_table[i].num_coefs; j++)
        m_table[i].coefs[j] = jsonStats["table"][std::to_string(i)]["coefs"][std::to_string(j)].get<OUT_TYPE>();
  }

  /* Constant interpolation from table point immediately below x */
  OUT_TYPE operator()(IN_TYPE x) override
  {
    return m_table[(unsigned)((x-m_minArg)/m_stepSize+0.5)];
  }
};

REGISTER_DOUBLE_AND_FLOAT_LUT_IMPLS(UniformConstantTaylorTable);

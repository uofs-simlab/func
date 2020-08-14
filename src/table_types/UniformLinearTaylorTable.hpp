/*
  Linear Taylor LUT with uniform sampling

  Usage example:
    UniformLinearTaylorTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "UniformLookupTable.hpp"

template <typename IN_TYPE, typename OUT_TYPE>
class UniformLinearTaylorTable final : public UniformLookupTable<IN_TYPE,OUT_TYPE>
{
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE);
  REGISTER_LUT(UniformLinearTaylorTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,2>[]> m_table;
  std::function<fvar<OUT_TYPE,1>(fvar<IN_TYPE,1>)> mp_boost_func;
  OUT_TYPE get_table_entry(unsigned int i, unsigned int j) override { return m_table[i].coefs[j]; }
  unsigned int get_num_coefs() override { return m_table[0].num_coefs; }

public:
  UniformLinearTaylorTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, UniformLookupTableParameters<IN_TYPE> par) :
    UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, par)
  {
    using boost::math::differentiation::make_fvar;

    /* Base class variables */
    m_name = STR(UniformLinearTaylorTable);
    m_order = 2;
    m_numTableEntries = m_numIntervals;
    m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;

    __IS_NULL(func_container->autodiff1_func);
    mp_boost_func = func_container->autodiff1_func;

    /* Allocate and set table */
    m_table.reset(new polynomial<OUT_TYPE,2>[m_numTableEntries]);
    for (int ii=0; ii<m_numIntervals; ++ii) {
      const IN_TYPE x = m_minArg + ii*m_stepSize;
      m_grid[ii]     = x;
      // get every derivative up to the first
      auto const derivs = (mp_boost_func)(make_fvar<IN_TYPE,1>(x));
      m_table[ii].coefs[0] = derivs.derivative(0);
      m_table[ii].coefs[1] = derivs.derivative(1);
    }
  }

  /* build this table from a file. Everything other than m_table is built by UniformLookupTable */
  UniformLinearTaylorTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, std::string filename) :
    UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, filename)
  {
    std::ifstream file_reader(filename);
    using nlohmann::json;
    json jsonStats;
    file_reader >> jsonStats;

    // double check the names match
    std::string temp_name = jsonStats["name"].get<std::string>();
    if(temp_name != "UniformLinearTaylorTable")
      throw std::invalid_argument("Error while reading " + filename + ": "
          "Cannot build a " + temp_name + " from a UniformLinearTaylorTable");

    m_table.reset(new polynomial<OUT_TYPE,2>[m_numTableEntries]);
    for(unsigned int i=0; i<m_numTableEntries; i++)
      for(unsigned int j=0; j<m_table[i].num_coefs; j++)
        m_table[i].coefs[j] = jsonStats["table"][std::to_string(i)]["coefs"][std::to_string(j)].get<OUT_TYPE>();
  }

  OUT_TYPE operator()(IN_TYPE x) override
  {
    // nondimensionalized x position
    OUT_TYPE dx  = (x-m_minArg);
    OUT_TYPE x1r = dx/m_stepSize+0.5;
    // index of previous table entry
    unsigned x1  = (unsigned) x1r;
    dx -= x1*m_stepSize;
    // linear Taylor series from grid point below x
    return m_table[x1].coefs[0]+m_table[x1].coefs[1]*dx;
  }

  std::function<fvar<OUT_TYPE,1>(fvar<IN_TYPE,1>)> boost_function(){ return mp_boost_func; }
};

REGISTER_DOUBLE_AND_FLOAT_LUT_IMPLS(UniformLinearTaylorTable);

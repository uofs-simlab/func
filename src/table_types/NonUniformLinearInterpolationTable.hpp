/*
  TODO add an enum template paramter to every non uniform table which changes operator()
  Potential implemetation:
    enum class hashType {standard, pseudo};
    template <typename T, TABLE_HASH = hashType::standard>
    NonUniformLinearInterpolationTable ...
    {
      ...
      T operator()(T x){ return operator()(x,TABLE_HASH); }
      T operator()(T x, hashType::standard){ ... } // this table's hash
      T operator()(T x, hashType::pseudo){ ... } // the pseudo table's hash
    }

  Then alias the following
    template<typename T>
    using NonUniformPseudoLinearInterpolationTable = NonUniformLinearInterpolationTable<T,hashType::pseudo>

  Linear Interpolation LUT with nonuniform sampling

  Usage example:
    NonUniformLinearInterpolationTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "NonUniformLookupTable.hpp"

template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE, class TRANSFER_FUNC_TYPE = TransferFunctionSinh<IN_TYPE>>
class NonUniformLinearInterpolationTable final : public NonUniformLookupTable<IN_TYPE,OUT_TYPE,TRANSFER_FUNC_TYPE>
{
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE);
  INHERIT_NONUNIFORM_LUT(IN_TYPE,OUT_TYPE);

  FUNC_REGISTER_LUT(NonUniformLinearInterpolationTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,1>[]> m_table;
  OUT_TYPE get_table_entry(unsigned int i, unsigned int j) override { return m_table[i].coefs[j]; }
  unsigned int get_num_coefs() override { return m_table[0].num_coefs; }
public:
  NonUniformLinearInterpolationTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container,
      UniformLookupTableParameters<IN_TYPE> par) :
    NonUniformLookupTable<IN_TYPE,OUT_TYPE,TRANSFER_FUNC_TYPE>(func_container, par)
  {
    /* Base class variables */
    m_name = FUNC_STR(NonUniformLinearInterpolationTable);
    m_order = 2;
    m_numTableEntries = m_numIntervals;
    m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;
    
    /* Allocate and set table */
    m_table.reset(new polynomial<OUT_TYPE,1>[m_numTableEntries]);
    for (int ii=0; ii<m_numIntervals; ++ii) {
      // transform the previously used uniform grid to a nonuniform grid
      const IN_TYPE x = m_transferFunction.g(m_minArg + ii*m_stepSize);
      m_grid[ii]  = x;
      m_table[ii].coefs[0] = m_func(x);
    }
  }

  /* build this table from a file. Everything other than m_table is built by UniformLookupTable */
  NonUniformLinearInterpolationTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, std::string filename) :
    NonUniformLookupTable<IN_TYPE,OUT_TYPE,TRANSFER_FUNC_TYPE>(func_container, filename)
  {
    std::ifstream file_reader(filename);
    using nlohmann::json;
    json jsonStats;
    file_reader >> jsonStats;

    // double check the names match
    std::string temp_name = jsonStats["name"].get<std::string>();
    if(temp_name != "NonUniformLinearInterpolationTable")
      throw std::invalid_argument("Error while reading " + filename + ": "
          "Cannot build a " + temp_name + " from a NonUniformLinearInterpolationTable");

    m_table.reset(new polynomial<OUT_TYPE,1>[m_numTableEntries]);
    for(unsigned int i=0; i<m_numTableEntries; i++)
      for(unsigned int j=0; j<m_table[i].num_coefs; j++)
        m_table[i].coefs[j] = jsonStats["table"][std::to_string(i)]["coefs"][std::to_string(j)].get<OUT_TYPE>();
  }

  OUT_TYPE operator()(IN_TYPE x) override
  {
    // set x0 = floor((g_inv(x)-m_minArg)/m_stepSize)
    // where each of the above member vars are encoded into g_inv
    unsigned x0 = m_transferFunction.g_inv(x);
    IN_TYPE h   = m_grid[x0+1] - m_grid[x0];
    OUT_TYPE dx = (x - m_grid[x0])/h;

    // value of table entries around x position
    OUT_TYPE y1  = m_table[x0].coefs[0];
    OUT_TYPE y2  = m_table[x0+1].coefs[0];
    // linear interpolation
    return y1+dx*(y2-y1);
  }
};

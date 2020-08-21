/*
  Cubic polynomial Interpolation LUT with nonuniform sampling

  Usage example:
    NonUniformPseudoCubicPrecomputedInterpolationTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "NonUniformLookupTable.hpp"

template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE, class TRANSFER_FUNC_TYPE = TransferFunctionSinh<IN_TYPE>>
class NonUniformPseudoCubicPrecomputedInterpolationTable final : public NonUniformLookupTable<IN_TYPE,OUT_TYPE,TRANSFER_FUNC_TYPE>
{
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE);
  INHERIT_NONUNIFORM_LUT(IN_TYPE,OUT_TYPE);

  FUNC_REGISTER_LUT(NonUniformPseudoCubicPrecomputedInterpolationTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,4>[]> m_table;
  OUT_TYPE get_table_entry(unsigned int i, unsigned int j) override { return m_table[i].coefs[j]; }
  unsigned int get_num_coefs() override { return m_table[0].num_coefs; }

public:
  NonUniformPseudoCubicPrecomputedInterpolationTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container,
      UniformLookupTableParameters<IN_TYPE> par) :
    NonUniformLookupTable<IN_TYPE,OUT_TYPE,TRANSFER_FUNC_TYPE>(func_container, par)
  {
    /* Base class variables */
    m_name = FUNC_STR(NonUniformPseudoCubicPrecomputedInterpolationTable);
    m_order = 4;
    m_numTableEntries = m_numIntervals + 1;
    m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;
    
    /* Allocate and set table */
    m_table.reset(new polynomial<OUT_TYPE,4>[m_numTableEntries]);
    for (int ii=0; ii<m_numIntervals; ++ii) {
      const IN_TYPE x = m_transferFunction.g(m_minArg + ii*m_stepSize);
      // the local stepsize of our nonuniform grid
      const IN_TYPE h = m_transferFunction.g(m_minArg + (ii+1)*m_stepSize) - x;
      // grid points
      m_grid[ii] = x;
      // polynomial coefficients
      const OUT_TYPE y0 = mp_func(x);
      const OUT_TYPE y1 = mp_func(x+h/3);
      const OUT_TYPE y2 = mp_func(x+2*h/3);
      const OUT_TYPE y3 = mp_func(x+h);
      m_table[ii].coefs[0] = y0;
      m_table[ii].coefs[1] = -11*y0/2+9*y1-9*y2/2+y3;
      m_table[ii].coefs[2] = 9*y0-45*y1/2+18*y2-9*y3/2;
      m_table[ii].coefs[3] = -9*y0/2+27*y1/2-27*y2/2+9*y3/2;
    }
  }

  /* build this table from a file. Everything other than m_table is built by UniformLookupTable */
  NonUniformPseudoCubicPrecomputedInterpolationTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, std::string filename) :
    NonUniformLookupTable<IN_TYPE,OUT_TYPE,TRANSFER_FUNC_TYPE>(func_container, filename)
  {
    std::ifstream file_reader(filename);
    using nlohmann::json;
    json jsonStats;
    file_reader >> jsonStats;

    // double check the names match
    std::string temp_name = jsonStats["name"].get<std::string>();
    if(temp_name != "NonUniformPseudoCubicPrecomputedInterpolationTable")
      throw std::invalid_argument("Error while reading " + filename + ": "
          "Cannot build a " + temp_name + " from a NonUniformPseudoCubicPrecomputedInterpolationTable");

    m_table.reset(new polynomial<OUT_TYPE,4>[m_numTableEntries]);
    for(unsigned int i=0; i<m_numTableEntries; i++)
      for(unsigned int j=0; j<m_table[i].num_coefs; j++)
        m_table[i].coefs[j] = jsonStats["table"][std::to_string(i)]["coefs"][std::to_string(j)].get<OUT_TYPE>();
  }


  OUT_TYPE operator()(IN_TYPE x) override
  {
    // find the subinterval x lives in
    OUT_TYPE dx = m_transferFunction.g_inv(x);
    unsigned x0 = (unsigned) dx;
    dx -= x0;

    // cubic interpolation
    return m_table[x0].coefs[0]+dx*(m_table[x0].coefs[1]+dx*(m_table[x0].coefs[2]+dx*m_table[x0].coefs[3]));
  }
};

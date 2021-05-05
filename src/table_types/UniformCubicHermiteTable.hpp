/*
  Cubic Interpolation LUT with uniform sampling (precomputed coefficients)

  Usage example:
    UniformCubicHermiteTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "UniformLookupTable.hpp"
#include "config.hpp"

#ifndef FUNC_USE_BOOST_AUTODIFF
#error "UniformCubicHermiteTable needs boost version >= 1.71"
#endif

template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE>
class UniformCubicHermiteTable final : public UniformLookupTable<IN_TYPE,OUT_TYPE>
{
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE);
  FUNC_REGISTER_LUT(UniformCubicHermiteTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,4>[]> m_table;
  std::function<adVar<OUT_TYPE,1>(adVar<IN_TYPE,1>)> mp_boost_func;
  OUT_TYPE get_table_entry(unsigned int i, unsigned int j) override { return m_table[i].coefs[j]; }
  unsigned int get_num_coefs() override { return m_table[0].num_coefs; }

public:
  UniformCubicHermiteTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, UniformLookupTableParameters<IN_TYPE> par) :
      UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, par)
  {
    using boost::math::differentiation::make_fvar;
    /* Base class default variables */
    m_name = FUNC_STR(UniformCubicHermiteTable);
    m_order = 4;
    m_numTableEntries = m_numIntervals+1;
    m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;

    __IS_NULL(func_container->autodiff1_func);
    mp_boost_func = func_container->autodiff1_func;

    /* Allocate and set table */
    m_table.reset(new polynomial<OUT_TYPE,4>[m_numTableEntries]);
    for (int ii=0; ii<m_numIntervals; ++ii) {
      const IN_TYPE x = m_minArg + ii*m_stepSize;
      m_grid[ii] = x;

      const auto derivs0 = (mp_boost_func)(make_fvar<IN_TYPE,1>(x));
      const OUT_TYPE y0    = derivs0.derivative(0);
      const OUT_TYPE m0    = derivs0.derivative(1);
      const auto derivs1 = (mp_boost_func)(make_fvar<IN_TYPE,1>(x+m_stepSize));
      const OUT_TYPE y1    = derivs1.derivative(0);
      const OUT_TYPE m1    = derivs1.derivative(1);

      m_table[ii].coefs[0] = y0;
      m_table[ii].coefs[1] = m_stepSize*m0;
      m_table[ii].coefs[2] = -3*y0+3*y1-(2*m0+m1)*m_stepSize;
      m_table[ii].coefs[3] = 2*y0-2*y1+(m0+m1)*m_stepSize;
    }
  }

  /* build this table from a file. Everything other than m_table is built by UniformLookupTable */
  UniformCubicHermiteTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, std::string filename) :
    UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, filename)
  {
    std::ifstream file_reader(filename);
    using nlohmann::json;
    json jsonStats;
    file_reader >> jsonStats;

    // double check the names match
    std::string temp_name = jsonStats["name"].get<std::string>();
    if(temp_name != "UniformCubicHermiteTable")
      throw std::invalid_argument("Error while reading " + filename + ": "
          "Cannot build a " + temp_name + " from a UniformCubicHermiteTable");

    m_table.reset(new polynomial<OUT_TYPE,4>[m_numTableEntries]);
    for(unsigned int i=0; i<m_numTableEntries; i++)
      for(unsigned int j=0; j<m_table[i].num_coefs; j++)
        m_table[i].coefs[j] = jsonStats["table"][std::to_string(i)]["coefs"][std::to_string(j)].get<OUT_TYPE>();
  }

  OUT_TYPE operator()(IN_TYPE x) override
  {
    // nondimensionalized x position, scaled by step size
    OUT_TYPE  dx = m_stepSize_inv*(x-m_minArg);
    // index of previous table entry
    unsigned x0  = (unsigned) dx;
    dx -= x0;
    // linear interpolation
    return m_table[x0].coefs[0]+dx*(m_table[x0].coefs[1]+dx*(m_table[x0].coefs[2]+dx*m_table[x0].coefs[3]));
  }

  std::function<adVar<OUT_TYPE,1>(adVar<OUT_TYPE,1>)> boost_function(){ return mp_boost_func; }
};

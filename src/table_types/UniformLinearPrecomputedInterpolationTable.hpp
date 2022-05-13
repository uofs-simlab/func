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
#include "UniformLookupTable.hpp"

template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE>
class UniformLinearPrecomputedInterpolationTable final : public UniformLookupTable<IN_TYPE,OUT_TYPE>
{
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE);
  REGISTER_LUT(UniformLinearPrecomputedInterpolationTable);

  __attribute__((aligned)) std::unique_ptr<polynomial<OUT_TYPE,2>[]> m_table;
public:
  UniformLinearPrecomputedInterpolationTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, UniformLookupTableParameters<IN_TYPE> par) :
    UniformLookupTable<IN_TYPE,OUT_TYPE>(func_container, par)
  {
    /* Base class default variables */
    m_name = STR(UniformLinearPrecomputedInterpolationTable);
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

REGISTER_DOUBLE_AND_FLOAT_LUT_IMPLS(UniformLinearPrecomputedInterpolationTable);

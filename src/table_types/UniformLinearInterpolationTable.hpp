/*
  Linear Interpolation LUT with uniform sampling

  Usage example:
    UniformLinearInterpolationTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "MetaTable.hpp"

template <typename TIN, typename TOUT = TIN>
class UniformLinearInterpolationTable final : public MetaTable<TIN,TOUT,1,HORNER>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_UNIFORM_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,1,HORNER);
  FUNC_REGISTER_LUT(UniformLinearInterpolationTable);
 
public:
  //#pragma omp declare simd
  UniformLinearInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, UniformLookupTableParameters<TIN> par) :
    MetaTable<TIN,TOUT,1,HORNER>(func_container, par)
  {
    /* Base class variables */
    m_name = "UniformLinearInterpolationTable";
    m_order = 1;
    m_numTableEntries = m_numIntervals;
    m_dataSize = (unsigned) sizeof(m_table[0]) * (m_numTableEntries);

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,1>[m_numTableEntries]);

    /* I think it would be nice to have simd table generation so something like the LUT generator can use actual
     * parallelism. The alignment of m_table might also give us a big speedup from simd */
    //#pragma omp simd aligned(m_table:sizeof(TOUT)) // needs the constructor to be declared simd
    // assuming each iteration will take about the same amount of time
    //#pragma omp parallel for schedule(static)
    for (int ii=0; ii<m_numIntervals; ++ii) {
      const TIN x = m_minArg + ii*m_stepSize;
      m_grid[ii]  = x;
      m_table[ii].coefs[0] = m_func(x);
    }
  }

  /* build this table from a file. Everything other than m_table is built by MetaTable */
  UniformLinearInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    MetaTable<TIN,TOUT,1,HORNER>(func_container, filename, "UniformLinearInterpolationTable") {}

  // operator() is slightly different from MetaTable's HORNER method
  TOUT operator()(TIN x) override
  {
    // nondimensionalized x position, scaled by step size
    TOUT dx = (x-m_minArg)/m_stepSize;
    // index of previous table entry
    unsigned x0 = (unsigned) dx;
    // value of table entries around x position
    dx -= x0;
    TOUT y1  = m_table[x0].coefs[0];
    TOUT y2  = m_table[x0+1].coefs[0];
    // linear interpolation
    return y1+dx*(y2-y1);
  }
};

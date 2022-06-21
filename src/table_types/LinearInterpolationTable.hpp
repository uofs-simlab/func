/*
  Linear Interpolation LUT. Coefficients are computed at lookup time.
  Approx 50% less memory usage compared to LinearPrecomputedInterpolationTable
  but the hash is slower.

  Usage example:
    LinearInterpolationTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "MetaTable.hpp"

template <typename TIN, typename TOUT=TIN, GridTypes GT=UNIFORM>
class LinearInterpolationTable final : public MetaTable<TIN,TOUT,1,HORNER,GT>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,1,HORNER,GT);

public:
  //#pragma omp declare simd
  LinearInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, LookupTableParameters<TIN> par) :
    MetaTable<TIN,TOUT,1,HORNER,GT>(func_container, par)
  {
    /* Base class variables */
    m_name  = grid_type_to_string<GT>() + "LinearInterpolationTable";
    m_order = 1;
    m_numTableEntries = m_numIntervals;
    m_dataSize = (unsigned) sizeof(m_table[0]) * (m_numTableEntries);

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,1>[m_numTableEntries]);
    for (unsigned int ii=0; ii<m_numIntervals; ++ii) {
      TIN x;
      // (possibly) transform the uniform grid into a nonuniform grid
      if (GT == UNIFORM)
        x = m_minArg + ii*m_stepSize;
      else
        x = m_transferFunction.g(m_minArg + ii*m_stepSize);

      m_grid[ii]  = x;
      m_table[ii].coefs[0] = m_func(x);
    }
  }

  /* build this table from a file. Everything other than m_table is built by MetaTable */
  LinearInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    MetaTable<TIN,TOUT,1,HORNER,GT>(func_container, filename,
        grid_type_to_string<GT>() + "LinearInterpolationTable") {}

  // operator() is slightly different from MetaTable's HORNER method
  TOUT operator()(TIN x) override
  {
    //enum GridTypes {UNIFORM, NONUNIFORM, NONUNIFORM_PSEUDO};
    TOUT dx;
    unsigned int x0;
    switch(GT){
    case UNIFORM:
      {
      // nondimensionalized x position, scaled by step size
      dx = (x-m_minArg)/m_stepSize;
      // index of previous table entry
      x0 = (unsigned) dx;
      // value of table entries around x position
      dx -= x0;
      break;
      }
    case NONUNIFORM:
      {
      // set x0 = floor((g_inv(x)-m_minArg)/m_stepSize)
      // where each of the above member vars are encoded into g_inv
      x0 = m_transferFunction.g_inv(x);
      TIN h   = m_grid[x0+1] - m_grid[x0];
      dx = (x - m_grid[x0])/h;
      break;
      }
    case NONUNIFORM_PSEUDO:
      {
      // set x0 = floor((g_inv(x)-m_minArg)/m_stepSize)
      // and dx = fractional((g_inv(x)-m_minArg)/m_stepSize)
      // where each of the above member vars are encoded into g_inv
      // The source of the pseudolinearity is from the way we compute dx
      dx = m_transferFunction.g_inv(x);
      x0 = (unsigned) dx;
      dx -= x0;
      break;
      }
    }

    TOUT y1  = m_table[x0].coefs[0];
    TOUT y2  = m_table[x0+1].coefs[0];
    // linear interpolation
    return y1+dx*(y2-y1);
  }
};

template <typename TIN, typename TOUT=TIN>
using UniformLinearInterpolationTable = LinearInterpolationTable<TIN,TOUT,UNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformLinearInterpolationTable = LinearInterpolationTable<TIN,TOUT,NONUNIFORM>;
template <typename TIN, typename TOUT=TIN>
using NonUniformPseudoLinearInterpolationTable = LinearInterpolationTable<TIN,TOUT,NONUNIFORM_PSEUDO>;

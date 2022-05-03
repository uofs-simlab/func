/*
  Linear Interpolation LUT with uniform sampling

  Usage example:
    LinearInterpolationTable look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
*/
#pragma once
#include "MetaTable.hpp"

enum IsUniform {UNIFORM, NONUNIFORM, NONUNIFORM_PSEUDO};

template <typename TIN, typename TOUT, IsUniform IU>
class LinearInterpolationTable final : public MetaTable<TIN,TOUT,1,HORNER>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_UNIFORM_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,1,HORNER);
  FUNC_REGISTER_LUT(LinearInterpolationTable);
 
public:
  //#pragma omp declare simd
  LinearInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, UniformLookupTableParameters<TIN> par) :
    MetaTable<TIN,TOUT,1,HORNER>(func_container, par)
  {
    /* Base class variables */
    m_name = "LinearInterpolationTable";
    //m_name = FUNC_STR(NonUniformLinearInterpolationTable);
    //m_name = FUNC_STR(NonUniformPseudoLinearInterpolationTable);
    m_order = 1;
    m_numTableEntries = m_numIntervals;
    m_dataSize = (unsigned) sizeof(m_table[0]) * (m_numTableEntries);

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,1>[m_numTableEntries]);
    for (int ii=0; ii<m_numIntervals; ++ii) {
      TIN x;
      // transform the uniform grid to a nonuniform grid
      if (IU = NONUNIFORM)
        x = m_transferFunction.g(m_minArg + ii*m_stepSize);
      else
        x = m_minArg + ii*m_stepSize;

      m_grid[ii]  = x;
      m_table[ii].coefs[0] = m_func(x);
    }
  }

  /* build this table from a file. Everything other than m_table is built by MetaTable */
  LinearInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    MetaTable<TIN,TOUT,1,HORNER>(func_container, filename, "LinearInterpolationTable") {}

  // operator() is slightly different from MetaTable's HORNER method
  TOUT operator()(TIN x) override
  {
    switch(IU){
    case UNIFORM:
      // nondimensionalized x position, scaled by step size
      TOUT dx = (x-m_minArg)/m_stepSize;
      // index of previous table entry
      unsigned x0 = (unsigned) dx;
      // value of table entries around x position
      dx -= x0;
      break;
    case NONUNIFORM_PSEUDO:
      // set x0 = floor((g_inv(x)-m_minArg)/m_stepSize)
      // and dx = fractional((g_inv(x)-m_minArg)/m_stepSize)
      // where each of the above member vars are encoded into g_inv
      // The source of the pseudolinearity is from the way we compute dx
      OUT_TYPE dx = m_transferFunction.g_inv(x);
      unsigned x0 = (unsigned) dx;
      dx -= x0;
      break;
    case NONUNIFORM:
      // set x0 = floor((g_inv(x)-m_minArg)/m_stepSize)
      // where each of the above member vars are encoded into g_inv
      unsigned x0 = m_transferFunction.g_inv(x);
      IN_TYPE h   = m_grid[x0+1] - m_grid[x0];
      OUT_TYPE dx = (x - m_grid[x0])/h;
      break;
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

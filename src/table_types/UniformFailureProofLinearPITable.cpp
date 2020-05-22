/* Implementation of a UniformPrecomputed Lookup table with linear interpolation */
#include "UniformFailureProofLinearPITable.hpp"

// RECORD_ARG and PRINT_ARGS are null operators if NDEBUG is specified at compile time
#ifdef NDEBUG
  #define RECORD_ARG(x)
  #define PRINT_ARGS(out)
#else
  #include <iostream>
  #define RECORD_ARG(x) m_args.push_back((x))
  // make sure we don't swallow the semicolon
   #define PRINT_ARGS(out)                                 \
    do {                                                  \
      if(!m_args.empty()){                                \
        out << "args outside table range:" << std::endl;  \
        out << m_args.front();                            \
        m_args.erase(m_args.begin());                     \
        for(auto x : m_args)                              \
          out << ", " << x;                               \
        out << std::endl;                                 \
      }                                                   \
    } while(0) 
#endif

#define IMPL_NAME UniformFailureProofLinearPITable
REGISTER_ULUT_IMPL(IMPL_NAME);

UniformFailureProofLinearPITable::UniformFailureProofLinearPITable(EvaluationFunctor<double,double> *func, UniformLookupTableParameters par) : UniformLookupTable(func, par)
{

  /* Base class default variables */
  m_name = STR(IMPL_NAME);
  m_order = 2;
  m_numTableEntries = 2*m_numIntervals+2;
  m_dataSize = (unsigned) sizeof(double) * m_numTableEntries;

  /* Allocate and set table */
  m_table.reset(new double[m_numTableEntries]);
  for (int ii=0;ii<m_numIntervals;++ii) {
    double x = m_minArg + ii*m_stepSize;
    m_grid[ii] = x;
    m_table[2*ii]   = (*mp_func)(x);
    x = m_minArg + (ii+1)*m_stepSize;
    m_table[2*ii+1] = (*mp_func)(x) - m_table[2*ii];
  }
}

double UniformFailureProofLinearPITable::operator()(double x)
{
  // check if x is in the range of the table
  if(x<m_minArg || x>m_maxArg){
    RECORD_ARG(x);
    return (*mp_func)(x);
  }

  // nondimensionalized x position, scaled by step size
  double   dx = m_stepSize_inv*(x-m_minArg);
  // index of previous table entry
  unsigned x0  = (unsigned) dx;
  dx -= x0;
  x0 *= 2;
  // linear interpolation
  return m_table[x0]+dx*m_table[x0+1];
}

UniformFailureProofLinearPITable::~UniformFailureProofLinearPITable()
{
  PRINT_ARGS(std::cerr);
}

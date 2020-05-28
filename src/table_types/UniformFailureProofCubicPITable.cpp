/* Implementation of a UniformPrecomputed Lookup table with linear interpolation */
#include "UniformFailureProofCubicPITable.hpp"

#ifdef NDEBUG
  #define RECORD_ARG(x)
  #define PRINT_ARGS(out)
#else
  #include <iostream>
  #include <vector>
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

#define IMPL_NAME UniformFailureProofCubicPITable
REGISTER_ULUT_IMPL(IMPL_NAME);

UniformFailureProofCubicPITable::UniformFailureProofCubicPITable(EvaluationFunctor<double,double> *func, UniformLookupTableParameters par) : UniformLookupTable(func, par)
{

  /* Base class default variables */
  m_name = STR(IMPL_NAME);
  m_order = 4;
  m_numTableEntries = 4*(m_numIntervals+1);
  m_dataSize = (unsigned) sizeof(m_table[0]) * m_numTableEntries;

  /* Allocate and set table */
  m_table.reset(new polynomial<4,32>[m_numTableEntries]);
  for (int ii=0;ii<m_numIntervals;++ii) {
    const double x = m_minArg + ii*m_stepSize;
    // grid points
    m_grid[ii] = x;
    // polynomial coefficients
    const double y0 = (*mp_func)(x);
    const double y1 = (*mp_func)(x+m_stepSize/3);
    const double y2 = (*mp_func)(x+2*m_stepSize/3);
    const double y3 = (*mp_func)(x+m_stepSize);
    m_table[ii].coefs[0] = y0;
    m_table[ii].coefs[1] = -11*y0/2+9*y1-9*y2/2+y3;
    m_table[ii].coefs[2] = 9*y0-45*y1/2+18*y2-9*y3/2;
    m_table[ii].coefs[3] = -9*y0/2+27*y1/2-27*y2/2+9*y3/2;
  }
}

double UniformFailureProofCubicPITable::operator()(double x)
{
  // check if x is in the range of the table
  if(x<m_minArg || x>m_maxArg){
    RECORD_ARG(x);
    return (*mp_func)(x);
  }

  // nondimensionalized x position, scaled by step size
  double   dx = m_stepSize_inv*(x-this->m_minArg);
  // index of previous table entry
  // unsigned x0  = (unsigned) floor(dx);
  unsigned x0  = (unsigned) dx;
  // value of table entries around x position
  dx -= x0;
  // cubic interpolation
  return m_table[x0].coefs[0]+dx*(m_table[x0].coefs[1]+dx*(m_table[x0].coefs[2]+dx*m_table[x0].coefs[3]));
}

UniformFailureProofCubicPITable::~UniformFailureProofCubicPITable()
{
  PRINT_ARGS(std::cerr);
}

/*
   N = number of coefficients used in underlying polynomial
   Defines several polynomial based arrays which we can use to
   make actual LUTs succinctly
*/
#pragma once
#include "UniformLookupTable.hpp"

#define INHERIT_META(TIN,TOUT,N,HT) \
  using MetaTable<TIN,TOUT,N,HT>::m_table;

/* Parallelization macro. Play around with this to see which OpenMP option for parallelizing a for loop is best
 * Might be nice to have simd table generation so something like the LUT generator can use actual
 * parallelism. The alignment of m_table might also give us a big speedup from simd */
//#pragma omp simd aligned(m_table:sizeof(TOUT)) // needs the constructor to be declared simd
// assuming each iteration will take about the same amount of time
//#pragma omp parallel for schedule(static)

enum HashTypes {HORNER, TAYLOR};

template <typename TIN, typename TOUT, unsigned int N, HashTypes HT>
class MetaTable : public UniformLookupTable<TIN,TOUT>
{
protected:
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_UNIFORM_LUT(TIN,TOUT);

  __attribute__((aligned)) std::unique_ptr<polynomial<TOUT,N>[]> m_table;
  TOUT get_table_entry(unsigned int i, unsigned int j) override { return m_table[i].coefs[j]; }
  unsigned int get_num_coefs() override { return N; }

public:
  MetaTable(FunctionContainer<TIN,TOUT> *func_container, UniformLookupTableParameters<TIN> par) :
    UniformLookupTable<TIN,TOUT>(func_container, par) { }

  /* build this table from a file. Everything other than m_table is built by UniformLookupTable */
  MetaTable(FunctionContainer<TIN,TOUT> *func_container, std::string filename, std::string tablename) :
    UniformLookupTable<TIN,TOUT>(func_container, filename)
  {
    std::ifstream file_reader(filename);
    using nlohmann::json;
    json jsonStats;
    file_reader >> jsonStats;

    // check that the names match
    m_name = jsonStats["name"].get<std::string>();
    if(m_name != tablename)
      throw std::invalid_argument("Error while reading " + filename + ": "
          "Cannot build a " + m_name + " from a " + tablename);

    m_table.reset(new polynomial<TOUT,N>[m_numTableEntries]);
    for(unsigned int i=0; i<m_numTableEntries; i++)
      for(unsigned int j=0; j<m_table[i].num_coefs; j++)
        m_table[i].coefs[j] = jsonStats["table"][std::to_string(i)]["coefs"][std::to_string(j)].get<TOUT>();

    // build the transfer function here as well? Could be in the more general LUT file
  }

  /* Provide the most common hash types. The compiler should simplify this when templates are specialized */
  TOUT operator()(TIN x) override
  {
    // TODO figure out a default case
    switch(HT){
    case TAYLOR:
      {
        // nondimensionalized x position
        TOUT dx  = (x-m_minArg);
        TOUT x0r = dx/m_stepSize+0.5;
        // index of previous table entry
        unsigned  x0 = (unsigned) x0r;
        dx -= x0*m_stepSize;

        // general degree horners method, evaluated from the inside out.
        TOUT sum = 0;
        for (int k=N-1; k>0; k--)
          sum = dx*(m_table[x0].coefs[k] + sum);
        return m_table[x0].coefs[0]+sum;
      }
    case HORNER:
      {
        // nondimensionalized x position, scaled by step size
        TOUT dx = (TOUT) m_stepSize_inv*(x-m_minArg);
        // index of previous table entry
        unsigned x0  = (unsigned) dx;
        // value of table entries around x position
        dx -= x0;
        
        // general degree horners method, evaluated from the inside out.
        TOUT sum = 0;
        for (int k=N-1; k>0; k--)
          sum = dx*(m_table[x0].coefs[k] + sum);
        return m_table[x0].coefs[0]+sum;
      }
    }
  }
};

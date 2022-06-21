/*
  TODO explain how this links up with MetaTable & each other implementation
  Intermediate abstract class for LUTs with uniformly spaced grid points.

  This class handles everything required to work with a piecewise
  evaluation implementation. For each such object, we have a
  registry, plus reading & writing of table data to json
  with print_details_json.

  Notes:
  - In the case where (max-min)/stepsize is not an integer then
  the real table max is greater than the user supplied max
*/
#pragma once
#include "FunctionContainer.hpp"
#include "EvaluationImplementation.hpp"
#include "TransferFunctionSinh.hpp"
#include "json.hpp"
#include "config.hpp" // FUNC_USE_SMALL_REGISTRY

#include <map>
#include <memory>
#include <vector>
#include <functional>
#include <fstream>

template <typename TIN>
struct LookupTableParameters
{
  TIN minArg;
  TIN maxArg;
  TIN stepSize;

  // support initializer lists
  LookupTableParameters(TIN min, TIN max, TIN step) :
    minArg(min), maxArg(max), stepSize(step) {}
  LookupTableParameters(){}
};

static constexpr unsigned int alignments[] = {0,1,2,4,4,8,8,8,8};

template <typename TOUT, unsigned int N>
struct alignas(sizeof(TOUT)*alignments[N]) polynomial{
  static const unsigned int num_coefs = N;
  TOUT coefs[N];
};

/* Use to inherit LookupTable's member variables
   without having to put "this->" everywhere. */
#define INHERIT_LUT(TIN,TOUT) \
  using LookupTable<TIN,TOUT>::m_grid; \
  using LookupTable<TIN,TOUT>::m_numIntervals; \
  using LookupTable<TIN,TOUT>::m_numTableEntries; \
  using LookupTable<TIN,TOUT>::m_stepSize; \
  using LookupTable<TIN,TOUT>::m_stepSize_inv; \
  using LookupTable<TIN,TOUT>::m_tableMaxArg


template <typename TIN, typename TOUT = TIN>
class LookupTable : public EvaluationImplementation<TIN,TOUT>
{
protected:
  INHERIT_EVALUATION_IMPL(TIN,TOUT);

  std::unique_ptr<TIN[]> m_grid;  // pointers to grid and evaluation data
  // a polynomial (above) array needs to be provided by each implementation

  // get the ith polynomial's jth coefficient
  virtual TOUT get_table_entry(unsigned int i, unsigned int j)=0;
  // get the number of coefficients used
  virtual unsigned int get_num_coefs()=0;
  virtual std::array<TIN,4> get_transfer_function_coefs()=0;

  unsigned int m_numIntervals;   // sizes of grid and evaluation data
  unsigned int m_numTableEntries;

  TIN  m_stepSize;    // fixed grid spacing (and its inverse)
  TIN  m_stepSize_inv;

  TIN  m_tableMaxArg; // > m_maxArg if (m_maxArg-m_minArg)/m_stepSize is non-integer

public:

  LookupTable(FunctionContainer<TIN,TOUT> *func_container,
      LookupTableParameters<TIN> par) :
    EvaluationImplementation<TIN,TOUT>(func_container->standard_func, "lookup_table")
  {
    /* Base class variables */
    m_minArg = par.minArg; m_maxArg = par.maxArg;

    /*
       Uniform Lookup Table variables
       - corresponds to the grid
    */
    m_stepSize     = par.stepSize;
    m_stepSize_inv = 1.0/m_stepSize;
    m_numIntervals = (unsigned) ceil(m_stepSize_inv*(m_maxArg-m_minArg))+1;
    /*
      If the step size does not exactly divide the arg domain, the max
      arg of the table is set to the nearest value above such that it
      does.
     */
    m_tableMaxArg = m_maxArg;
    if ( m_tableMaxArg < m_minArg+m_stepSize*(m_numIntervals-1) )
      m_tableMaxArg = m_minArg+m_stepSize*(m_numIntervals-1);
    m_grid.reset(new TIN[m_numIntervals]); // TODO maybe only initialize if grid is nonuniform?
  }

  /* Set every generic member variable from a json file */
  LookupTable(FunctionContainer<TIN,TOUT> *func_container,
      std::string filename) :
    EvaluationImplementation<TIN,TOUT>(func_container->standard_func, "lookup_table")
  {
    // Assuming we're still using the same function container as before
    // We'll end up making 2 jsonStats objects because the MetaTable needs
    // to check the table name and construct m_table. Not a big deal
    std::ifstream file_reader(filename);
    using nlohmann::json;
    json jsonStats;
    file_reader >> jsonStats;

    // name checking happens in MetaTable
    m_minArg = jsonStats["minArg"].get<TIN>();
    m_maxArg = jsonStats["maxArg"].get<TIN>();
    m_stepSize = jsonStats["stepSize"].get<TIN>();
    m_stepSize_inv = 1.0/m_stepSize;
    m_order = jsonStats["order"].get<unsigned int>();
    m_dataSize = jsonStats["dataSize"].get<unsigned int>();
    m_numIntervals = jsonStats["numIntervals"].get<unsigned int>();
    m_numTableEntries = jsonStats["numTableEntries"].get<unsigned int>();

    // recompute m_tableMaxArg
    m_tableMaxArg = m_maxArg;
    if ( m_tableMaxArg < m_minArg+m_stepSize*(m_numIntervals-1) )
      m_tableMaxArg = m_minArg+m_stepSize*(m_numIntervals-1);

    // Recompute m_table, m_grid (the array of polynomials), and the
    // transfer function in MetaTable
  }

  virtual ~LookupTable(){};

  /* public access of protected data */
  TIN step_size(){ return m_stepSize; };
  unsigned num_table_entries(){ return m_numTableEntries; };
  unsigned num_intervals(){ return m_numIntervals; };

  void print_details(std::ostream &out)
  {
    out << m_name << " " << m_minArg << " " << m_maxArg << " "
        << m_stepSize << " " << m_numIntervals << " ";
  }

  virtual std::pair<TIN,TIN> arg_bounds_of_interval(unsigned intervalNumber)
  {
    return std::make_pair(m_minArg + intervalNumber*m_stepSize,m_minArg + (intervalNumber+1)*m_stepSize);
  }

  /* Write table data to the provided ostream in the form of json. Virtual in case derived
     classes need to save their own member vars. eg. NonuniformLUTs will store their
     grid and rebuild their transfer function */
  virtual void print_details_json(std::ostream& out)
  {
    // TODO if FunC gets a namespace we should consider renaming our json
    // functions to to_json() and from_json() functions as seen here
    // https://github.com/nlohmann/json
    using nlohmann::json;
    json jsonStats;

    jsonStats["_comment"] = "FunC lookup table data";
    jsonStats["name"] = m_name;
    jsonStats["minArg"] = m_minArg;
    jsonStats["maxArg"] = m_maxArg;
    jsonStats["stepSize"] = m_stepSize;
    jsonStats["order"] = m_order;
    jsonStats["dataSize"] = m_dataSize;
    jsonStats["numIntervals"] = m_numIntervals;
    jsonStats["numTableEntries"] = m_numTableEntries;
    jsonStats["transfer_function_coefs"] = get_transfer_function_coefs();

    // save the polynomial coefs of each lookup table
    // Note: m_order is used as the number of polynomial coefs
    for(unsigned int i=0; i<m_numTableEntries; i++)
      for(unsigned int j=0; j<get_num_coefs(); j++)
        jsonStats["table"][std::to_string(i)]["coefs"][std::to_string(j)] = get_table_entry(i,j);

    out << jsonStats.dump(2) << std::endl;
  }
};

// Legacy func typedef
template <typename TIN, typename TOUT>
using UniformLookupTable = LookupTable<TIN,TOUT>;

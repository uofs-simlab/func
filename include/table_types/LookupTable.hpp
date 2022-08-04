/*
  This class is an interface for any piecewise evaluation
  implementation which we can use with LookupTableFactory.

  Notes:
  - In the case where (max-min)/stepsize is not an integer then
  the real table max is greater than the user supplied max

*/
#pragma once
#include "FunctionContainer.hpp"
#include "EvaluationImplementation.hpp"

#include <memory>
#include <functional>
#include <stdexcept>

namespace func {

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


/* Use to inherit LookupTable's member variables
   without having to put "this->" everywhere. */
#define INHERIT_LUT(TIN,TOUT) \
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

  // a polynomial (defn above) array must be provided by each implementation. 

  unsigned int m_numIntervals;   // sizes of grid and evaluation data
  unsigned int m_numTableEntries;

  TIN  m_stepSize;    // fixed grid spacing (and its inverse)
  TIN  m_stepSize_inv;
  TIN  m_tableMaxArg; // > m_maxArg if (m_maxArg-m_minArg)/m_stepSize is non-integer

public:
  /* Generate a LUT from func_container.
   * TODO How helpful is the error message if func_container == nullptr? */
  LookupTable(FunctionContainer<TIN,TOUT> *func_container,
      LookupTableParameters<TIN> par) :
    EvaluationImplementation<TIN,TOUT>(func_container->standard_func)
  {
    if(func_container->standard_func == nullptr)
      throw std::invalid_argument("function not defined in given FunctionContainer");

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
    m_tableMaxArg = m_minArg+m_stepSize*(m_numIntervals-1);
  }

  // TODO
  // virtual + and * methods for curried LUTs (implementation in MetaTable ideally)
  /* return an approximation to f^(N)(x) */
  //template <unsigned int N>
  //virtual TOUT diff(TIN x) = 0;

  LookupTable() = default;
  virtual ~LookupTable(){};
  // std::unique_ptr m_grid implicitly deletes the copy ctor so we have to explicitly
  // ask for the default copy ctor
  LookupTable(LookupTable&&) = default;

  /* public access of protected data */
  TIN step_size() const { return m_stepSize; };
  unsigned num_table_entries() const { return m_numTableEntries; };
  unsigned num_intervals() const { return m_numIntervals; };

  virtual void print_details(std::ostream &out) override
  {
    out << m_name << " " << m_minArg << " " << m_maxArg << " "
        << m_stepSize << " " << m_numIntervals << " ";
  }

  virtual std::pair<TIN,TIN> arg_bounds_of_interval(unsigned intervalNumber)
  {
    return std::make_pair(m_minArg + intervalNumber*m_stepSize,m_minArg + (intervalNumber+1)*m_stepSize);
  }
};

/* TODO Can we partialially specialize operator() when TOUT = std::unique_ptr<LookupTable<TIN, TOUT2>> to fully evaluate everything?
 * Or potentially, when TOUT is any callable type (just to simplify user syntax) */
} // namespace func

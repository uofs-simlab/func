/*
  A wrapper for a standard func table class. If an argument is outside a
  table's range then the original function is used and the arg is recorded.

  Usage example:
    FailureProofTable<double> failsafe(unique_ptr<UniformLookupTable<double>>(
      new UniformCubicPrecomputedInterpolationTable<double>(&function,0,10,0.0001))
    );
    double val = failsafe(0.87354);

  Notes:
  - ownership of the original table is moved to this class upon construction
  (not a problem since tables don't have move constructors for unique_ptr)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
  - specify the NDEBUG flag to turn off argument recording for args outside
  the table's range.
  - Args outside the table's range are printed to stderr upon table destruction
  (so you could redirect that output to a file with '2> foo.txt')
*/
#pragma once
#include "EvaluationImplementation.hpp"
#include "UniformLookupTable.hpp"

#ifndef NDEBUG
  #include <iostream>
  #include <vector>
  #include <mutex>
  // make sure we don't swallow the semicolon
  #define RECORD_ARG(x)                   \
    do{                                   \
      const std::lock_guard<std::mutex>   \
        lock(m_args_mutex);               \
      m_args.push_back((x));              \
    } while(0)
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
#else
  #define RECORD_ARG(x)
  #define PRINT_ARGS(out)
#endif

template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE>
class FailureProofTable final : public EvaluationImplementation<IN_TYPE,OUT_TYPE> {
  #ifndef NDEBUG
    std::vector<IN_TYPE> m_args;
    std::mutex m_args_mutex;
  #endif

  std::unique_ptr<UniformLookupTable<IN_TYPE,OUT_TYPE>> mp_LUT;
public:
  /* Steal the given LUTs identity TODO remove std::move() ?*/
  FailureProofTable(std::unique_ptr<UniformLookupTable<IN_TYPE,OUT_TYPE>> LUT) :
    mp_LUT(std::move(LUT)), EvaluationImplementation<IN_TYPE,OUT_TYPE>(mp_LUT->function(), "FailureProof" + mp_LUT->name())
  {
    this->m_minArg   = mp_LUT->min_arg();
    this->m_maxArg   = mp_LUT->max_arg();
    this->m_order    = mp_LUT->order();
    this->m_dataSize = mp_LUT->size();
  }

  /* if x isn't contained within the tables bounds,
   * then resort to evaluating the original function */
  OUT_TYPE operator()(IN_TYPE x) override
  {
    // check if x is in the range of the table
    if(x < this->m_minArg || x > this->m_maxArg){
      RECORD_ARG(x);
      return this->mp_func(x);
    }
    return (*(this->mp_LUT))(x);
  }

  void print_details(std::ostream &out) override 
  {
    out << this->m_name << " " << this->m_minArg << " " << this->m_maxArg << " "
        << mp_LUT->step_size() << " " << mp_LUT->num_intervals() << " ";
  }
  
  /* print out any args that were recorded as out of bounds */
  ~FailureProofTable()
  {
    PRINT_ARGS(std::cerr);
  }
};

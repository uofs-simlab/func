/*
  A wrapper for a standard func table class. If an argument is outside a
  table's range then the original function is used and the arg is recorded.

  Usage example:
    UniformFailureProofTable failsafe(unique_ptr<UniformLookupTable>(
      new UniformCubicPrecomputedInterpolationTable(&function,0,10,0.0001))
    );
    double val = failsafe(0.87354);

  Notes:
  - ownership of the original table is moved to this class upon construction
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
  #include <vector>
  #include <mutex>
#endif

class UniformFailureProofTable final : public EvaluationImplementation {
  #ifndef NDEBUG
    std::vector<double> m_args;
    std::mutex m_args_mutex;
  #endif

  // the "UniformLookupTable" part can be changed when 
  // non-uniform tables are added
  std::unique_ptr<UniformLookupTable> mp_LUT;
public:
  UniformFailureProofTable(std::unique_ptr<UniformLookupTable>);
  ~UniformFailureProofTable();
  double operator()(double) override;
  void print_details(std::ostream &out) override;
};

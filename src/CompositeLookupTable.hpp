/*
  A wrapper for several FunC lookup tables.

  Usage example:
    CompositeLookupTable failsafe(unique_ptr<UniformLookupTable>(
      new UniformCubicPrecomputedInterpolationTable(&function,0,10,0.0001))
    );
    double val = failsafe(0.87354);

  Notes:
  - Takes in a variable number of previously constructed lookup tables 
  inside a std::initializer_list
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
  - throws an exception if args are outside table ranges
  - operator() is much faster when repeatedly evaluating
  from the same table's range
*/
#pragma once
#include "EvaluationImplementation.hpp"
#include "UniformLookupTable.hpp"
#include <vector> // store LUTs & discontinuities
#include <memory> // shared_ptr
#include <utility> // std::pair

class CompositeLookupTable final : public EvaluationImplementation {

  std::vector<std::shared_ptr<UniformLookupTable>> mv_LUT;
  std::vector<std::pair<double,double>> mv_discontinuities;
  unsigned int mostRecentlyUsed_idx;

  // find which table to sample from
  double binarySearch(double, int i, int min_idx, int max_idx);
public:
  // call with CompositeLookupTable({EvaluationImplementation, ... });
  CompositeLookupTable(std::initializer_list<std::shared_ptr<UniformLookupTable>>);

  // allow the user to skip the curly braces
  // CompositeLookupTable(EvaluationImplementation, ...)
  //template <typename... SharedPtrToULUT>
  //CompositeLookupTable(SharedPtrToULUT... strings);
  ~CompositeLookupTable();

  // function to return all "holes" in the domain
  std::vector<std::pair<double,double>> discontinuities(){ return mv_discontinuities; }

  double operator()(double) override;

  void print_details(std::ostream &out) override;
};

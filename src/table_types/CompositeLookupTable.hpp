/*
  A wrapper for several FunC lookup tables. Good for approximating piecewise
  functions, and automatic table generation of discontinuous functions
  given all their singularities. Can also be used as a naive non-uniform 
  lookup table. The hash for this table is O(logn) or O(n) depending on
  certain parameters, where n is the number of FunC LUTs.

  TODO move special points to be specific to this class and decide how
  they'll affect table generation.

  Usage example:
    CompositeLookupTable comp_table(unique_ptr<UniformLookupTable>(
      new UniformCubicPrecomputedInterpolationTable(&function,0,10,0.0001))
    );
    double val = comp_table(0.87354);

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

  // This can allow the user to skip the curly braces.
  // CompositeLookupTable(EvaluationImplementation, ...)
  //template <typename... SharedPtrToULUT>
  //CompositeLookupTable(SharedPtrToULUT... strings);
  // Ideal future constructor:
  // CompositeLookupTable(FunctionContainer *func_container, double min, double max, SpecialPoint ... points);
  ~CompositeLookupTable();

  // function to return all "holes" in the domain
  std::vector<std::pair<double,double>> discontinuities(){ return mv_discontinuities; }

  double operator()(double) override;

  void print_details(std::ostream &out) override;
};

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
#include "SpecialPoint.hpp"
#include <vector> // store LUTs, names, special points
#include <memory> // shared_ptr
#include <utility> // std::pair

template <typename IN_TYPE, typename OUT_TYPE>
class CompositeLookupTable final : public EvaluationImplementation<IN_TYPE,OUT_TYPE> {
  // collection of FunC lookup tables used to sample from
  std::vector<std::shared_ptr<UniformLookupTable<IN_TYPE,OUT_TYPE>>> mv_LUT;

  // names of each lookup table used
  std::vector<std::string>  mv_LUT_names;

  // describe function behaviour at the endpoints
  std::vector<SpecialPoint> mv_special_points;

  // index of the last table sampled from
  unsigned int mostRecentlyUsed_idx;
  IN_TYPE smallest_interval;

  // find which table to sample from
  OUT_TYPE binarySearch(IN_TYPE, int, int min_idx, int max_idx);
  OUT_TYPE linearSearch(IN_TYPE, int, bool is_left);

  // Private constructor does all the work but needs to be called with all the 
  // special points in curly braces which is a bit awkward
  CompositeLookupTable(FunctionContainer *func_container, double global_tol, std::initializer_list<SpecialPoint>);

public:
  // public constructors. Build this table with a global tolerance and a collection of special points
  template <typename... SPECIAL_POINTS>
  CompositeLookupTable(FunctionContainer *func_container, double global_tol, SPECIAL_POINTS ... points);

  // or with a vector of n table names, a vector of n step sizes, and a vector of n+1 special points.
  CompositeLookupTable(FunctionContainer *func_container, std::vector<std::string> names,
      std::vector<IN_TYPE> stepSizes, std::vector<SpecialPoint> special_points);

  ~CompositeLookupTable();

  // return the vector of special points in the domain
  std::vector<SpecialPoint<IN_TYPE,OUT_TYPE>> special_points(){ return mv_special_points; };

  OUT_TYPE operator()(IN_TYPE x) override;

  void print_details(std::ostream &out) override;
};

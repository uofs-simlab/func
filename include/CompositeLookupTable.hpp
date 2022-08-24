/*
  A wrapper for several FunC lookup tables. Good for approximating piecewise
  functions, and automatic table generation of discontinuous functions
  given all their singularities. Can also be used as a naive non-uniform
  lookup table. The hash for this table is O(logn) or O(n) depending on
  certain parameters, where n is the number of FunC LUTs.

  Usage example:
    CompositeLookupTable<double> comp_table(&fc, 1e-2, {-1,2},{4,5});
    // or
    // CompositeLookupTable<double> comp_table(&fc, "UniformLinearInterpolationTable", {{-1,2},{4,5}});
    double val = comp_table(0.87354);

  Notes:
  - Takes in a variable number of previously constructed lookup tables
  inside a std::initializer_list
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
  - throws an exception if args are outside table ranges
  - operator() is much faster when repeatedly evaluating
  from the same table's range


  TODO this class should support to/from_json. Could perhaps use the
  unique_ptr<LookupTable> version of from_json in LookupTableFactory

  TODO I think we should use the standard library's map where each LUT is hashed
  based exclusively on its left endpoint. That would be much more maintainable &
  we can still do the "most recently used LUT" thing

  TODO I vote that the special points are only used for table generation. I think it's
  likely too slow & inelegant to have the points do anything nontrivial in operator()
  */
#pragma once
#include "EvaluationImplementation.hpp"
#include "LookupTable.hpp"
#include "LookupTableGenerator.hpp"
#include <vector> // store LUTs, names, special points
#include <memory> // shared_ptr
#include <utility> // std::pair
#include <stdexcept> // domain_error, invalid_argument
#include <string> // to_string

// FUNC_LENIENCE_FACTOR is number of times we try linear search before switching to a binary search
// TODO I doubt we ever need this to not be 0 ...
#define FUNC_LENIENCE_FACTOR 0

namespace func {

/* A subclass used to define function behaviour at table endpoints and breakpoints
 * TODO We should make the pade tables handle +/- infty special values with a special ctor
 * I think "approaches" is a bit pointless because users can just put an if statement in the function's definition
 * */
template <typename TIN, typename TOUT = TIN>
class SpecialPoint
{
protected:
  // coordinate such that (x,y) = (x,f(x))
  std::pair<TIN,TOUT> m_point;

public:
  // use these enums to explain why this point is special
  enum DiscontType { Discont=0, FirstDiscont=1, SecondDiscont=2, ThirdDiscont=3, None=8 };
  enum LimitType { Equals=2, Approaches=0, Inf = 1, NegInf = -1 };
protected:
  DiscontType m_discType;
  LimitType m_limType;

public:
  SpecialPoint(std::pair<TIN,TOUT> pt, DiscontType dt = None, LimitType lt = Equals) :
    m_point(pt), m_discType(dt), m_limType(lt) {}

  // TODO add support for initializer lists
  SpecialPoint(TIN x, TOUT y, DiscontType dt = None, LimitType lt = Equals) :
    SpecialPoint(std::make_pair(x,y), dt, lt) {}


  // public getters
  TIN get_x(){ return m_point.first; }

  TOUT get_y(){
    if(m_limType % 2 != 0)
      return m_limType*std::numeric_limits<TOUT>::max();
    return m_point.second;
  }
  DiscontType discType(){ return m_discType; }
  LimitType limType(){ return m_limType; }
};


template <typename TIN, typename TOUT = TIN>
class CompositeLookupTable final : public EvaluationImplementation<TIN,TOUT> {
  INHERIT_EVALUATION_IMPL(TIN,TOUT);

  // collection of FunC lookup tables used to sample from
  std::vector<std::shared_ptr<LookupTable<TIN,TOUT>>> mv_LUT;

  LookupTableFactory<TIN,TOUT> factory;

  // names of each lookup table used
  std::vector<std::string> mv_LUT_names;

  // describe function behaviour at the endpoints
  std::vector<SpecialPoint<TIN,TOUT>> mv_special_points;

  // index of the last table sampled from
  unsigned int mostRecentlyUsed_idx;
  TIN m_len_smallest_interval;

  // Based on this type of discontinuity in p1 and p2, return a table key
  // NOTE: ASSUMES p1 < p2
  std::string chooseName(SpecialPoint<TIN,TOUT> p1, SpecialPoint<TIN,TOUT> p2);

  // find which table to sample from
  TOUT binarySearch(TIN, int, int min_idx, int max_idx);
  TOUT linearSearchLeft(TIN, int);
  TOUT linearSearchRight(TIN, int);

  // Private constructor does all the work but needs to be called with all the
  // special points in curly braces which is a bit awkward
  CompositeLookupTable(FunctionContainer<TIN,TOUT> *func_container, double global_tol,
      std::initializer_list<SpecialPoint<TIN,TOUT>>);

public:
  // public constructors. Build this table with a global tolerance and a collection of special points
  // Note: This constructor may not converge if you use function roots as special points b/c of the metric
  // we use to quantify relative error
  template <typename... SPECIAL_POINTS>
  CompositeLookupTable(FunctionContainer<TIN,TOUT> *func_container, double global_tol, SPECIAL_POINTS ... points);

  // or with a vector of n table names, a vector of n step sizes, and a vector of n+1 special points. Order
  // determines which tables are used on which subintervals
  CompositeLookupTable(FunctionContainer<TIN,TOUT> *func_container, std::vector<std::string> names,
      std::vector<TIN> stepSizes, std::vector<SpecialPoint<TIN,TOUT>> special_points);

  ~CompositeLookupTable(){}

  TOUT operator()(TIN x) override;

  // return the vector of special points in the domain
  std::vector<SpecialPoint<TIN,TOUT>> special_points(){ return mv_special_points; };

  void print_details(std::ostream &out) override
  {
    out << this->m_name << " ";
    for(auto lut : mv_LUT)
      lut->print_details(out);
  }

  // TODO just use to_json from each lut
  void print_details_json(std::ostream & /* out */) override
  {}

  /* TODO this function is the only advantage we get from using shared_ptr. I don't think this
   * function should be used in 99% of use cases (just use a LUTGenerator instead). 
   *
   * Is this worth keeping? switching to unique_ptr might make life easier */
  std::shared_ptr<LookupTable<TIN,TOUT>> get_table(unsigned int table_idx)
  {
    return mv_LUT[table_idx];
  }
};

/* ---- Class implementation ---- */
template <typename TIN, typename TOUT>
CompositeLookupTable<TIN,TOUT>::CompositeLookupTable(
    FunctionContainer<TIN,TOUT> *func_container,
    std::vector<std::string> names,
    std::vector<TIN> stepSizes,
    std::vector<SpecialPoint<TIN,TOUT>> special_points) :
  EvaluationImplementation<TIN,TOUT>((func_container != nullptr) ? func_container->standard_func : nullptr, "CompositeLookupTable"),
  mv_LUT_names(names), mv_special_points(special_points)
{
  if(m_func == nullptr)
    throw std::invalid_argument("Error in func::CompositeLUT: requires a FunctionContainer.");

  // check if names, stepSizes, and special_points are the right sizes
  if(names.size() != stepSizes.size())
    throw std::invalid_argument("Error in func::CompositeLUT: The " + std::to_string(names.size()) + " given table(s) need(s) "
        "a corresponding stepsize but " + std::to_string(stepSizes.size()) + " stepsizes were given");

  if(names.size() != special_points.size() - 1)
    throw std::invalid_argument("Error in func::CompositeLUT: Need exactly " + std::to_string(names.size() + 1) +
        " SpecialPoints but only " + std::to_string(special_points.size()) + " were given");

  // make sure special_points is ordered (the parameters don't make sense otherwise)
  for(unsigned int i=0; i<special_points.size()-1; i++)
    if(mv_special_points[i].get_x() > mv_special_points[i+1].get_x())
      throw std::invalid_argument("Error in func::CompositeLUT: The x values in the given vector of special points must be ordered "
          "but special_points[" + std::to_string(i) + "].get_x() > special_points[" + std::to_string(i+1) + "].get_x()");

  // naive initial smallest interval
  m_len_smallest_interval = std::numeric_limits<TIN>::max();
  mostRecentlyUsed_idx = names.size()/2;
  this->m_dataSize = 0;
  this->m_order = 0;

  const TIN closeness = std::numeric_limits<TIN>::epsilon();

  /* -- actually build the given tables and update cumulative member vars -- */
  for(unsigned int i=0; i<names.size(); i++){
    // build a table from
    LookupTableParameters<TIN> par;

    // if the discType is "+/- Inf" or "approaches" then we will let m_func take care of that point
    par.minArg = mv_special_points[i].get_x();
    if(mv_special_points[i].limType() != SpecialPoint<TIN,TOUT>::Equals)
      par.minArg += closeness;

    par.maxArg = mv_special_points[i+1].get_x();
    if(mv_special_points[i+1].limType() != SpecialPoint<TIN,TOUT>::Equals)
      par.maxArg -= closeness;

    auto width = par.maxArg - par.minArg;

    par.stepSize = width/ceil(width/stepSizes[i]);
    mv_LUT.push_back(factory.create(mv_LUT_names[i], func_container, par));

    // update the smallest interval and data size
    if(width < m_len_smallest_interval)
      m_len_smallest_interval = width;
    this->m_dataSize += mv_LUT[i]->size();
  }

  // set the global min/max
  this->m_minArg = mv_LUT.front()->min_arg();
  this->m_maxArg = mv_LUT.back()->max_arg();
}


template <typename TIN, typename TOUT>
CompositeLookupTable<TIN,TOUT>::CompositeLookupTable(
    FunctionContainer<TIN,TOUT> *func_container,
    double global_tol,
    std::initializer_list<SpecialPoint<TIN,TOUT>> special_points) :
  EvaluationImplementation<TIN,TOUT>(func_container->standard_func, "CompositeLookupTable"),
  mv_special_points(special_points)
{
  /* TODO decide how the special points affect table generation.
     We want to preserve any discontinuity present in the given function.
     Pade tables for poles? Also have to be conscious of how we approach points
     equality is easy, but the other two will require care */
  // make sure special_points is ordered
  if(special_points.size() < 2)
    throw std::invalid_argument("At least 2 special points must be provided");

  for(unsigned int i=0; i<special_points.size()-1; i++)
    if(mv_special_points[i].get_x() > mv_special_points[i+1].get_x())
      throw std::invalid_argument("The x values in the given special points must be ordered "
          "but special_points[" + std::to_string(i) + "].get_x() > special_points["
          + std::to_string(i+1) + "].get_x()");

  this->m_minArg = mv_special_points.front().get_x();
  this->m_maxArg = mv_special_points.back().get_x();

  // We have everything but the names of tables needed to use our LookupTableGenerator
  for(unsigned int i=0; i < mv_special_points.size()-1; i++){
    TIN ith_min_arg = mv_special_points[i];
    TIN ith_max_arg = mv_special_points[i+1];
    LookupTableGenerator<TIN,TOUT> gen(func_container, ith_min_arg, ith_max_arg);

    // now we just need to decide which of our tables to use
    std::string ith_table_name = chooseName(mv_special_points[i], mv_special_points[i+1]);
    mv_LUT.push_back(gen.generate_by_tol(ith_table_name, global_tol));
  }
}

// Call the above constructor with an initializer list of special points
template <typename TIN, typename TOUT>
template <typename ... SPECIAL_POINTS>
CompositeLookupTable<TIN,TOUT>::CompositeLookupTable(
    FunctionContainer<TIN,TOUT> *func_container,
    double global_tol,
    SPECIAL_POINTS ... points) :
  CompositeLookupTable(func_container,global_tol,{points ...})
{}

template <typename TIN, typename TOUT>
std::string CompositeLookupTable<TIN,TOUT>::chooseName(SpecialPoint<TIN,TOUT> /* p1 */, SpecialPoint<TIN,TOUT> /* p2 */)
{
  /* Decide what table to use on this interval
     TODO because this table's operator() uses so many conditional statements
     we won't save time by using lower order tables with tons of subintervals.
     That would require much more memory and there's no way their coefs will stay in cache for long.
     Thus, use the highest order precomputed interpolation table less than or equal to cubic, and use
     Pade tables if -DUSE_ARMADILLO is defined */
  return std::string("UniformCubicPrecomputedInterpolationTable");
//  if(p1.discType() > 2 && p2.discType() > 2)
//    return "UniformCubicPrecomputedInterpolationTable";
//  // at least one of the points have discType() <= 2
//  if((p1.limType() == Equals || p1.limType() == Approaches) && (p2.limType() == Equals || p2.limType() == Approaches))
//    return "UniformCubicPrecomputedInterpolationTable";
//
//  // at least one of the points explodes
//  if((p1.discType() == 2 && p2.discType() > 2) || (p1.discType() > 2 && p2.discType() == 2))
//    return "UniformQuadraticPrecomputedInterpolationTable";
//  // at least one of the points have discType() <= 1
//  if((p1.discType() == 1 && p2.discType() > 2) || (p1.discType() > 2 && p2.discType() == 1)){
//#ifdef USE_ARMADILLO
//    return "UniformPadeTable<2,1>";
//#endif
//    return "UniformLinearPrecomputedInterpolationTable";
//  }
}

// TODO reference special points for interval bounds
template <typename TIN, typename TOUT>
TOUT CompositeLookupTable<TIN,TOUT>::binarySearch(TIN x, int i, int min_idx, int max_idx)
{
  // Binary search for the correct interval starting with i (most
  // recently used table index). Best for seemly random table evaluations.
  if(x < mv_LUT[i]->min_arg()){
    if(i == min_idx)
      return m_func(x); // none of our LUTs have x in their domain.
    return binarySearch(x, (int)(i+min_idx)/2, min_idx, i);
  }
  else if(x > mv_LUT[i]->max_arg()){
    if(i == max_idx)
      return m_func(x); // none of our LUTs have x in their domain.
    return binarySearch(x, (int)(i+max_idx)/2, i, max_idx);
  }

  // we're in the right interval. Check the special points
  //if(abs(x - mv_LUT[i]->min_arg()) < std::numeric_limits<TIN>::epsilon())
  //  return mv_special_points[i].get_y(); // use the left special point's value
  //if(abs(x - mv_LUT[i]->max_arg()) < std::numeric_limits<TIN>::epsilon())
  //  return mv_special_points[i+1].get_y(); // use the right special point's value

  // return the LUT value
  return (*mv_LUT[i])(x);
}

template <typename TIN, typename TOUT>
TOUT CompositeLookupTable<TIN,TOUT>::linearSearchLeft(TIN x, int i)
{
  // Assuming this function is called with all the correct parameters and
  // mv_LUT[i] is defined
  if(x < mv_LUT[i]->min_arg())
    return linearSearchLeft(x, i-1);
  if(abs(x - mv_LUT[i]->min_arg()) < std::numeric_limits<TIN>::epsilon())
    return mv_special_points[i].get_y(); // use the left special point's value
  return (*mv_LUT[i])(x);
}

template <typename TIN, typename TOUT>
TOUT CompositeLookupTable<TIN,TOUT>::linearSearchRight(TIN x, int i)
{
  // Assuming this function is called with all the correct parameters and
  // mv_LUT[i] is defined
  if(x > mv_LUT[i]->max_arg())
    return linearSearchRight(x, i+1);
  if(abs(x - mv_LUT[i]->max_arg()) < std::numeric_limits<TIN>::epsilon())
    return mv_special_points[i+1].get_y(); // use the right special point's value
  return (*mv_LUT[i])(x);
}

template <typename TIN, typename TOUT>
TOUT CompositeLookupTable<TIN,TOUT>::operator()(TIN x)
{
  // If x is close, use a linear search. Otherwise, use a binary search
  std::shared_ptr<LookupTable<TIN,TOUT>> recentTable = mv_LUT[mostRecentlyUsed_idx];

  if(x < recentTable->min_arg() - FUNC_LENIENCE_FACTOR*m_len_smallest_interval){
    // x is far away, do binary search on the left
    return binarySearch(x, (int) mostRecentlyUsed_idx/2, 0, mostRecentlyUsed_idx);
  }
  else if(x < recentTable->min_arg()){
    // x is near, do a linear search on the left
    return linearSearchLeft(x, mostRecentlyUsed_idx);
  }
  else if(x < recentTable->max_arg()){
    // x is here
    return (*recentTable)(x);
  }
  else if(x > recentTable->max_arg() + FUNC_LENIENCE_FACTOR*m_len_smallest_interval){
    // x is near, do a linear search on the right
    return linearSearchRight(x, mostRecentlyUsed_idx);
  }
  // x is far, do a binary search on the right
  return binarySearch(x, (int)(mostRecentlyUsed_idx+mv_LUT.size()-1)/2,
      mostRecentlyUsed_idx, mv_LUT.size()-1);
}

} // namespace func

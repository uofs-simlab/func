/*
  A wrapper for several FunC lookup tables. Good for approximating piecewise
  functions, and automatic table generation of discontinuous functions
  given all their singularities. Can also be used as a naive non-uniform 
  lookup table. The hash for this table is O(logn) or O(n) depending on
  certain parameters, where n is the number of FunC LUTs.

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
#include "UniformLookupTableGenerator.hpp"
#include <vector> // store LUTs, names, special points
#include <memory> // shared_ptr
#include <utility> // std::pair
#include <stdexcept> // domain_error, invalid_argument
#include <string> // to_string


/* A subclass used to define function behaviour at table endpoints and breakpoints */
template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE>
class SpecialPoint
{
  // x,y coordinate
  std::pair<IN_TYPE,OUT_TYPE> m_point;

  // explain why this point is special
  enum DiscontType { None=-1, Discont=0, FirstDiscont=1, SecondDiscont=2, ThirdDiscont=3 };
  enum LimitType { Equals, Approaches, Inf };
  DiscontType m_discType;
  LimitType m_limType;

public:
  SpecialPoint(IN_TYPE x, OUT_TYPE y, DiscontType dt, LimitType lt) :
    m_point(std::make_pair(x,y)), m_discType(dt), m_limType(lt) {}

  SpecialPoint(std::pair<IN_TYPE,OUT_TYPE> pt, DiscontType dt, LimitType lt) :
    m_point(pt), m_discType(dt), m_limType(lt) {}

  // public getters
  std::pair<IN_TYPE,OUT_TYPE> point(){ return m_point; }
  DiscontType discType(){ return m_discType; }
  LimitType limType(){ return m_limType; }
};


template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE>
class CompositeLookupTable final : public EvaluationImplementation<IN_TYPE,OUT_TYPE> {
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);

  // collection of FunC lookup tables used to sample from
  std::vector<std::shared_ptr<UniformLookupTable<IN_TYPE,OUT_TYPE>>> mv_LUT;

  // names of each lookup table used
  std::vector<std::string> mv_LUT_names;

  // describe function behaviour at the endpoints
  std::vector<SpecialPoint<IN_TYPE,OUT_TYPE>> mv_special_points;

  // index of the last table sampled from
  unsigned int mostRecentlyUsed_idx;
  IN_TYPE m_len_smallest_interval;

  // find which table to sample from
  OUT_TYPE binarySearch(IN_TYPE, int, int min_idx, int max_idx);
  OUT_TYPE linearSearch(IN_TYPE, int, bool is_left);

  // Private constructor does all the work but needs to be called with all the 
  // special points in curly braces which is a bit awkward
  CompositeLookupTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, double global_tol,
      std::initializer_list<SpecialPoint<IN_TYPE,OUT_TYPE>>);

public:
  // public constructors. Build this table with a global tolerance and a collection of special points
  template <typename... SPECIAL_POINTS>
  CompositeLookupTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, double global_tol, SPECIAL_POINTS ... points);

  // or with a vector of n table names, a vector of n step sizes, and a vector of n+1 special points. Order
  // determines which tables are used on which subintervals
  CompositeLookupTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, std::vector<std::string> names,
      std::vector<IN_TYPE> stepSizes, std::vector<SpecialPoint<IN_TYPE,OUT_TYPE>> special_points);

  ~CompositeLookupTable(){}

  // return the vector of special points in the domain
  std::vector<SpecialPoint<IN_TYPE,OUT_TYPE>> special_points(){ return mv_special_points; };

  OUT_TYPE operator()(IN_TYPE x) override;

  void print_details(std::ostream &out) override;
  // TODO Add a way to print details about each subinterval
};

/* ---- Class implementation ---- */

// TODO template or enum specifying pure linear search

template <typename IN_TYPE, typename OUT_TYPE>
inline CompositeLookupTable<IN_TYPE,OUT_TYPE>::CompositeLookupTable(
    FunctionContainer<IN_TYPE,OUT_TYPE> *func_container, 
    std::vector<std::string> names,
    std::vector<IN_TYPE> stepSizes,
    std::vector<SpecialPoint<IN_TYPE,OUT_TYPE>> special_points) :
  EvaluationImplementation<IN_TYPE,OUT_TYPE>(func_container->standard_func, "CompositeLookupTable"),
  mv_LUT_names(names), mv_special_points(special_points)
{
  // check if names, stepSizes, and special_points are the right sizes
  if(names.size() != stepSizes.size())
    throw std::invalid_argument("The " + std::to_string(names.size()) + " given table(s) need(s) "
        "a corresponding stepsize but " + std::to_string(stepSizes.size()) + " stepsizes were given");

  if(names.size() != special_points.size() + 1)
    throw std::invalid_argument("Function behaviour for the " + std::to_string(names.size() + 1) + 
        " breakpoints and endpoints need to be defined with SpecialPoints but "
        "only " + std::to_string(special_points.size()) + "SpecialPoints were given");

  // make sure special_points is ordered (the parameters don't make sense otherwise)
  for(unsigned int i=0; i<names.size(); i++)
    if(mv_special_points[i].point().first > mv_special_points[i+1].point().first)
      throw std::invalid_argument("The x values in the given vector of special points must be ordered "
          "but special_points[" + std::to_string(i) + "].point().first > special_points["
          + std::to_string(i+1) + "].point().first");

  // naive initial smallest interval
  m_len_smallest_interval = std::numeric_limits<IN_TYPE>::max();
  mostRecentlyUsed_idx = names.size()/2;
  this->m_dataSize = 0;
  this->m_order = 0;

  /* -- actually build the given tables and update cumulative member vars -- */
  for(unsigned int i=0; i<names.size(); i++){
    // build a table from 
    UniformLookupTableParameters<IN_TYPE> par;
    par.minArg = mv_special_points[i].point().first;
    par.maxArg = mv_special_points[i+1].point().first;
    par.stepSize = stepSizes[i];
    mv_LUT.push_back(UniformLookupTableFactory<IN_TYPE,OUT_TYPE>::Create(mv_LUT_names[i], func_container, par));

    // update the smallest interval and data size
    if(par.maxArg - par.minArg < m_len_smallest_interval)
      m_len_smallest_interval = par.maxArg - par.minArg;
    this->m_dataSize += mv_LUT[i]->size();
  }

  // set the global min/max
  this->m_minArg = mv_LUT.front()->min_arg();
  this->m_maxArg = mv_LUT.back()->max_arg();
}


template <typename IN_TYPE, typename OUT_TYPE>
inline CompositeLookupTable<IN_TYPE,OUT_TYPE>::CompositeLookupTable(
    FunctionContainer<IN_TYPE,OUT_TYPE> *func_container,
    double global_tol,
    std::initializer_list<SpecialPoint<IN_TYPE,OUT_TYPE>>)
{
  // TODO how do special points affect table generation
  // Will want to preserve any discontinuity present in the given function
  // Pade tables for explosions? Also have to be conscious of how we approach points
  // equality is easy, but the other two will require care
}

// Call the above constructor with an initializer list of special points
template <typename IN_TYPE, typename OUT_TYPE>
template <typename ... SPECIAL_POINTS>
inline CompositeLookupTable<IN_TYPE,OUT_TYPE>::CompositeLookupTable(
    FunctionContainer<IN_TYPE,OUT_TYPE> *func_container,
    double global_tol,
    SPECIAL_POINTS ... points) :
  CompositeLookupTable(func_container,global_tol,{points ...})
{}

template <typename IN_TYPE, typename OUT_TYPE>
inline OUT_TYPE CompositeLookupTable<IN_TYPE,OUT_TYPE>::binarySearch(IN_TYPE x, int i, int min_idx, int max_idx)
{
  // Binary search for the correct interval starting with i (most
  // recently used table index). Best for seemly random table evaluations.
  if(x < mv_LUT[i]->min_arg()){
    if(i == min_idx)
      throw std::domain_error(std::string("Composite table undefined for x=") +
          std::to_string(x));
    return binarySearch(x, (int)(i+min_idx)/2, min_idx, i);
  }
  else if(x > mv_LUT[i]->max_arg()){
    if(i == max_idx)
      throw std::domain_error(std::string("Composite table undefined for x=") +
          std::to_string(x));
    return binarySearch(x, (int)(i+max_idx)/2, i, max_idx);
  }
  return (*mv_LUT[i])(x);
}

template <typename IN_TYPE, typename OUT_TYPE>
inline OUT_TYPE CompositeLookupTable<IN_TYPE,OUT_TYPE>::linearSearch(IN_TYPE x, int i, bool is_left)
{
  // Assuming this function is called with all the correct parameters and 
  // mv_LUT[i] is defined
  if(is_left){
    if(x < mv_LUT[i]->min_arg())
      return linearSearch(x, i-1, true);
    return (*mv_LUT[i])(x);
  }
  if(x > mv_LUT[i]->max_arg())
    return linearSearch(x, i+1, false);
  return (*mv_LUT[i])(x);
}

template <typename IN_TYPE, typename OUT_TYPE>
inline OUT_TYPE CompositeLookupTable<IN_TYPE,OUT_TYPE>::operator()(IN_TYPE x)
{
  // If x is close, use a linear search. Otherwise, use a binary search
  std::shared_ptr<UniformLookupTable<IN_TYPE,OUT_TYPE>> recentTable = mv_LUT[mostRecentlyUsed_idx];

  if(x < recentTable->min_arg() - 2*m_len_smallest_interval){
    // x is far away, do binary search on the left
    return binarySearch(x, (int) mostRecentlyUsed_idx/2, 0, mostRecentlyUsed_idx);
  }
  else if(x < recentTable->min_arg()){
    // x is near, do a linear search on the left
    return linearSearch(x, mostRecentlyUsed_idx, true);
  }
  else if(x < recentTable->max_arg()){
    // x is here
    return (*recentTable)(x);
  }
  else if(x > recentTable->max_arg() + 2*m_len_smallest_interval){
    // x is near, do a linear search on the right
    return linearSearch(x, mostRecentlyUsed_idx, false);
  }
  // x is far, do a binary search on the right
  return binarySearch(x, (int)(mostRecentlyUsed_idx+mv_LUT.size()-1)/2,
      mostRecentlyUsed_idx, mv_LUT.size()-1);
}

template <typename IN_TYPE, typename OUT_TYPE>
inline void CompositeLookupTable<IN_TYPE, OUT_TYPE>::print_details(std::ostream &out)
{
  out << this->m_name << " ";
  for(auto lut : mv_LUT)
    lut->print_details(out);
}

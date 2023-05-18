/*
  A wrapper for several FunC lookup tables. Good for approximating piecewise
  functions, and automatic table generation of discontinuous functions
  given all their singularities. Can also be used as a naive non-uniform
  lookup table. The hash for this table is O(logn) or O(n) depending on
  certain parameters, where n is the number of FunC LUTs.

  Usage example:
    CompositeLookupTable<double> comp_table(&fc, 1e-2, {-1,2},{4,5});
    // or
    // std::vector<std::string> names {"UniformLinearInterpolationTable", "UniformCubicInterpolationTable"};
    // std::vector<double> steps {0.15,0.3};
    // std::vector<func::SpecialPoint<double>> points {{-1.0,f(-1.0)}, {0.0,f(0.0)}, {1.0,f(1.0)}};
    // func::CompositeLookupTable<double> comp_table(&fc, names, steps, points);

    double val = comp_table(0.87354);

  Notes:
  - Takes in a variable number of previously constructed lookup tables
  inside a std::initializer_list
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
  - operator() remembers the most recently used LUT and skips the binary search when repeatedly evaluating from the same table's range


  TODO this class should support to/from_json. We can use the
  unique_ptr<LookupTable> version of from_json in LookupTableFactory
  */
#pragma once
#include "LookupTable.hpp"
#include "LookupTableGenerator.hpp"
#include <map> // store LUTs
#include <vector> // arguments to ctor are vectors
#include <memory> // shared_ptr
#include <utility> // std::tuple
#include <stdexcept> // invalid_argument
#include <numeric>


namespace func {

template <typename TIN, typename TOUT = TIN>
class CompositeLookupTable final : public LookupTableTIN,TOUT> {
  // collection of FunC lookup tables used to sample from
  std::map<TIN,std::shared_ptr<LookupTable<TIN,TOUT>>> m_lutmap;
  std::shared_ptr<LookupTable<TIN,TOUT>> recentLUT;
  std::function<TOUT(TIN)> m_f;
  //LookupTableFactory<TIN,TOUT> factory; TODO maybe useful for building from a file...

public:
  /* vector of n table names, a vector of n step sizes, and a vector of n pairs.
   * Order of names determines which tables are used on which subintervals */
  CompositeLookupTable(const FunctionContainer<TIN,TOUT>& func_container, const std::vector<std::tuple<std::string,TIN,TIN,TIN>>& name_l_r_steps);
  //CompositeLookupTable(const FunctionContainer<TIN,TOUT>& func_container, const std::vector<std::tuple<std::string,TIN,TIN,TIN,TIN>>& name_l_r_atol_rtols);
  ~CompositeLookupTable(){}

  TOUT operator()(TIN x) final
  {
    /* check if x is within the bounds of the most recently used LUT */
    if((recentLUT.min_arg() < x) && (x < recentLUT.max_arg()))
      return (*recentLUT)(x);

    /* Find the LUT whose right endpoint is strictly greater than x
     * TODO is this problematic for the max arg of the compositeLUT? Probably..... */
    auto iter = m_lutmap(x).upper_bound();
    if(iter != m_lutmap.end()){
      auto lut = iter->first;
      if(lut.min_arg() < x){
        recentLUT = lut;
        return (*lut)(x);
      }else{
        return m_func(x);
      }
    }
  }

  std::shared_ptr<LookupTable<TIN,TOUT>> get_table(TIN x){ return m_lutmap(x).lower_bound->first; }
};

/* ---- Class implementation ---- */
template <typename TIN, typename TOUT>
CompositeLookupTable<TIN,TOUT>::CompositeLookupTable(const FunctionContainer<TIN,TOUT>& func_container, std::vector<std::tuple<std::string,TIN,TIN,TIN>> name_l_r_steps) :
  function(func_container.standard_func)
{
  if(m_func == nullptr)
    throw std::invalid_argument("Error in func::CompositeLUT: requires a func::FunctionContainer.");

  // TODO check that each LUT is over distinct intervals

  // TODO there isn't really a great way to set order
  this->m_order = 0;

  /* -- actually build the given tables and update cumulative member vars -- */
  for(auto name_l_r_step : name_l_r_steps){
    // build each given LUT
    auto name  = std::get<0>(name_l_r_step);
    auto left  = std::get<1>(name_l_r_step);
    auto right = std::get<2>(name_l_r_step);
    auto step  = std::get<3>(name_l_r_step);
    LookupTableGenerator<TIN,TOUT> gen(func_container, left, right);
    m_lutmap.insert(right, gen.generate_by_step(name, step));
  }

  // set the global min/max
  this->m_dataSize = std::accumulate();
  this->m_minArg = m_lutmap.begin()->first;
  this->m_maxArg = m_lutmap.end()->second->max_arg();
}

} // namespace func

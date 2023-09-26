/*
  A wrapper for several FunC lookup tables. Good for approximating a single function with multiple LUTs.
  Can also be used as a naive non-uniform lookup table. The hash for this table is O(logn) where n is the number of LUTs.
  Each LUT is considered a single "subinterval"

  Usage example:

    double val = comp_table(0.87354);

  Notes:
  - Takes in a variable number of previously constructed lookup tables
  inside a std::initializer_list
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
  - operator() caches the most recently used LUT and skips the binary search when repeatedly evaluating from the same table's range
  - 


  TODO this class should support to/from_json. We could use the unique_ptr<LookupTable> version of from_json in LookupTableFactory
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
class CompositeLookupTable final : public LookupTable<TIN,TOUT> {
  // collection of FunC lookup tables used to sample from
  using lookup_map_t = std::map<TIN,std::shared_ptr<LookupTable<TIN,TOUT>>>;
  lookup_map_t m_lutmap;
  mutable std::shared_ptr<LookupTable<TIN,TOUT>> recentLUT;
  std::function<TOUT(TIN)> m_fun;
  //LookupTableFactory<TIN,TOUT> factory; TODO might be useful for building a CompositeTable from a file...

public:
  /* vector of n table names, a vector of n step sizes, and a vector of n pairs.
   * Order of names determines which tables are used on which subintervals */
  CompositeLookupTable(const FunctionContainer<TIN,TOUT>& func_container, const std::vector<std::tuple<std::string,TIN,TIN,TIN>>& name_l_r_steps) :
    m_fun(func_container.standard_fun)
  {
    if(m_fun == nullptr)
      throw std::invalid_argument("Error in func::CompositeLUT: given a null func::FunctionContainer.");

    /* TODO we could check that each LUT is over distinct intervals, but that'll probably require setting an order for name_l_r_steps
     * and it's probably not that bad if the output from operator() is not unique for each input */

    /* -- build each of the given LUTs -- */
    for(auto&& name_l_r_step : name_l_r_steps){
      std::string name; TIN left, right, step;
      std::tie(name, left, right, step) = name_l_r_step;

      double N = std::ceil((right-left)/step);
      LookupTableGenerator<TIN,TOUT> gen(func_container, left, right);
      m_lutmap.insert({right, std::move(gen.generate_by_step(name, (right-left)/N))});
    }
    /* initialize the cached LUT */
    recentLUT = m_lutmap.begin()->second;
  }

  CompositeLookupTable(const FunctionContainer<TIN,TOUT>& func_container, const std::vector<std::tuple<std::string,TIN,TIN,TIN,TIN>>& name_l_r_atol_rtols) :
    m_fun(func_container.standard_fun)
  {
    if(m_fun == nullptr)
      throw std::invalid_argument("Error in func::CompositeLUT: given a null func::FunctionContainer.");

    /* -- build each of the given LUTs -- */
    for(auto&& name_l_r_atol_rtol : name_l_r_atol_rtols){
      std::string name; TIN left, right, atol, rtol;
      std::tie(name, left, right, atol, rtol) = name_l_r_atol_rtol;

      LookupTableGenerator<TIN,TOUT> gen(func_container, left, right);
      unsigned int N = gen.generate_by_tol(name, atol, rtol)->num_subintervals();
      m_lutmap.insert({right, std::move(gen.generate_by_step(name, (right-left)/static_cast<double>(N) ))});
    }
    /* initialize the cached LUT */
    recentLUT = m_lutmap.begin()->second;
  }

  ~CompositeLookupTable(){}

  TOUT operator()(TIN x) const final
  {
    /* check if x is within the bounds of the most recently used LUT */
    if((recentLUT->min_arg() < x) && (x < recentLUT->max_arg())) return (*recentLUT)(x);

    /* Find the LUT whose right endpoint is strictly greater than x
     * TODO is this problematic for the max arg of the compositeLUT? */
    auto iter = m_lutmap.upper_bound(x);
    if(iter != m_lutmap.end()){
      auto lut = iter->second;
      if(lut->min_arg() < x){
        recentLUT = lut;
        return (*lut)(x);
      }
    }
    return m_fun(x);
  }

  std::string name() const final { return "CompositeLookupTable"; }
  TIN min_arg() const final { return m_lutmap.begin()->second->min_arg(); }
  TIN max_arg() const final { return m_lutmap.rbegin()->first; }
  unsigned int order() const final { return 0u; } // TODO maybe return min order?

  /* sum the sizes of each LookupTable */
  std::size_t size() const final {
    return std::accumulate(m_lutmap.begin(), m_lutmap.end(), std::size_t{0},
      [](std::size_t total, const typename lookup_map_t::value_type& L)
      { return total + L.second->size(); });
  }

  unsigned int num_subintervals() const final {
    return std::accumulate(m_lutmap.begin(), m_lutmap.end(), 0u,
      [](unsigned int total, const typename lookup_map_t::value_type& L)
      { return total + L.second->num_subintervals(); });
  }

  /* return the min step size of each LookupTable */
  TIN step_size() const final {
    return std::min_element(m_lutmap.begin(), m_lutmap.end(),
      [](const typename lookup_map_t::value_type& L1, const typename lookup_map_t::value_type& L2)
      { return L1.second->step_size() <= L2.second->step_size(); })->second->step_size(); }

  std::pair<TIN,TIN> bounds_of_subinterval(unsigned int intervalNumber) const final {
    long int N = intervalNumber; // need a signed integer to safely do subtraction
    auto it = m_lutmap.begin();
    for(; it != m_lutmap.end(); ++it){
      if(N - it->second->num_subintervals() < 0)
        break;
      N -= it->second->num_subintervals();
    }
    if(it == m_lutmap.end())
      throw std::invalid_argument(std::string("Error in func::CompositeLookupTable: requested intervalNumber is ") + std::to_string(N) +
          " which is greater than num_subintervals = " + std::to_string(num_subintervals()));
    return it->second->bounds_of_subinterval(N);
  }

  void print_json(std::ostream& out) const final { (void) out; /* TODO call to_json() */ };
  std::shared_ptr<LookupTable<TIN,TOUT>> get_table(TIN x){ return m_lutmap.lower_bound(x)->second; }
};

} // namespace func

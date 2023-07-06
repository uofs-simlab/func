/*
  MetaTable handles any _piecewise polynomial based_ interpolation:

  Note:
  - In the case where (max-min)/stepsize is not an integer then
  the actual maximum allowable argument is greater than max
  - if stepSize divides max-min exactly then operator(max) will be out of array bounds
  so every implementation has an extra (unnecessary in most cases) table entry.

  Info on template parameters:

  N = number of coefficients used in underlying piecewise polynomials
  Provided Horner's method which is the most common table evaluation method in FunC

  UNIFORM: Every subinterval is the same length so the hash is super fast; however,
  many more subintervals might be needed
  NONUNIFORM: Use a transfer function to create a nonuniform grid with an O(1) hash.
*/
#pragma once
#include "LookupTable.hpp"
#include "TransferFunction.hpp"
#include "Polynomial.hpp"

#include <array>
#include <stdexcept>
#include "json.hpp"

/* Give inheriting classes access to member variables without
   having to use "this->" excessively. These "using" statements must have protected access */
#define INHERIT_META(N,TIN,TOUT,GT) \
  using MetaTable<N,TIN,TOUT,GT>::m_order; \
  using MetaTable<N,TIN,TOUT,GT>::m_name; \
  using MetaTable<N,TIN,TOUT,GT>::m_dataSize; \
  using MetaTable<N,TIN,TOUT,GT>::m_minArg; \
  using MetaTable<N,TIN,TOUT,GT>::m_maxArg; \
  using MetaTable<N,TIN,TOUT,GT>::m_numIntervals; \
  using MetaTable<N,TIN,TOUT,GT>::m_numTableEntries; \
  using MetaTable<N,TIN,TOUT,GT>::m_stepSize; \
  using MetaTable<N,TIN,TOUT,GT>::m_stepSize_inv; \
  using MetaTable<N,TIN,TOUT,GT>::m_tableMaxArg; \
  using MetaTable<N,TIN,TOUT,GT>::m_table; \
  using MetaTable<N,TIN,TOUT,GT>::m_transferFunction

/* Parallelization macro. Notes:
 * - we could align m_table to sizeof(TOUT) to get some more speedup? */
#define FUNC_BUILDPAR _Pragma("omp parallel for schedule(dynamic)")

namespace func {

static constexpr unsigned int alignments[] = {0,1,2,4,4,8,8,8,8,16,16,16,16,16,16,16,16};

enum class GridTypes {UNIFORM, NONUNIFORM};

template <GridTypes GT>
inline std::string grid_type_to_string() {
  FUNC_IF_CONSTEXPR(GT == GridTypes::UNIFORM){
    return "Uniform";
  }else{
    return "NonUniform";
  } 
}

template <unsigned int N, typename TIN, typename TOUT=TIN, GridTypes GT=GridTypes::UNIFORM>
class MetaTable : public LookupTable<TIN,TOUT>
{
protected:
  std::string m_name;             // name of implementation type
  TIN m_minArg, m_maxArg;         // bounds of evaluation
  TIN m_stepSize, m_stepSize_inv; // fixed grid spacing (and its inverse)
  TIN m_tableMaxArg; // > m_maxArg if (m_maxArg-m_minArg)/m_stepSize is non-integer

  unsigned int m_order;           // order of accuracy of implementation
  std::size_t m_dataSize;        // size of relevant data for impl evaluation
  unsigned int m_numIntervals;    // = (m_tableMaxArg - m_minArg)/m_stepSize;
  unsigned int m_numTableEntries; // length of m_table (usually = m_numIntervals + 1)
  __attribute__((aligned)) std::unique_ptr<polynomial<TOUT,N>[]> m_table; // holds polynomials coefficients
  TransferFunction<TIN> m_transferFunction; // used to make nonuniform grids (default constructable)

public:
  /* using a std::unique_ptr member variables implicitly deletes the default
   * move ctor so we must explicitly ask for the default move ctor */
  MetaTable(MetaTable&&) = default;
  MetaTable() = default;

  /* Set every generic member variable from a json file */
  MetaTable(const FunctionContainer<TIN,TOUT>& func_container, const LookupTableParameters<TIN>& par) :
    m_minArg(par.minArg), m_maxArg(par.maxArg), m_stepSize(par.stepSize) {
    /* If the step size does not exactly divide the arg domain, the max arg of the table is set
     * to the nearest value above such that it does. */
    if(m_stepSize <= static_cast<TIN>(0.0))
      throw std::invalid_argument("func::MetaTable was given a nonpositive stepSize. stepSize must be positive.");

    m_stepSize_inv = static_cast<TIN>(1.0)/m_stepSize;
    m_numIntervals = static_cast<unsigned>(ceil(m_stepSize_inv*(m_maxArg-m_minArg)));
    m_tableMaxArg = m_minArg+m_stepSize*m_numIntervals; // always >= m_maxArg

    // We need a valid FunctionContainer to generate any LUT
    if(func_container.standard_fun == nullptr)
      throw std::invalid_argument("Error in func::MetaTable. Function not defined in given FunctionContainer");

    /* build the transfer function for nonuniform grids
     * TODO could make transfer function coefficients a field in LookupTableParameters */
    FUNC_IF_CONSTEXPR(GT != GridTypes::UNIFORM)
      m_transferFunction = TransferFunction<TIN>(func_container,m_minArg,m_tableMaxArg,m_stepSize);
  }

  /* build this table from a json file */
  MetaTable(const nlohmann::json& jsonStats) {
    if(jsonStats.empty())
      throw std::invalid_argument("Error in func::MetaTable: The provided nlohmann::json is empty");
    from_json(jsonStats, *this);
  }

  /* public access to protected member vars */
  std::string name() const final { return m_name; }
  TIN min_arg() const final { return m_minArg; }
  TIN max_arg() const final { return m_maxArg; }
  TIN tablemax_arg() const { return m_tableMaxArg; }
  unsigned int order() const final { return m_order; }
  std::size_t size() const final { return m_dataSize; }
  unsigned int num_subintervals() const final { return m_numIntervals; };
  TIN step_size() const final { return m_stepSize; };
  std::pair<TIN,TIN> bounds_of_subinterval(unsigned intervalNumber) const final {
    FUNC_IF_CONSTEXPR(GT == GridTypes::UNIFORM)
      return std::make_pair(m_minArg + intervalNumber*m_stepSize,m_minArg + (intervalNumber+1)*m_stepSize);
    else
      return std::make_pair(m_transferFunction(m_minArg + intervalNumber*m_stepSize),
          m_transferFunction(m_minArg + (intervalNumber+1)*m_stepSize));
  }
  void print_json(std::ostream& out) const final {
    nlohmann::json jsonStats = *this; // call to_json(jsonStats, this)
    out << jsonStats.dump(2) << std::endl;
  };

  unsigned int num_table_entries() const { return m_numTableEntries; };
  unsigned int ncoefs_per_entry() const { return N; }
  TOUT table_entry(unsigned int i, unsigned int j) const { return m_table[i].coefs[j]; }
  std::array<TIN,4> transfer_function_coefs() const { return m_transferFunction.get_coefs(); }


  /* give the nlohmann json functions access to every member variable (must use SFINAE because we want people
   * to construct LUTs over arbitrary types regardless of whether they have to/from json functions) */
  template <unsigned int N1, typename TIN1, typename TOUT1, GridTypes GT1,
         typename std::enable_if<std::is_constructible<nlohmann::json,TIN1 >::value && 
                                 std::is_constructible<nlohmann::json,TOUT1>::value, bool>::type>
  friend void from_json(const nlohmann::json& jsonStats, MetaTable<N1,TIN1,TOUT1,GT1>& lut);


  /* find which polynomial p_k to evaluate. Also, each p_k:[0,1]->R so we must set dx=(x-x_k)/(x_{k+1}-x_k) */
  template <GridTypes GT1, typename std::enable_if<GT1 == GridTypes::UNIFORM,bool>::type = true>
  inline std::pair<unsigned int, TIN> hash(TIN x) const {
    // nondimensionalized x position, scaled by step size
    TIN dx = m_stepSize_inv*(x-m_minArg);
    // index of table entry
    unsigned int x0 = static_cast<unsigned int>(dx);
    dx -= x0; // dx\in[0,1)
    return std::make_pair(x0, dx);
  }

  /* the polynomials for nonuniform LUTs map [x_k,x_{k+1}]->R so we don't have to change x */
  template <GridTypes GT1, typename std::enable_if<GT1 == GridTypes::NONUNIFORM,bool>::type = true>
  inline std::pair<unsigned int, TIN> hash(TIN x) const {
    // perform interval search in 6 FLOPs and 4 std::array<TIN,4> accesses
    unsigned int x0 = static_cast<unsigned int>(m_transferFunction.inverse(x));
    return std::make_pair(x0, x); // don't subtract dx by x0 because every polynomial was already rescaled accordingly
  }

  /* TODO Pade & LinearRawInterpTable must override this operator. The vtable MIGHT be optimized out,
   * but if each implementation provides their own operator() and diff() then we can be sure a
   * vtable isn't slowing a LUT down. */
  //#pragma omp declare simd // warning: GCC does not currently support mixed size types for 'simd' functions
  TOUT operator()(TIN x) const override {
    unsigned int x0; TIN dx;
    std::tie(x0,dx) = hash<GT>(x);
    
    // general degree horners method, evaluated from the inside out.
    TOUT sum = static_cast<TOUT>(0); // must initialize sum with a nonzero value to make LUT addition easy
    for(unsigned int k=N-1; k>0; k--)
      sum = (m_table[x0].coefs[k] + sum)*dx;
    return m_table[x0].coefs[0]+sum;
  }
  //TOUT diff(unsigned int N, TIN x) final {}
};



/* Reading & writing functions for any LUT derived from MetaTable.
 * Enables the convenient "get" syntax from nlohmann::json for specific implementations.
   eg:
  ```c++
  nlohmann::json jsonStats;
  std::ifstream(filename) >> jsonStats;
  auto lut = jsonStats.get<func::UniformLinearInterpolationTable<TIN,TOUT>>();
  ```
 * Uses SFINAE to automatically disable these functions if TIN or TOUT do not support nlohmann's json
 *
 * TODO can we just edit the json library to make the compile time errors into runtime errors? SFINAE takes forever to compile...
 * */
template <unsigned int N, typename TIN, typename TOUT, GridTypes GT,
         typename std::enable_if<std::is_constructible<nlohmann::json,TIN >::value && 
                                 std::is_constructible<nlohmann::json,TOUT>::value, bool>::type = true>
void to_json(nlohmann::json& jsonStats, const MetaTable<N,TIN,TOUT,GT>& lut) {
  jsonStats["_comment"] = "FunC lookup table data";
  jsonStats["name"] = lut.name();
  jsonStats["minArg"] = lut.min_arg();
  jsonStats["maxArg"] = lut.max_arg();
  jsonStats["order"] = lut.order();
  jsonStats["dataSize"] = lut.size();
  jsonStats["stepSize"] = lut.step_size();
  jsonStats["numTableEntries"] = lut.num_table_entries();
  jsonStats["numIntervals"] = lut.num_subintervals();
  jsonStats["tableMaxArg"] = lut.tablemax_arg();

  // things that are important for nonuniform tables:
  jsonStats["transfer_function_coefs"] = lut.transfer_function_coefs();

  // save the polynomial coefs of each lookup table
  // Note: m_order is used as the number of polynomial coefs
  for(unsigned int i=0; i<lut.num_table_entries(); i++)
    for(unsigned int j=0; j<lut.ncoefs_per_entry(); j++)
      jsonStats["table"][std::to_string(i)]["coefs"][std::to_string(j)] = lut.table_entry(i,j);
}

template <unsigned int N, typename TIN, typename TOUT, GridTypes GT,
         typename std::enable_if<!(std::is_constructible<nlohmann::json,TIN >::value && 
                                   std::is_constructible<nlohmann::json,TOUT>::value), bool>::type = true>
void to_json(nlohmann::json& jsonStats, const MetaTable<N,TIN,TOUT,GT>& lut) {
  (void) jsonStats;
  (void) lut;
  throw std::invalid_argument(std::string(typeid(TIN).name()) + " or " + std::string(typeid(TOUT).name()) + " does not have an implementation of nlohmann's to_json");
}

/* this variant of from_json will be called for any specific implementation of a LUT
 * inhereting from MetaTable */
template <unsigned int N, typename TIN, typename TOUT, GridTypes GT,
         typename std::enable_if<std::is_constructible<nlohmann::json,TIN >::value && 
                                 std::is_constructible<nlohmann::json,TOUT>::value, bool>::type = true>
void from_json(const nlohmann::json& jsonStats, MetaTable<N,TIN,TOUT,GT>& lut) {
  // name checking happens in MetaTable's constructor
  jsonStats.at("name").get_to(lut.m_name);
  jsonStats.at("minArg").get_to(lut.m_minArg);
  jsonStats.at("maxArg").get_to(lut.m_maxArg);
  jsonStats.at("stepSize").get_to(lut.m_stepSize);
  lut.m_stepSize_inv = 1.0/lut.m_stepSize;

  jsonStats.at("stepSize").get_to(lut.m_stepSize);
  jsonStats.at("order").get_to(lut.m_order);
  jsonStats.at("dataSize").get_to(lut.m_dataSize);
  jsonStats.at("numIntervals").get_to(lut.m_numIntervals);
  jsonStats.at("numTableEntries").get_to(lut.m_numTableEntries);
  jsonStats.at("tableMaxArg").get_to(lut.m_tableMaxArg);

  // Recompute m_table (the array of polynomials) and the transfer function
  lut.m_table.reset(new polynomial<TOUT,N>[lut.m_numTableEntries]);
  for(unsigned int i=0; i<lut.m_numTableEntries; i++)
    for(unsigned int j=0; j<lut.ncoefs_per_entry(); j++)
      jsonStats.at("table").at(std::to_string(i)).at("coefs").at(std::to_string(j)).get_to(lut.m_table[i].coefs[j]);

  // rebuild the transfer function
  lut.m_transferFunction = TransferFunction<TIN>(jsonStats["transfer_function_coefs"].get<std::array<TIN,4>>());
}

template <unsigned int N, typename TIN, typename TOUT, GridTypes GT,
         typename std::enable_if<!(std::is_constructible<nlohmann::json,TIN >::value && 
                                   std::is_constructible<nlohmann::json,TOUT>::value), bool>::type = true>
void from_json(const nlohmann::json& jsonStats, MetaTable<N,TIN,TOUT,GT>& lut) {
  (void) jsonStats;
  (void) lut;
  throw std::invalid_argument(std::string(typeid(TIN).name()) + " or " + std::string(typeid(TOUT).name()) + " does not implement nlohmann's to_json");
}

} // namespace func

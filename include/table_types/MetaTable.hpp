/*
   MetaTable handles any _piecewise polynomial based_ interpolation:
   - Reduce code redundancy by factoring out common differences between table types
   (eg method of generating a nonuniform grid type, setup/reading polynomial coefficients)
   which can be pieced together based on template parameters.

   NOTE: if stepSize divides max-min exactly then operator(max) will be out of array bounds!
   - Every table has an extra (unnecessary in most cases) table entry to avoid this
   problem.
   
   N = number of coefficients used in underlying piecewise polynomials
   Provided Horner's method which is the most common table evaluation method in FunC

   UNIFORM: Every subinterval is the same length so the hash is super fast; however,
    many more subintervals might be needed
   NONUNIFORM: Use a transfer function to create a nonuniform grid with an O(1) hash.
   NONUNIFORM_PSEUDO: same as NONUNIFORM but uses a faster, less accurate hash.

   TODO we could make another template which makes FunC save derivative coefs
   on the same cache line as the original coefs (hopefully with implementation
   independent of the LUT in use)
*/
#pragma once
#include "MetaTable.hpp"
#include "TransferFunctionSinh.hpp"


#include <array>
#include <stdexcept>
#include "json.hpp"

/* Give inheriting classes access to member variables without
   having to use "this->" excessively. These "using" statements must have protected access */
#define INHERIT_META(N,TIN,TOUT,GT) \
  using MetaTable<N,TIN,TOUT,GT>::m_func; \
  using MetaTable<N,TIN,TOUT,GT>::m_order; \
  using MetaTable<N,TIN,TOUT,GT>::m_name; \
  using MetaTable<N,TIN,TOUT,GT>::m_dataSize; \
  using MetaTable<N,TIN,TOUT,GT>::m_minArg; \
  using MetaTable<N,TIN,TOUT,GT>::m_maxArg \
  using MetaTable<N,TIN,TOUT,GT>::m_numIntervals; \
  using MetaTable<N,TIN,TOUT,GT>::m_numTableEntries; \
  using MetaTable<N,TIN,TOUT,GT>::m_stepSize; \
  using MetaTable<N,TIN,TOUT,GT>::m_stepSize_inv; \
  using MetaTable<N,TIN,TOUT,GT>::m_tableMaxArg \
  using MetaTable<N,TIN,TOUT,GT>::m_table; \
  using MetaTable<N,TIN,TOUT,GT>::m_grid; \
  using MetaTable<N,TIN,TOUT,GT>::m_transferFunction

/* Parallelization macro. Notes:
 * - m_table is aligned to sizeof(TOUT) so that might give some speedup.
 * */
#define FUNC_BUILDPAR _Pragma("omp parallel for schedule(dynamic)")
//#define FUNC_BUILDPAR


namespace func {

static constexpr unsigned int alignments[] = {0,1,2,4,4,8,8,8,8,16,16,16,16,16,16,16,16};

/* Our convention for writing polynomials is:
 *  p(x) = m_table[x0].coefs[0] + m_table[x0].coefs[1]*x + ... + m_table[x0].coefs[N-1]*x^(N-1)
 * operator() does a Horner's style evaluation>
 *
 * Note: LUTs can have other things in their polynomial struct
 * (2D LUTs may have coefs for x & y dimensions of each subsquare,
 * and LUTs could have that polynomial's derivative's coefs, etc). */
template <typename TOUT, unsigned int N>
struct alignas(sizeof(TOUT)*alignments[N]) polynomial {
  static const unsigned int ncoefs_per_interval = N;
  TOUT coefs[N];
};

enum class GridTypes {UNIFORM, NONUNIFORM, NONUNIFORM_PSEUDO};

template <GridTypes GT>
constexpr std::string grid_type_to_string() {
  switch(GT){
    case GridTypes::UNIFORM:
      return "Uniform";
    case GridTypes::NONUNIFORM:
      return "NonUniform";
    case GridTypes::NONUNIFORM_PSEUDO:
      return "NonUniformPseudo";
    default: { throw std::logic_error("Broken switch case in func::MetaTable"); }
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
  unsigned int m_dataSize;        // size of relevant data for impl evaluation

  unsigned int m_numIntervals;    // = (m_tableMaxArg - m_minArg)/m_stepSize;
  unsigned int m_numTableEntries; // length of m_grid and m_table (usually = m_numIntervals + 1)
  std::unique_ptr<TIN[]> m_grid;  // necessary for nonuniform tables
  __attribute__((aligned)) std::unique_ptr<polynomial<TOUT,N>[]> m_table; // holds polynomials coefficients
  TransferFunctionSinh<TIN> m_transferFunction; // used to make nonuniform grids (default constructable)

public:
  /* using a std::unique_ptr member variables implicitly deletes the default
   * move ctor so we must explicitly ask for the default move ctor */
  MetaTable(MetaTable&&) = default;
  MetaTable() = default;

  /* Set every generic member variable from a json file */
  MetaTable(const FunctionContainer<TIN,TOUT>& func_container, const LookupTableParameters<TIN>& par) :
    m_minArg(par.minArg), m_maxArg(par.maxArg), m_stepSize(par.stepSize)
  {
    /* If the step size does not exactly divide the arg domain, the max arg of the table is set
     * to the nearest value above such that it does. */
    if(m_stepSize <= static_cast<TIN>(0.0))
      throw std::invalid_argument("func::MetaTable was given a nonpositive stepSize. stepSize must be positive.");

    m_stepSize_inv = static_cast<TIN>(1.0)/m_stepSize;
    m_numIntervals = static_cast<unsigned>(ceil(m_stepSize_inv*(m_maxArg-m_minArg)));
    m_tableMaxArg = m_minArg+m_stepSize*numIntervals; // always >= m_maxArg

    // We need a valid FunctionContainer to generate any LUT
    if(m_func == nullptr)
      throw std::invalid_argument("Error in func::MetaTable. Function not defined in given FunctionContainer");

    // initialize the transfer function to something useful if the grid is nonuniform
    if(GT != GridTypes::UNIFORM)
      m_transferFunction = TransferFunctionSinh<TIN>(func_container,m_minArg,m_tableMaxArg,m_stepSize);
  }

  /* build this table from a json file */
  MetaTable(const nlohmann::json& jsonStats)
  {
    if(jsonStats.empty())
      throw std::invalid_argument("Error in func::MetaTable: The provided nlohmann::json is empty");
    from_json(jsonStats, *this);
  }

  void print_details_json(std::ostream& out) override {
    nlohmann::json jsonStats = *this; // call to_json(jsonStats, this)
    out << jsonStats.dump(2) << std::endl;
  }

  /* public access to protected member vars */
  TIN min_arg() const final { return m_minArg; }
  TIN max_arg() const final { return m_maxArg; }
  unsigned int order() const final { return m_order; }
  unsigned int size() const final { return m_dataSize; }
  std::string name() const final { return m_name; }
  TIN step_size() const final { return m_stepSize; };
  unsigned int num_table_entries() const final { return m_numTableEntries; };
  unsigned int num_intervals() const final { return m_numIntervals; };
  std::pair<TIN,TIN> bounds_of_subinterval(unsigned intervalNumber) const final
  {
    return std::make_pair(m_minArg + intervalNumber*m_stepSize,m_minArg + (intervalNumber+1)*m_stepSize);
  }

  unsigned int ncoefs_per_interval() const { return N; }
  TOUT table_entry(unsigned int i, unsigned int j) const { return m_table[i].coefs[j]; }
  TIN grid_entry(unsigned int i) const { return m_grid[i]; }
  std::array<TIN,4> transfer_function_coefs() const { return m_transferFunction.get_coefs(); }


  /* give the nlohmann json functions access to every member variable (note the silly template names)
   * TODO don't have to friend if to_json just uses the member functions */
  template <typename TIN1, typename TOUT1, unsigned int N1, GridTypes GT1,
         typename std::enable_if<std::is_constructible<nlohmann::json,TIN1 >::value && 
                                 std::is_constructible<nlohmann::json,TOUT1>::value, bool>::type>
  friend void to_json(nlohmann::json& jsonStats, const MetaTable<TIN1,TOUT1,N1,GT1>& lut);

  template <typename TIN1, typename TOUT1, unsigned int N1, GridTypes GT1,
         typename std::enable_if<std::is_constructible<nlohmann::json,TIN1 >::value && 
                                 std::is_constructible<nlohmann::json,TOUT1>::value, bool>::type>
  friend void from_json(const nlohmann::json& jsonStats, MetaTable<TIN1,TOUT1,N1,GT1>& lut);





  /* use SFINAE (COMPILES SLOWLY) to emulate partial template specialization
   * TODO a switch case might compile down to the exact same code. Profile this because switch case would compile way faster */
  template <typename TIN1, typename TOUT1, GridTypes GT1,
    typename std::enable_if<GT1 == GridTypes::UNIFORM,bool>::type = true>
  inline std::pair<unsigned int, TOUT1> hash(TIN1 x){
    unsigned int x0; TOUT1 dx;
    // nondimensionalized x position, scaled by step size
    dx  = static_cast<TOUT1>(m_stepSize_inv*(x-m_minArg));
    // index of previous table entry
    x0  = static_cast<unsigned>(dx);
    // value of table entries around x position
    dx -= x0;
    return std::make_pair(x0, dx);
  }

  template <typename TIN1, typename TOUT1, GridTypes GT1,
    typename std::enable_if<GT1 == GridTypes::NONUNIFORM,bool>::type = true>
  inline std::pair<unsigned int, TOUT1> hash(TIN1 x){
    unsigned int x0; TOUT1 dx;
    // find the subinterval x lives in
    x0 = static_cast<unsigned>(m_transferFunction.g_inv(x));
    // find where x is within that interval
    TIN1 h = m_grid[x0+1] - m_grid[x0];
    dx     = (x - m_grid[x0])/h;
    return std::make_pair(x0, dx);
  }

  template <typename TIN1, typename TOUT1, GridTypes GT1,
    typename std::enable_if<GT1 == GridTypes::NONUNIFORM_PSEUDO,bool>::type = true>
  inline std::pair<unsigned int, TOUT1> hash(TIN1 x){
    unsigned int x0; TOUT1 dx;
    // find the subinterval x lives in
    dx  = static_cast<TOUT1>(m_transferFunction.g_inv(x));
    // save some time by taking the fractional part of dx as x's location in this interval
    x0  = static_cast<unsigned>(dx);
    dx -= x0;
    return std::make_pair(x0, dx);
  }

  /* Implementations must provide their own operator() and diff() */
  /*
  #pragma omp declare simd
  TOUT operator()(TIN x) final
  {
    unsigned int x0; TOUT dx;
    std::tie(x0,dx) = hash(x);
    
    // general degree horners method, evaluated from the inside out.
    TOUT sum = 0;
    for(int k=N-1; k>0; k--)
      sum = dx*(m_table[x0].coefs[k] + sum);
    return m_table[x0].coefs[0]+sum;
  }
  TOUT diff(unsigned int N, TIN x) final {}
  */
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
template <typename TIN, typename TOUT, unsigned int N, GridTypes GT,
         typename std::enable_if<std::is_constructible<nlohmann::json,TIN >::value && 
                                 std::is_constructible<nlohmann::json,TOUT>::value, bool>::type = true>
void to_json(nlohmann::json& jsonStats, const MetaTable<TIN,TOUT,N,GT>& lut)
{
  jsonStats["_comment"] = "FunC lookup table data";
  jsonStats["name"] = lut.m_name;
  jsonStats["minArg"] = lut.m_minArg;
  jsonStats["maxArg"] = lut.m_maxArg;
  jsonStats["order"] = lut.m_order;
  jsonStats["dataSize"] = lut.m_dataSize;
  jsonStats["stepSize"] = lut.m_stepSize;
  jsonStats["numTableEntries"] = lut.m_numTableEntries;
  jsonStats["numIntervals"] = lut.m_numIntervals;
  jsonStats["tableMaxArg"] = lut.m_tableMaxArg;

  // things that are important for nonuniform tables:
  jsonStats["transfer_function_coefs"] = lut.transfer_function_coefs();
  for(unsigned int i=0; i<lut.m_numTableEntries; i++)
    jsonStats["grid"][std::to_string(i)] = lut.m_grid[i];

  // save the polynomial coefs of each lookup table
  // Note: m_order is used as the number of polynomial coefs
  for(unsigned int i=0; i<lut.m_numTableEntries; i++)
    for(unsigned int j=0; j<lut.ncoefs_per_interval(); j++)
      jsonStats["table"][std::to_string(i)]["coefs"][std::to_string(j)] = lut.table_entry(i,j);
}

template <typename TIN, typename TOUT, unsigned int N, GridTypes GT,
         typename std::enable_if<!(std::is_constructible<nlohmann::json,TIN >::value && 
                                   std::is_constructible<nlohmann::json,TOUT>::value), bool>::type = true>
void to_json(nlohmann::json& jsonStats, const MetaTable<TIN,TOUT,N,GT>& lut)
{
  (void) jsonStats;
  (void) lut;
  throw std::invalid_argument(std::string(typeid(TIN).name()) + " or " + std::string(typeid(TOUT).name()) + " does not implement nlohmann's to_json");
}

/* this variant of from_json will be called for any specific implementation of a LUT
 * inhereting from MetaTable */
template <typename TIN, typename TOUT, unsigned int N, GridTypes GT,
         typename std::enable_if<std::is_constructible<nlohmann::json,TIN >::value && 
                                 std::is_constructible<nlohmann::json,TOUT>::value, bool>::type = true>
void from_json(const nlohmann::json& jsonStats, MetaTable<TIN,TOUT,N,GT>& lut)
{
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

  // read array of points used
  lut.m_grid.reset(new TIN[lut.m_numTableEntries]);
  for(unsigned int i=0; i<lut.m_numTableEntries; i++)
    jsonStats.at("grid").at(std::to_string(i)).get_to(lut.m_grid[i]);

  // Recompute m_table (the array of polynomials) and the transfer function
  lut.m_table.reset(new polynomial<TOUT,N>[lut.m_numTableEntries]);
  for(unsigned int i=0; i<lut.m_numTableEntries; i++)
    for(unsigned int j=0; j<lut.ncoefs_per_interval(); j++)
      jsonStats.at("table").at(std::to_string(i)).at("coefs").at(std::to_string(j)).get_to(lut.m_table[i].coefs[j]);

  // rebuild the transfer function
  lut.m_transferFunction = TransferFunctionSinh<TIN>(jsonStats["transfer_function_coefs"].get<std::array<TIN,4>>());
}

template <typename TIN, typename TOUT, unsigned int N, GridTypes GT,
         typename std::enable_if<!(std::is_constructible<nlohmann::json,TIN >::value && 
                                   std::is_constructible<nlohmann::json,TOUT>::value), bool>::type = true>
void from_json(const nlohmann::json& jsonStats, MetaTable<TIN,TOUT,N,GT>& lut)
{
  (void) jsonStats;
  (void) lut;
  throw std::invalid_argument(std::string(typeid(TIN).name()) + " or " + std::string(typeid(TOUT).name()) + " does not implement nlohmann's to_json");
}

} // namespace func

#pragma once
#include "LookupTable.hpp"
#include "TransferFunction.hpp"
#include "Polynomial.hpp"

#include <array>
#include <stdexcept>
#include "json.hpp"

/* Classes inheriting MetaTable need this macro to access member variables without
   writing "this->" excessively. These "using" statements must have protected access */
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

enum class GridTypes {UNIFORM, NONUNIFORM};

template <GridTypes GT>
inline std::string grid_type_to_string() {
  FUNC_IF_CONSTEXPR(GT == GridTypes::UNIFORM){
    return "Uniform";
  }else{
    return "NonUniform";
  } 
}

/**
  \brief MetaTable handles any piecewise _polynomial_ based interpolation. Highly templated
  \ingroup MetaTable

  Note:
  - In the case where (max-min)/stepsize is not an integer then the actual
    maximum allowable argument `m_tableMaxArg` is `std::ceil(m_stepSize_inv*(m_maxArg-m_minArg)))`
  - If stepsize divides max-min exactly, then hash(max) is one too large. To
    resolve this issue, every implementation must have an extra (unnecessary in
    most cases) m_table entry
  - If TIN approximates a field and TOUT approximates a vector space over TIN,
    then MetaTable<N,TIN,TOUT,GT> approximates a vector space over TIN.
    MetaTable provides methods for addition/subtraction & scalar
    multiplication/division. This feature makes the N-D LUTs possible. 
  - Provides functions to read/write LUTs from a JSON file by implementing
    nlohmann's to/from_json.

  \tparam N is the number of coefficients required to represent each polynomial. So, N = deg(p_k) + 1.
  \tparam TIN is the type of the inputs passed to f.
  \tparam TOUT is the type that f outputs (and the type of the coefficients of each p_k).
  \tparam GT is the type of partition that this LUT uses. The options are
  - UNIFORM: Every subinterval is the same length so the hash takes 4 FLOPs &
    zero comparisons; however, unnecessary subintervals might be needed
  - NONUNIFORM: Use a TransferFunction to create a nonuniform partition of [a,b] with an O(1) hash.
*/
template <unsigned int N, typename TIN, typename TOUT=TIN, GridTypes GT=GridTypes::UNIFORM>
class MetaTable : public LookupTable<TIN,TOUT>
{
protected:
  std::string m_name;             //!< name of implementation type
  TIN m_minArg, m_maxArg;         //!< bounds of evaluation
  TIN m_stepSize, m_stepSize_inv; //!< fixed grid spacing (and its inverse)
  TIN m_tableMaxArg; //!< \> m_maxArg if (m_maxArg-m_minArg)/m_stepSize is non-integer

  unsigned int m_order;           //!< order of accuracy of implementation
  std::size_t  m_dataSize;        //!< size of relevant data for impl evaluation
  unsigned int m_numIntervals;    //!< = (m_tableMaxArg - m_minArg)/m_stepSize;
  unsigned int m_numTableEntries; //!< length of m_table (usually = m_numIntervals + 1)
  __attribute__((aligned)) std::unique_ptr<polynomial<TOUT,N>[]> m_table; //!< holds polynomials coefficients
  TransferFunction<TIN> m_transferFunction; //!< used to make nonuniform grids (default constructable)

  /** \brief Copy-swap pattern, necessary for `operator=`
   *  \post L is overwritten by the data in `this*` and `this*` loses access to its data. */
  void swap(MetaTable<N,TIN,TOUT,GT>& L) noexcept {
    L.m_name = m_name;
    L.m_minArg = m_minArg;
    L.m_maxArg = m_maxArg;
    L.m_order = m_order;
    L.m_dataSize = m_dataSize;
    L.m_numIntervals = m_numIntervals;
    L.m_numTableEntries = m_numTableEntries;
    L.m_stepSize = m_stepSize;
    L.m_stepSize_inv = m_stepSize_inv;
    L.m_tableMaxArg = m_tableMaxArg;
    L.m_transferFunction = m_transferFunction;
    L.m_table = std::move(m_table);
    return;
  }

public:
  /* using a std::unique_ptr member variables implicitly deletes the default
   * move ctor so we must explicitly ask for the default move ctor */
  MetaTable() = default;

  /** deepcopy constructor */
  MetaTable(const MetaTable<N,TIN,TOUT,GT>& L) : m_name(L.m_name), m_minArg(L.m_minArg), m_maxArg(L.m_maxArg),
    m_stepSize(L.m_stepSize), m_stepSize_inv(L.m_stepSize_inv), m_tableMaxArg(L.m_tableMaxArg),
    m_order(L.m_order), m_dataSize(L.m_dataSize), m_numIntervals(L.m_numIntervals), m_numTableEntries(L.m_numTableEntries)
  {
    m_table.reset(new polynomial<TOUT,N>[m_numTableEntries]);
    #pragma omp simd collapse(2)
    for(unsigned int ii=0; ii<m_numTableEntries; ++ii)
      for(unsigned int jj=0; jj<N; ++jj)
        m_table[ii].coefs[jj] = L.m_table[ii].coefs[jj];
  }

  /** Implemented with the copy-swap pattern */
  MetaTable<N,TIN,TOUT,GT>& operator=(MetaTable<N,TIN,TOUT,GT> L){
    L.swap(*this);
    return *this;
  }

  /** Use a json file to set every generic member variable */
  MetaTable(const FunctionContainer<TIN,TOUT>& func_container, const LookupTableParameters<TIN>& par, const nlohmann::json& jsonStats) :
    m_minArg(par.minArg), m_maxArg(par.maxArg), m_stepSize(par.stepSize) {
    /* build this table from a json file */
    if(!jsonStats.empty()){
      from_json(jsonStats, *this);
      return;
    }

    /* only positive stepsizes allowed */
    if(m_stepSize <= static_cast<TIN>(0.0))
      throw std::invalid_argument("func::MetaTable was given a nonpositive stepSize. stepSize must be positive.");

    /* If the step size does not exactly divide the arg domain, the max arg of the table is set
     * to the nearest value above such that it does. */
    m_stepSize_inv = static_cast<TIN>(1.0)/m_stepSize;
    using std::ceil;
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

  /* public access to protected member vars */
  std::string name() const final { return m_name; }
  TIN min_arg() const final { return m_minArg; }
  TIN max_arg() const final { return m_maxArg; }
  TIN tablemax_arg() const { return m_tableMaxArg; }
  unsigned int order() const final { return m_order; }
  std::size_t size() const final { return m_dataSize; }
  unsigned int num_subintervals() const final { return m_numIntervals; }
  TIN step_size() const final { return m_stepSize; }
  std::pair<TIN,TIN> bounds_of_subinterval(unsigned int intervalNumber) const final {
    // initialize to bounds of uniform LUT and adjust if LUT is nonuniform
    TIN intervalMin = m_minArg + intervalNumber*m_stepSize;
    TIN intervalMax = m_minArg + (intervalNumber+1)*m_stepSize;
    FUNC_IF_CONSTEXPR(GT != GridTypes::UNIFORM){
      intervalMin = m_transferFunction(intervalMin);
      intervalMax = m_transferFunction(intervalMax);
    }
    // m_tableMaxArg can be greater than m_maxArg
    if(intervalMax > m_maxArg) intervalMax = m_maxArg;
    return std::make_pair(intervalMin,intervalMax);
  }
  void print_json(std::ostream& out) const final {
    nlohmann::json jsonStats = *this; // call to_json(jsonStats, this)
    out << jsonStats.dump(2) << std::endl;
  }

  unsigned int num_table_entries() const { return m_numTableEntries; }
  unsigned int ncoefs_per_entry() const { return N; }
  TOUT table_entry(unsigned int i, unsigned int j) const { return m_table[i].coefs[j]; }
  std::array<TIN,4> transfer_function_coefs() const { return m_transferFunction.get_coefs(); }


  /* give the nlohmann json functions access to every member variable (must use SFINAE because we want people
   * to construct LUTs over arbitrary types regardless of whether they have to/from json functions) */
  template <unsigned int N1, typename TIN1, typename TOUT1, GridTypes GT1,
         typename std::enable_if<std::is_constructible<nlohmann::json,TIN1 >::value && 
                                 std::is_constructible<nlohmann::json,TOUT1>::value, bool>::type>
  friend void from_json(const nlohmann::json& jsonStats, MetaTable<N1,TIN1,TOUT1,GT1>& lut);

  /* LUTs form a vector space over TIN if TOUT forms a vector space over TIN
   * Note these operations could be heavily optimized with template expressions
   * but that doesn't work with the auto keyword (see operator()(...)) and that'd be a lot of work... 
   *
   * Maybe we could implement polynomial using an existing linear algebra library with fast
   * (possibly template expression based) operator+ and operator*? */
  MetaTable<N,TIN,TOUT,GT>& operator+=(const MetaTable<N,TIN,TOUT,GT>& other) {
    if((m_numTableEntries != other.m_numTableEntries) || (m_minArg != other.m_minArg) || (m_maxArg != other.m_maxArg))
      throw std::invalid_argument("Error in func::MetaTable: cannot add two LUTs with different subintervals");

    /* TODO Transfer functions don't work with LUTs of LUTs yet I think...
     * Try building nonuniform LUTs of LUTs by lerping transfer functions? */

    #pragma omp simd collapse(2)
    for(unsigned int ii=0; ii<m_numTableEntries; ++ii)
      for(unsigned int jj=0; jj<N; ++jj)
        m_table[ii].coefs[jj] += other.m_table[ii].coefs[jj];
    return *this;
  }
  MetaTable<N,TIN,TOUT,GT>& operator-=(const MetaTable<N,TIN,TOUT,GT>& other) {
    if((m_numTableEntries != other.m_numTableEntries) || (m_minArg != other.m_minArg) || (m_maxArg != other.m_maxArg))
      throw std::invalid_argument("Error in func::MetaTable: cannot subtract two LUTs with different subintervals");

    #pragma omp simd collapse(2)
    for(unsigned int ii=0; ii<m_numTableEntries; ++ii)
      for(unsigned int jj=0; jj<N; ++jj)
        m_table[ii].coefs[jj] -= other.m_table[ii].coefs[jj];
    return *this;
  }
  MetaTable<N,TIN,TOUT,GT>& operator*=(const TIN& a) {
    #pragma omp simd collapse(2)
    for(unsigned int ii=0; ii<m_numTableEntries; ++ii)
      for(unsigned int jj=0; jj<N; ++jj)
        m_table[ii].coefs[jj] *= a;
    return *this;
  }
  MetaTable<N,TIN,TOUT,GT>& operator/=(const TIN& a) {
    #pragma omp simd collapse(2)
    for(unsigned int ii=0; ii<m_numTableEntries; ++ii)
      for(unsigned int jj=0; jj<N; ++jj)
        m_table[ii].coefs[jj] /= a;
    return *this;
  }


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

  /* the polynomials for nonuniform LUTs map [x_k,x_{k+1}]->R so we don't have to change x.
   * Calls m_transferFunction.inverse(x), using 6 FLOPs and 4 std::array<TIN,4> access */
  template <GridTypes GT1, typename std::enable_if<GT1 == GridTypes::NONUNIFORM,bool>::type = true>
  inline std::pair<unsigned int, TIN> hash(TIN x) const {
    unsigned int x0 = static_cast<unsigned int>(m_transferFunction.inverse(x));
    return std::make_pair(x0, x); // we don't subtract dx by x0 because every polynomial was rescaled during construction
  }

  /** Find the subinterval [x_k,x_{k+1}) that x belongs to, fetch polynomial
   * coefficients from m_table[k], and use Horner's method to compute p_k(x).
   *
   * TODO Pade & LinearRawInterpTable must override this operator. Maybe
   * operator() will be faster if each implementation provides their own
   * operator() and diff(). If the vtable isn't optimized out then perchance
   * removing the use of virtual will remove that overhead. */
  //#pragma omp declare simd // warning: GCC does not currently support mixed size types for 'simd' functions
  TOUT operator()(TIN x) const override {
    unsigned int x0; TIN dx;
    std::tie(x0,dx) = hash<GT>(x);
    
    // general degree horners method, evaluated from the inside out.
    TOUT sum = m_table[x0].coefs[N-1];
    for(unsigned int k=N-1; k>0; k--){
      sum *= dx;
      sum += m_table[x0].coefs[k-1];
    }
    return sum;
  }

  /** \brief If LUT coefficients (of type TOUT) are callable, then apply each to `args` as soon as possible.
   *
   * Notes:
   * - We must use an `auto` return type because there's no straightforward way
   *   to statically deduce the return type ahead of time.
   * - This function does not perform well with template expressions (e.g. LUTs
   *   of arma::mat) so use the other operator() in that case.
   * - This variadic operator() cannot be in the LUT interface because it is
   *   templated. */
  template<typename... TIN2>
  inline auto operator()(TIN x, TIN2 ... args) const {
    unsigned int x0; TIN dx;
    std::tie(x0,dx) = hash<GT>(x);

    auto sum = m_table[x0].coefs[N-1](args...);
    for(unsigned int k=N-1; k>0; k--){
      sum *= dx;
      sum += m_table[x0].coefs[k-1](args...);
    }
    return sum;
  }

  /** \brief Return the sth derivative of L at x: p_k^(s)(x)
   * TODO make this function virtual and override in Pade and linear raw tables */
  TOUT diff(unsigned int s, TIN x) const {
    unsigned int x0; TIN dx;
    std::tie(x0,dx) = hash<GT>(x);

    auto sum = static_cast<TIN>(permutation(N-1,s))*m_table[x0].coefs[N-1];
    for(unsigned int k=N-1; k>s; k--){
      sum *= dx;
      sum += static_cast<TIN>(permutation(k-1,s))*m_table[x0].coefs[k-1];
    }
    return static_cast<TIN>(pow(m_stepSize_inv,s))*sum;
  }

  /** \brief Same pattern as the variadic operator() but for the s^th derivative. */
  template<typename... TIN2>
  inline auto diff(unsigned int s, TIN x, TIN2... args) const {
    unsigned int x0; TIN dx;
    std::tie(x0,dx) = hash<GT>(x);

    auto sum = static_cast<TIN>(permutation(N-1,s))*(m_table[x0].coefs[N-1].diff(args...));
    for(unsigned int k=N-1; k>s; k--){
      sum *= dx;
      sum += static_cast<TIN>(permutation(k-1,s))*(m_table[x0].coefs[k-1].diff(args...));
    }
    return static_cast<TIN>(pow(m_stepSize_inv,s))*sum;
  }
};

template <unsigned int N, typename TIN, typename TOUT, GridTypes GT>
MetaTable<N,TIN,TOUT,GT> operator+(MetaTable<N,TIN,TOUT,GT> lhs, const MetaTable<N,TIN,TOUT,GT>& rhs){
    return lhs += rhs;
}
template <unsigned int N, typename TIN, typename TOUT, GridTypes GT>
MetaTable<N,TIN,TOUT,GT> operator-(MetaTable<N,TIN,TOUT,GT> lhs, const MetaTable<N,TIN,TOUT,GT>& rhs){
    return lhs -= rhs;
}
template <unsigned int N, typename TIN, typename TOUT, GridTypes GT>
MetaTable<N,TIN,TOUT,GT> operator*(MetaTable<N,TIN,TOUT,GT> lhs, const TIN& scalar){
    return lhs *= scalar;
}
template <unsigned int N, typename TIN, typename TOUT, GridTypes GT>
MetaTable<N,TIN,TOUT,GT> operator*(const TIN& scalar, MetaTable<N,TIN,TOUT,GT> lhs){
    return lhs *= scalar;
}
template <unsigned int N, typename TIN, typename TOUT, GridTypes GT>
MetaTable<N,TIN,TOUT,GT> operator/(MetaTable<N,TIN,TOUT,GT> lhs, const TIN& scalar){
    return lhs /= scalar;
}

/** \brief Reading & writing functions for any LUT derived from MetaTable.
 * Enables the convenient "get" syntax from nlohmann::json for any class that implements MetaTable.
 *
   eg:
  \code{.cpp}
  nlohmann::json jsonStats;
  std::ifstream(filename) >> jsonStats; // call to_json
  auto lut = jsonStats.get<func::UniformLinearInterpolationTable<TIN,TOUT>>(); // call from_json
  \endcode
 * Uses SFINAE to automatically disable these functions if TIN or TOUT do not
 * support nlohmann's json (then users can have LUTs over abstract types
 * without having to implement to/from_json first). SFINAE is also necessary
 * because nlohmann uses static_asserts to check if TIN/TOUT have to/from_json
 *
 * TODO This use of SFINAE does (as expected) increase the compile
 * time by several percentage points. Can we edit the json library instead to
 * make the compile errors into runtime errors?
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

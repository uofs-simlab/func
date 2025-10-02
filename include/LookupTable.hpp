#pragma once
#include <string>
#include <stdexcept>
#include <functional>
#include <ostream> // operator<<
#include "FunctionContainer.hpp"

namespace func {

/**
 * \brief A struct containing data necessary/useful for constructing a LUT
 */
template <typename TIN, typename TOUT = TIN>
struct LookupTableParameters
{
  TIN minArg;
  TIN maxArg;
  TIN stepSize;
  /// special_points (roots, critical points, inflection points)
  std::vector<std::tuple<TIN,unsigned,TOUT>> special_points;

  LookupTableParameters(TIN min, TIN max, TIN step, std::vector<std::tuple<TIN,unsigned,TOUT>> pts) :
    minArg(min), maxArg(max), stepSize(step), special_points(pts) {}
  LookupTableParameters(TIN min, TIN max, TIN step) :
    LookupTableParameters(min, max, step, {}) {}
  LookupTableParameters(){}
};

/**
  \brief Abstract interface for representing an approximation to a user provided mathematical function.

  LookupTable possesses no member variables, or runnable code.
  Implementations of this class handle actual data (reading, writing, hashing, etc).
  The LookupTable interface is necessary for our LookupTableFactory,
  LookupTableGenerator, and LookupTableComparator.

  \warning { We make no promises about checking array bounds (as this notably reduces performance) }
*/
template <typename TIN, typename TOUT = TIN>
class LookupTable
{
public:
  using input_type = TIN;
  using output_type = TOUT;

  /** Every implementation of LookupTable will have a constructor that looks like this:
   * LookupTable(const FunctionContainer<TIN,TOUT>& func_container, const LookupTableParameters<TIN>& par) {} */
  LookupTable() = default;
  virtual ~LookupTable(){};

  virtual TOUT operator()(TIN x) const = 0;
  //virtual TOUT diff(unsigned int N, TIN x) = 0;

  /* public interface for access to protected data */
  virtual std::string name() const = 0;
  virtual TIN min_arg() const = 0;
  virtual TIN max_arg() const = 0;
  virtual unsigned int order() const = 0;
  virtual std::size_t size() const = 0;
  virtual unsigned int num_subintervals() const = 0;
  virtual TIN step_size() const = 0;
  virtual std::pair<TIN,TIN> bounds_of_subinterval(unsigned int intervalNumber) const = 0;

  /** \note Every class implementing LookupTable should call their implementation
   * of nlohmann's to_json from print_json */
  virtual void print_json(std::ostream& out) const = 0;
};

/** \brief Print basic info about a LookupTable if the user attempts to write L to std::out */
template <typename TIN, typename TOUT = TIN>
std::ostream& operator<<(std::ostream& out, const LookupTable<TIN,TOUT>& L){
  out << L.name() << " " << L.min_arg() << " " << L.max_arg() << " "
  << L.step_size() << " " << L.num_subintervals() << " ";
  return out;
}

/** \brief Return the output from operator<< as a string */
template <typename TIN, typename TOUT = TIN>
inline std::string to_string(const LookupTable<TIN,TOUT>& L) {
  std::ostringstream ss;
  ss << L;
  return ss.str();
}

} // namespace func

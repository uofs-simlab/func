/*
  This class is an interface that is useful for our LookupTableFactory.

  Actual data (reading, writing, hashing, etc) is handled
  by specific implementations of this class.

  Notes:
  - In the case where (max-min)/stepsize is not an integer then
  the actual maximum allowable argument is greater than max
*/
/*
  Abstract base class for any FunC object approximating a
  user provided function. This lets us make a vector with
  both a DirectEvaluation and various other LookupTables
  in ImplementationComparator.hpp.

  Derived classes generally do the following:
  - set up anything needed for 'evaluation'
  - determine size of data needed for 'evaluation'
  - override the brackets operator to perform 'evaluation'
  - write data to json with print_json
  - cleanup in destructor
*/

#pragma once
#include <string>
#include <stdexcept>
#include <functional>
#include <ostream> // operator<<
#include "FunctionContainer.hpp"

namespace func {

template <typename TIN>
struct LookupTableParameters
{
  TIN minArg;
  TIN maxArg;
  TIN stepSize;

  LookupTableParameters(TIN min, TIN max, TIN step) :
    minArg(min), maxArg(max), stepSize(step) {}
  LookupTableParameters(){}
};

template <typename TIN, typename TOUT = TIN>
class LookupTable
{
public:
  /* any implementation will have a constructor that looks like this:
   * LookupTable(const FunctionContainer<TIN,TOUT>& func_container, const LookupTableParameters<TIN>& par) {} */
  LookupTable() = default;
  virtual ~LookupTable(){};

  virtual TOUT operator()(TIN x) const = 0;
  //virtual TOUT diff(unsigned int N, TIN x) = 0;

  /* public interface for access to protected data */
  //virtual std::function<TOUT(TIN)> function() const = 0;
  virtual std::string name() const = 0;

  virtual TIN min_arg() const = 0;
  virtual TIN max_arg() const = 0;
  virtual unsigned int order() const = 0;
  virtual unsigned int size() const = 0;
  virtual unsigned int num_subintervals() const = 0;
  virtual TIN step_size() const = 0;
  virtual std::pair<TIN,TIN> bounds_of_subinterval(unsigned int intervalNumber) const = 0;

  /* any implementation of print_json should call their implementation of nlohmann's to_json */
  virtual void print_json(std::ostream& out) const = 0;
};

/* print basic info about a LookupTable */
template <typename TIN, typename TOUT = TIN>
std::ostream& operator<<(std::ostream& out, const LookupTable<TIN,TOUT>& L){
  out << L.name() << " " << L.min_arg() << " " << L.max_arg() << " "
  << L.step_size() << " " << L.num_subintervals() << " ";
  return out;
}

/* wraps operator<< */
template <typename TIN, typename TOUT = TIN>
inline std::string to_string(const LookupTable<TIN,TOUT>& L) {
  std::ostringstream ss;
  ss << L;
  return ss.str();
}

} // namespace func

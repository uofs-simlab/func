/*
  std::function wrapper with optional support for plotting the usage
  of a function's domain in a histogram.
  Replace your original function's usage with this class and compile
  with -DFUNC_DEBUG in order to figure out what bounds you should
  set your table intervals to.

  Usage example:
    DirectEvaluation de(&function,0,10);
    double f = de(0.87354);
    // sim code
    de.print_details(std::cout); // prints max/min recorded args if FUNC_DEBUG is defined

  Notes:
  - Wraps a std::function
  - When compiled with -DFUNC_DEBUG, min and max args are used for histogram's bounds
  which record arguments passed to the direct evaluation. No error if out of bounds.
  View the histogram with print_details(). Histogram code is in ArgumentRecord.hpp
  */
#pragma once

#include "FunctionContainer.hpp"
#include "config.hpp" // FUNC_USE_BOOST
#include "json.hpp"
#include "StdRng.hpp"

#include <fstream> //ifstream

#ifdef FUNC_DEBUG
  #include "ArgumentRecord.hpp"
#endif

namespace func {

template <typename TIN, typename TOUT = TIN>
class DirectEvaluation final : public LookupTable<TIN,TOUT>
{
  static constexpr const char * classname = "DirectEvaluation";
  std::function<TOUT(TIN)> m_func;
#ifdef FUNC_DEBUG
  std::unique_ptr<mutable ArgumentRecord<TIN>> mp_recorder;
  StdRng<TIN> m_sampler{0,1}; // uniformly distrubuted random numbers in [0,1]
  // TODO save/load error to json?
  TIN  m_rerr;
  TOUT m_aerr;
#endif
public:

  /* Setup argument recording if FUNC_DEBUG is defined used at compile time */
  DirectEvaluation(const FunctionContainer<TIN,TOUT>& func_container,
      TIN min = 0.0, TIN max = 1.0, unsigned int histSize = 10, TOUT aerr = 0.0, TIN rerr = 0.0) :
    m_func(func_container.standard_fun)
  {
    if(m_func == nullptr)
      throw std::invalid_argument("Error in func::DirectEvaluation: given a FunctionContainer with null function.");

    #ifdef FUNC_DEBUG
      mp_recorder = std::unique_ptr<ArgumentRecord<TIN>>(new ArgumentRecord<TIN>(min, max, histSize));
      m_rerr = rerr;
      m_aerr = aerr;
    #endif
    /* ignore debugging args if -DFUNC_DEBUG isn't specified */
    (void) min; (void) max; (void) histSize; (void) rerr; (void) aerr;
  }

  /* rebuild this class and it's arg record from a file */
  //DirectEvaluation(const FunctionContainer<TIN,TOUT>& func_container, std::string filename) :
  //  m_func(func_container.standard_fun)
  //{
  //  nlohmann::json jsonStats;
  //  std::ifstream(filename) >> jsonStats;
  //#ifdef FUNC_DEBUG
  //  /* reconstruct the arg record */
  //  mp_recorder = std::make_unique<ArgumentRecord<TIN>>(jsonStats);
  //#endif
  //}

  /* Evaluate the underlying std::function and optionally record the arg */
  TOUT operator()(TIN x) const final
  {
    #ifdef FUNC_DEBUG
      mp_recorder->record_arg(x);
      return m_aerr*m_sampler.get_point() + m_func(x)*(static_cast<TIN>(1.0) + m_rerr*m_sampler.get_point());
    #endif
    return m_func(x);
  }

  std::string name() const final { return classname; }
  TIN min_arg() const final { return -std::numeric_limits<TIN>::infinity(); };
  TIN max_arg() const final { return std::numeric_limits<TIN>::infinity(); };
  unsigned int order() const final { return std::numeric_limits<unsigned int>::infinity(); };
  unsigned int size() const final { return 0u; };
  unsigned int num_subintervals() const final { return 0u; };
  TIN step_size() const final { return static_cast<TIN>(0); };
  std::pair<TIN,TIN> bounds_of_subinterval(unsigned int intervalNumber) const final { (void) intervalNumber; return std::make_pair(min_arg(),max_arg()); };

  /* TODO make to_json() */
  void print_json(std::ostream& out) const final {
    out << name() << "\n";
    #ifdef FUNC_DEBUG
    //out << mp_recorder->print_json(out);
    #endif
  }

  ~DirectEvaluation(){};
};

template <typename TIN, typename TOUT = TIN>
std::ostream& operator<<(std::ostream& out, const DirectEvaluation<TIN,TOUT>& D){
  out << D.name() << "\n"; // TODO call superclass operator<<
#ifdef FUNC_DEBUG
  out << *D.mp_recorder;
#endif
  return out;
}

// TODO make to_json()
} // namespace func

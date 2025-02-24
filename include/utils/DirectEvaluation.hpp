#pragma once

#include "FunctionContainer.hpp"
#include "config.hpp" // FUNC_USE_BOOST
#include "json.hpp"
#include "StdRng.hpp"

#include <fstream> // ifstream, ostream

#ifdef FUNC_DEBUG
  #include "ArgumentRecord.hpp"
#endif

namespace func {

/**
  \brief Wrap a std::function and optionally plot that function's domain usage
   with an ArgumentRecord (builds a histogram). To determine useful LUT bounds,
   users should replace their mathematical function with this class and compile
   with -DFUNC_DEBUG.

  \ingroup Utils

  Usage example:
  \code{.cpp}
  DirectEvaluation<double> de({MyFunction},0,10);
  double fx = de(0.87354);
  // sim code calling de goes here
  de.print_details(std::cout); // prints max/min recorded args if FUNC_DEBUG is defined
  \endcode

  \notes When compiled with -DFUNC_DEBUG, the ArgumentRecord uses min and max
  constructor arguments to construct a histogram's bounds. This histogram
  record arguments passed to the DirectEvaluation and there is no issue if
  sampled arguments are out of bounds.
  \note View the histogram with print_details(), or construct DirectEvaluation
  with a pointer to ostream and get output upon destruction.
  */
template <typename TIN, typename TOUT = TIN>
class DirectEvaluation final : public LookupTable<TIN,TOUT>
{
  std::function<TOUT(TIN)> m_func;
#ifdef FUNC_DEBUG
  mutable std::unique_ptr<ArgumentRecord<TIN>> mp_recorder;
  mutable StdRng<TIN> m_sampler{0,1}; // uniformly distrubuted random numbers in [0,1]
  // TODO save/load error to json?
  TIN  m_rerr;
  TOUT m_aerr;
#endif
public:

  /* Setup argument recording if FUNC_DEBUG is defined used at compile time */
  DirectEvaluation(const FunctionContainer<TIN,TOUT>& func_container,
      TIN min = 0.0, TIN max = 1.0, unsigned int histSize = 10, TOUT aerr = 0.0, TIN rerr = 0.0, std::ostream* streamer = nullptr) :
    m_func(func_container.standard_fun)
  {
    if(m_func == nullptr)
      throw std::invalid_argument("Error in func::DirectEvaluation: given a FunctionContainer with null function.");

    #ifdef FUNC_DEBUG
      mp_recorder = std::unique_ptr<ArgumentRecord<TIN>>(new ArgumentRecord<TIN>(min, max, histSize,streamer));
      m_rerr = rerr;
      m_aerr = aerr;
    #endif
    /* ignore debugging args if -DFUNC_DEBUG isn't specified */
    (void) min; (void) max; (void) histSize; (void) rerr; (void) aerr; (void) streamer;
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

  /** \brief Evaluate the underlying std::function and optionally record the argument x */
  TOUT operator()(TIN x) const final
  {
    #ifdef FUNC_DEBUG
      mp_recorder->record_arg(x);
      return m_aerr*m_sampler.get_point() + m_func(x)*(static_cast<TIN>(1.0) + m_rerr*m_sampler.get_point());
    #endif
    return m_func(x);
  }

  std::string name() const final { return "DirectEvaluation"; }
  TIN min_arg() const final { return -std::numeric_limits<TIN>::infinity(); };
  TIN max_arg() const final { return std::numeric_limits<TIN>::infinity(); };
  unsigned int order() const final { return std::numeric_limits<unsigned int>::infinity(); };
  std::size_t size() const final { return 0u; };
  unsigned int num_subintervals() const final { return 0u; };
  TIN step_size() const final { return static_cast<TIN>(0); };
  std::pair<TIN,TIN> bounds_of_subinterval(unsigned int intervalNumber) const final { (void) intervalNumber; return std::make_pair(min_arg(),max_arg()); };

  void print_json(std::ostream& out) const final {
    (void) out;
    #ifdef FUNC_DEBUG
    out << mp_recorder->print_json(out);
    #endif
  }

  ~DirectEvaluation(){};
};

template <typename TIN, typename TOUT = TIN>
std::ostream& operator<<(std::ostream& out, const DirectEvaluation<TIN,TOUT>& D){
  out << D.name() << "\n"; // TODO call superclass operator<<
  D.print_hist(out);
  return out;
}

// TODO make to_json()
} // namespace func

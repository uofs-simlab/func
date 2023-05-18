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
  - Used to store a std::function as an EvaluationImplementation which
  made ImplementationComparator easier to implement and use.
  - requires min and max args for consistency with the EvaluationImplementation.
  Once compiled with -DFUNC_DEBUG, they're used as a rough guess for the
  histogram's bounds. print_details() will tell you if they were a bad guess.
  - Record where the function is being evaluated by specifying the -DFUNC_DEBUG
  flag at compile time. ArgumentRecorder is an extension of this class which does
  all the work recording function arguments.
*/
#pragma once

#include "EvaluationImplementation.hpp"
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
#ifdef FUNC_DEBUG
  std::unique_ptr<ArgumentRecord<TIN>> mp_recorder;
  StdRng<TIN> m_sampler{0,1}; // uniformly distrubuted random numbers in [0,1]
  // TODO save/load error to json?
  TIN  m_rerr;
  TOUT m_aerr;
#endif
public:

  /* Setup argument recording if FUNC_DEBUG is defined used at compile time */
  DirectEvaluation(const FunctionContainer<TIN,TOUT>& func_container,
      TIN min = 0, TIN max = 1, unsigned int histSize = 10, TIN aerr = 0, TIN rerr = 0) :
    m_name("DirectEvaluation")
  {
    #ifdef FUNC_DEBUG
      mp_recorder = std::unique_ptr<ArgumentRecord<TIN>>(new ArgumentRecord<TIN>(min, max, histSize));
      m_rerr = rerr;
      m_aerr = aerr;
    #endif
    /* ignore debugging args if -DFUNC_DEBUG isn't specified */
    (void) histSize; (void) rerr; (void) aerr;
  }

  /* rebuild this class and it's arg record from a file */
  DirectEvaluation(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    EvaluationImplementation<TIN,TOUT>(func_container->standard_func, "DirectEvaluation")
  {
    nlohmann::json jsonStats;
    std::ifstream(filename) >> jsonStats;
    m_name = jsonStats["name"].get<std::string>();

  #ifdef FUNC_DEBUG
    // reconstruct our arg record
    mp_recorder = std::unique_ptr<ArgumentRecord<TIN>>(new ArgumentRecord<TIN>(jsonStats));
  #endif
  }

  // Evaluate the underlying std::function and optionally record the arg
  TOUT operator()(TIN x) final
  {
    #ifdef FUNC_DEBUG
      mp_recorder->record_arg(x);
      return m_aerr*m_sampler.get_point() + m_func(x)*(static_cast<TIN>(1.0) + m_rerr*m_sampler.get_point());
    #endif
    return m_func(x);
  }

  //void print_details(std::ostream& out) override;
  //void print_details_json(std::ostream& out) override;
  ~DirectEvaluation(){};
};

/* Implementations of functions that print to the console */
template <typename TIN, typename TOUT>
inline void DirectEvaluation<TIN,TOUT>::print_details(std::ostream& out)
{
  out<< m_name;
  #ifdef FUNC_DEBUG
    out << std::endl;
    mp_recorder->print_details(out);
  #endif
}

template <typename TIN, typename TOUT>
inline void DirectEvaluation<TIN,TOUT>::print_details_json(std::ostream& out)
{
  nlohmann::json jsonStats;

  jsonStats["_comment"] = "FunC DirectEvaluation data";
  jsonStats["name"] = m_name;
  jsonStats["minArg"] = m_minArg;
  jsonStats["maxArg"] = m_maxArg;
  #ifdef FUNC_DEBUG
    // have our ArgumentRecord add it's own data
    mp_recorder->print_details_json(jsonStats);
  #endif

  out << jsonStats.dump(2) << std::endl;
}
// TODO make to_json()
} // namespace func

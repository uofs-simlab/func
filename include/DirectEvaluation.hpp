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
class DirectEvaluation final : public EvaluationImplementation<TIN,TOUT>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  #ifdef FUNC_DEBUG
    std::unique_ptr<ArgumentRecord<TIN>> mp_recorder;
    StdRng<TOUT> mp_sampler{-1,1};
    TOUT m_tol; // TODO save/load to file?
  #endif
public:

  /* set base class args and set up argument recording if 
     FUNC_DEBUG is defined used at compile time */
  DirectEvaluation(FunctionContainer<TIN,TOUT> *func_container,
      TIN min = 0, TIN max = 1, unsigned int histSize = 10, TOUT tol = 0) :
    EvaluationImplementation<TIN,TOUT>(func_container->standard_func, "DirectEvaluation")
  {
    m_minArg = min, m_maxArg = max, m_dataSize = 0;
    #ifdef FUNC_DEBUG
      mp_recorder = std::unique_ptr<ArgumentRecord<TIN>>(new ArgumentRecord<TIN>(min, max, histSize));
      m_tol = tol;
    #endif
    (void) histSize; // ignore debugging args if -DFUNC_DEBUG isn't specified
    (void) tol;
  }

  /* rebuild this class and it's arg record from a file */
  DirectEvaluation(FunctionContainer<TIN,TOUT> *func_container, std::string filename) :
    EvaluationImplementation<TIN,TOUT>(func_container->standard_func, "DirectEvaluation")
  {
    nlohmann::json jsonStats;
    std::ifstream(filename) >> jsonStats;
    m_name = jsonStats["name"].get<std::string>();
    m_minArg = jsonStats["minArg"].get<TIN>();
    m_maxArg = jsonStats["maxArg"].get<TIN>();

  #ifdef FUNC_DEBUG
    // reconstruct our arg record
    mp_recorder = std::unique_ptr<ArgumentRecord<TIN>>(new ArgumentRecord<TIN>(jsonStats));
  #endif
  }

  // Evaluate the underlying std::function and optionally record the arg
  TOUT operator()(TIN x) override
  {
    #ifdef FUNC_DEBUG
      mp_recorder->record_arg(x);
      return m_func(x)*(((TOUT) 1) + m_tol*mp_sampler.get_point());
    #endif
    return m_func(x);
  }

  void print_details(std::ostream& out) override;
  void print_details_json(std::ostream& out) override;
  ~DirectEvaluation(){};
};

/* Implementations of functions that print to the console */
template <typename TIN, typename TOUT>
inline void DirectEvaluation<TIN,TOUT>::print_details(std::ostream& out)
{
  out<< m_name << " " << m_minArg << " " << m_maxArg;
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

/*
  std::function wrapper with optional support
  for plotting the usage of a function's domain

  Usage example:
    DirectEvaluation de(&function,0,10);
    double f = de(0.87354);

  Notes:
  - basically just a wrapper around a function
  - requires a max and min value (for consistency with EvaluationImplementation)
  - Record where the function is being evaluated by specifying the -DFUNC_RECORD
  flag at compile time. ArgumentRecorder is an extension of this class which does
  all the work recording function arguments.
*/
#pragma once

#include "EvaluationImplementation.hpp"
#include "FunctionContainer.hpp"
#include "json.hpp"
#include <cmath>

#ifdef FUNC_RECORD
  #include "ArgumentRecord.hpp"
#endif

template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE>
class DirectEvaluation final : public EvaluationImplementation<IN_TYPE,OUT_TYPE>
{
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  #ifdef FUNC_RECORD
    std::unique_ptr<ArgumentRecord> mp_recorder;
  #endif
public:

  /* set base class args and set up argument recording if 
     -DFUNC_RECORD is used at compile time */
  DirectEvaluation(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container,
      IN_TYPE min = 0, IN_TYPE max = 1, unsigned int histSize = 10) :
    EvaluationImplementation<IN_TYPE,OUT_TYPE>(func_container->standard_func, "DirectEvaluation")
  {
    m_minArg = min, m_maxArg = max, m_dataSize = 0;
    #ifdef FUNC_RECORD
      mp_recorder = std::unique_ptr<ArgumentRecord>(new ArgumentRecord(histSize, min, max));
    #endif
    (void) histSize; // ignore histSize if -DFUNC_RECORD isn't specified
  }

  // Evaluate the underlying std::function and optionally record the arg
  OUT_TYPE operator()(IN_TYPE x) override
  {
    #ifdef FUNC_RECORD
      mp_recorder->record_arg(x);
    #endif
    return mp_func(x);
  }

  void print_details(std::ostream& out) override;
  void print_details_json(std::ostream& out);
  ~DirectEvaluation(){};
};

/* Implementations of functions that print to the console */
template <typename IN_TYPE, typename OUT_TYPE>
inline void DirectEvaluation<IN_TYPE,OUT_TYPE>::print_details(std::ostream& out)
{
  out<< m_name << " " << m_minArg << " " << m_maxArg;
  #ifdef FUNC_RECORD
    out << std::endl;
    mp_recorder->print_details(out);
  #endif
}

template <typename IN_TYPE, typename OUT_TYPE>
inline void DirectEvaluation<IN_TYPE,OUT_TYPE>::print_details_json(std::ostream& out)
{
  using nlohmann::json;
  json jsonStats;

  jsonStats["name"] = m_name;
  jsonStats["minArg"] = m_minArg;
  jsonStats["maxArg"] = m_maxArg;
  #ifdef FUNC_RECORD
    mp_recorder->print_details_json(out);
  #endif
  out << jsonStats.dump(2) << std::endl;
}

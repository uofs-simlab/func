/*
  std::function wrapper with optional support for plotting the usage
  of a function's domain in a histogram.
  Replace your original function's usage with this class and compile
  with -DFUNC_RECORD in order to figure out what bounds you should
  set your table intervals to.

  Usage example:
    DirectEvaluation de(&function,0,10);
    double f = de(0.87354);
    // sim code
    de.print_details(); // prints max/min recorded args if FUNC_RECORD is defined

  Notes:
  - Used to store a std::function as an EvaluationImplementation which
  made ImplementationComparator easier to implement and use.
  - requires min and max args for consistency with the EvaluationImplementation.
  Once compiled with -DFUNC_RECORD, they're used as a rough guess for the
  histogram's bounds. print_details() will tell you if they were a bad guess.
  - Record where the function is being evaluated by specifying the -DFUNC_RECORD
  flag at compile time. ArgumentRecorder is an extension of this class which does
  all the work recording function arguments.
*/
#pragma once

#include "EvaluationImplementation.hpp"
#include "FunctionContainer.hpp"
#include "json.hpp"
#include <fstream> //ifstream

#ifdef FUNC_RECORD
  #include "ArgumentRecord.hpp"
#endif

template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE>
class DirectEvaluation final : public EvaluationImplementation<IN_TYPE,OUT_TYPE>
{
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  #ifdef FUNC_RECORD
    std::unique_ptr<ArgumentRecord<IN_TYPE>> mp_recorder;
  #endif
public:

  /* set base class args and set up argument recording if 
     FUNC_RECORD is defined used at compile time */
  DirectEvaluation(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container,
      IN_TYPE min = 0, IN_TYPE max = 1, unsigned int histSize = 10) :
    EvaluationImplementation<IN_TYPE,OUT_TYPE>(func_container->standard_func, "DirectEvaluation")
  {
    m_minArg = min, m_maxArg = max, m_dataSize = 0;
    #ifdef FUNC_RECORD
      mp_recorder = std::unique_ptr<ArgumentRecord<IN_TYPE>>(new ArgumentRecord<IN_TYPE>(min, max, histSize));
    #endif
    (void) histSize; // ignore histSize if -DFUNC_RECORD isn't specified
  }

  /* rebuild this class and it's arg record from a file */
  DirectEvaluation(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container,
    std::string filename)
  {
    std::ifstream file_reader(filename);
    using nlohmann::json;
    json jsonStats;
    file_reader >> jsonStats;
    m_name = jsonStats["name"].get<std::string>();
    m_minArg = jsonStats["minArg"].get<IN_TYPE>();
    m_maxArg = jsonStats["maxArg"].get<IN_TYPE>();

  #ifdef FUNC_RECORD
    // reconstruct our arg record
    mp_recorder = std::unique_ptr<ArgumentRecord<IN_TYPE>>(new ArgumentRecord<IN_TYPE>(jsonStats));
  #endif
  }

  // Evaluate the underlying std::function and optionally record the arg
  OUT_TYPE operator()(IN_TYPE x) override
  {
    #ifdef FUNC_RECORD
      mp_recorder->record_arg(x);
    #endif
    return m_func(x);
  }

  void print_details(std::ostream& out) override;
  void print_details_json(std::ostream& out) override;
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

  jsonStats["_comment"] = "FunC DirectEvaluation data";
  jsonStats["name"] = m_name;
  jsonStats["minArg"] = m_minArg;
  jsonStats["maxArg"] = m_maxArg;
  #ifdef FUNC_RECORD
    // have our ArgumentRecord add it's own data
    mp_recorder->print_details_json(jsonStats);
  #endif

  out << jsonStats.dump(2) << std::endl;
}

/*
  std::function wrapper with optional support
  for plotting the usage of a function's domain

  Usage example:
    DirectEvaluation de(&function,0,10);
    double f = de(0.87354);

  Notes:
  - basically just a wrapper around a function
  - requires a max and min value (for consistency with EvaluationImplementation)
  - Record where the function is being evaluated by specifying the -DFUNC_RECORD flag at compile time
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
private:
  #ifdef FUNC_RECORD
    std::unique_ptr<ArgumentRecord> mp_recorder;
  #endif
public:

  // set base class args and setup argument recording if 
  // -DFUNC_RECORD is used at compile time
  DirectEvaluation(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container,
      IN_TYPE min = 0, IN_TYPE max = 1, unsigned int histSize = 10) :
    EvaluationImplementation<IN_TYPE,OUT_TYPE>(func_container->standard_func, "DirectEvaluation")
  {
    this->m_minArg = min, this->m_maxArg = max, this->m_dataSize = 0;
    #ifdef FUNC_RECORD
      mp_recorder = std::unique_ptr<ArgumentRecord>(new ArgumentRecord(histSize, min, max));
    #endif
    (void) histSize; // don't let the compiler complain if -DFUNC_RECORD isn't specified
  }

  // Evaluate the underlying std::function and optionally record the arg
  OUT_TYPE operator()(IN_TYPE x) override
  {
    #ifdef FUNC_RECORD
      mp_recorder->record_arg(x);
    #endif
    return this->mp_func(x);
  }

  void print_details_json(std::ostream& out)
  {
    using nlohmann::json;
    json jsonStats;

    jsonStats["name"] = this->m_name;
    jsonStats["minArg"] = this->m_minArg;
    jsonStats["maxArg"] = this->m_maxArg;
    #ifdef FUNC_RECORD
      mp_recorder->print_details_json(out);
    #endif
    out << jsonStats.dump(2) << std::endl;
  }

  void print_details(std::ostream& out) override
  {
    out<< this->m_name << " " << this->m_minArg << " " << this->m_maxArg;
    #ifdef FUNC_RECORD
      out << std::endl;
      mp_recorder->print_details(out);
    #endif
  }

  ~DirectEvaluation(){};
};

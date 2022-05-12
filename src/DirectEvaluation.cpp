/* Implementation of a function wrapper with optional support for plotting usage of the function's domain*/
#include "DirectEvaluation.hpp"
#include "json.hpp"

#ifdef FUNC_RECORD
#include "ArgumentRecord.hpp"
#endif

DirectEvaluation::DirectEvaluation(FunctionContainer *func_container, double min, double max, unsigned int histSize) : 
  EvaluationImplementation(func_container->double_func, "DirectEvaluation")
{
  /* Base class default variables */
  this->m_minArg = min; this->m_maxArg = max;
  this->m_dataSize = 0;

  #ifdef FUNC_RECORD
    mp_recorder=std::unique_ptr<ArgumentRecord>(new ArgumentRecord(histSize, min, max));
  #endif
  (void) histSize; // don't let the compiler complain if -DFUNC_RECORD isn't specified
}

void DirectEvaluation::print_details_json(std::ostream& out)
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

double DirectEvaluation::operator()(double x)
{
  #ifdef FUNC_RECORD
    mp_recorder->record_arg(x);
  #endif
  return (*mp_func)(x);
}

void DirectEvaluation::print_details(std::ostream& out)
{
  out<< m_name << " " << m_minArg << " " << m_maxArg;
  #ifdef FUNC_RECORD
    out << std::endl;
    mp_recorder->print_details(out);
  #endif
}

DirectEvaluation::~DirectEvaluation(){}

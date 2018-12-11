/* Implementation of a Uniform Lookup table with linear interpolation */
#include "DirectEvaluation.hpp"

DirectEvaluation::DirectEvaluation(EvaluationFunctor<double,double> *func, double min, double max) : EvaluationImplementation(func, "DirectEvaluation")
{
  /* Base class default variables */
  this->m_minArg = min; this->m_maxArg = max;
  this->m_dataSize = 0;
}

double DirectEvaluation::operator()(double x)
{
  return (*mp_func)(x);
}

void DirectEvaluation::print_details(std::ostream& out)
{
  out << m_name << " " << m_minArg << " " << m_maxArg << " ";
}

DirectEvaluation::~DirectEvaluation()
{
}

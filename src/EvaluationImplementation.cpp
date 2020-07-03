/* Implementation of EvaluationImplementation */
#include "EvaluationImplementation.hpp"

EvaluationImplementation::EvaluationImplementation(std::function<double(double)> func, std::string name) : m_name(name), mp_func(func), m_minArg(0), m_maxArg(0){}
/* public access of protected data */
double EvaluationImplementation::min_arg(){ return m_minArg; }
double EvaluationImplementation::max_arg(){ return m_maxArg; }
unsigned EvaluationImplementation::order(){ return m_order; }
unsigned EvaluationImplementation::size(){ return m_dataSize; }
std::string EvaluationImplementation::name(){ return m_name; }
std::function<double(double)> EvaluationImplementation::function(){ return mp_func; }

/* Implementation of a function wrapper with optional support for plotting usage of the function's domain*/
#include "DirectEvaluation.hpp"
#include "json.hpp"

DirectEvaluation::DirectEvaluation(EvaluationFunctor<double,double> *func, double min, double max, unsigned int histSize) : 
  EvaluationImplementation(func, "DirectEvaluation")
{
  /* Base class default variables */
  this->m_minArg = min; this->m_maxArg = max;
  this->m_dataSize = 0;

  #ifdef FUNC_RECORD
  /* variables needed for recording function arguments */
  this->m_histogram=new unsigned int[histSize];
  std::fill(m_histogram, m_histogram+histSize, 0);
  this->m_histSize=histSize;
  this->m_peak=0;
  #endif
  (void) histSize; // don't let the compiler complain if -DFUNC_RECORD isn't specified
}

void DirectEvaluation::record_arg(double x)
{
#ifdef FUNC_RECORD
  unsigned int index = calcIndex(x);
  m_histogram[index]++;
  unsigned int peak_index = calcIndex(m_peak);
  if(m_histogram[index]>m_histogram[peak_index])
    m_peak=x;
#endif
}

unsigned int DirectEvaluation::calcIndex(double x)
{
  return ((unsigned int) (m_histSize*(x-m_minArg)/(m_maxArg+1-m_minArg)))%m_histSize;
}

void DirectEvaluation::print_details_json(std::ostream& out)
{
  using nlohmann::json;
  json jsonStats;

  jsonStats["name"] = m_name;
  jsonStats["minArg"] = m_minArg;
  jsonStats["maxArg"] = m_maxArg;
#ifdef FUNC_RECORD
  jsonStats["_comment"] = "Histogram of function evaluations.";
  jsonStats["histogramSize"] = m_histSize;
  jsonStats["histogram"]=histogramToString();
  jsonStats["peak"] = m_histogram[calcIndex(m_peak)];
  jsonStats["peakLoc"] = calcIndex(m_peak);
  //other statistics
#endif
  out << jsonStats.dump(2) << std::endl;
}

double DirectEvaluation::operator()(double x)
{
  // FUNC_RECORD_ARG is a null operator is -DFUNC_RECORD is not 
  // specified at compile time 
  FUNC_RECORD_ARG(x);
  return (*mp_func)(x);
}

std::string DirectEvaluation::histogramToString()
{
#ifdef FUNC_RECORD
  std::string hist_str;
  hist_str=hist_str+std::to_string(m_minArg)+'\n';
  for(unsigned int i=0; i<m_histSize; i++){
    for(unsigned int j=0; j<(unsigned int)(15*m_histogram[i]/m_histogram[calcIndex(m_peak)]); j++)
      hist_str=hist_str+'*';
    hist_str=hist_str+'\n';
  }
  hist_str=hist_str+std::to_string(m_maxArg);
  return hist_str;
#else
  return "";
#endif
}

void DirectEvaluation::print_details(std::ostream& out)
{
#ifdef FUNC_RECORD
  // output a histogram if the function has been evaluated at least once
  unsigned int peak_idx=calcIndex(m_peak);
  if(m_histogram[peak_idx]==0){
    out << m_name << " " << m_minArg << " " << m_maxArg << " " << std::endl;
    return;
  }
  out << m_name << std::endl;
  out<<"histogram of evaluations:"<<std::endl;
  out<<histogramToString()<<std::endl;
  out<<"The peak is located at x=" << m_peak << " with ";
  out<< m_histogram[peak_idx] << " evaluations." << std::endl;
#else 
  out<< m_name << " " << m_minArg << " " << m_maxArg << " " << std::endl;
#endif
}

DirectEvaluation::~DirectEvaluation()
{
#ifdef FUNC_RECORD
  delete m_histogram;
#endif
}

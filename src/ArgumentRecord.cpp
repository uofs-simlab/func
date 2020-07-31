#include "ArgumentRecord.hpp"
#include "json.hpp"
#define CALC_INDEX(X) (((unsigned int) (m_histSize*(X-m_min)/(m_max+1-m_min)))%m_histSize)

ArgumentRecord::ArgumentRecord(unsigned int histSize, double min, double max) : 
  m_histSize(histSize), m_min(min), m_max(max)
{
  /* variables needed for recording function arguments */
  this->mp_histogram=new unsigned int[histSize];
  this->mp_histogram_mutex.reset(new std::mutex[histSize]);

  // init the array to zeros so we can easily increment each entry
  std::fill(mp_histogram, mp_histogram+histSize, 0);
  this->m_peak=0;
}

void ArgumentRecord::record_arg(double x)
{
  unsigned int index = CALC_INDEX(x);
  // lock only exists for the scope of this function
  std::lock_guard<std::mutex> lock1(mp_histogram_mutex[index]);
  mp_histogram[index]++;

  std::lock_guard<std::mutex> lock2(m_peak_mutex);
  unsigned int peak_index = CALC_INDEX(m_peak);
  if(mp_histogram[index]>mp_histogram[peak_index])
    m_peak=x;
}

void ArgumentRecord::print_details_json(std::ostream& out)
{
  using nlohmann::json;
  json jsonStats;

  jsonStats["_comment"] = "Histogram of function evaluations.";
  jsonStats["histogramSize"] = m_histSize;
  for(unsigned int i=0; i<m_histSize; i++)
    jsonStats["histogram["+std::to_string(i)+"]"]=mp_histogram[i];
  jsonStats["peak_idx"] = CALC_INDEX(m_peak);
  jsonStats["max"] = m_max;
  jsonStats["min"] = m_min;
  //other statistics
  out << jsonStats.dump(2) << std::endl;
}

/* Create a horizontal string of the histogram such that the tallest row is 15 asterisks tall */
std::string ArgumentRecord::to_string()
{
  // avoid division by zero by printing nothing if the histogram is empty
  unsigned int peak_index=CALC_INDEX(m_peak);
  if(mp_histogram[peak_index]==0)
    return "";
  
  std::string hist_str;
  hist_str=std::to_string(m_min)+'\n';
  for(unsigned int i=0; i<m_histSize; i++){
    for(unsigned int j=0; j<(unsigned int)(15*mp_histogram[i]/mp_histogram[peak_index]); j++)
      hist_str=hist_str+'*';
    hist_str=hist_str+'\n';
  }
  hist_str=hist_str+std::to_string(m_max);
  return hist_str;
}

void ArgumentRecord::print_details(std::ostream& out)
{
  out<<"histogram: "<<std::endl;
  out<<this->to_string()<<std::endl;
  out<<"The peak is located at x=" << m_peak << " with ";
  out<< mp_histogram[CALC_INDEX(m_peak)] << " evaluations." << std::endl;
}

// free up the space used by the histogram
ArgumentRecord::~ArgumentRecord(){ delete mp_histogram; }

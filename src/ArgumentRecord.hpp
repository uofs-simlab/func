/*
  Helper class which acts as an extension to any existing EvaluationImplementation. 
  Wraps an array of unsigned int which acts as a histogram
  for recording the usage of a function's domain.

  Currently, this class doesn't end up in the current version 
  of FunC unless the -DFUNC_RECORD flag is specified at compile time.
  See DirectEvaluation for example usage.
*/
#pragma once
#include <string> // to_string()
#include <memory>
#include <mutex>
#include "json.hpp"

// macro used to decide where an arg should be placed in the histogram
#define COMPUTE_INDEX(X) (((unsigned int) (m_histSize*(X-m_min)/(m_max+1-m_min)))%m_histSize)

template <typename IN_TYPE>
class ArgumentRecord
{
  // Histogram used to record locations of function evaluations
  // and other helper helper vars
  std::unique_ptr<unsigned int[]> mp_histogram;
  std::unique_ptr<std::mutex[]>   mp_histogram_mutex;
  unsigned int m_histSize;

  // min and max should always be the same as the EvaluationImpl
  // containing the ArgumentRecord
  IN_TYPE m_min;
  IN_TYPE m_max;

  // vars containing any statistics each with a corresponding
  // mutex to keep things threadsafe
  IN_TYPE m_peak;
  std::mutex m_peak_mutex;

  public:
    ArgumentRecord(unsigned int histSize, IN_TYPE min, IN_TYPE max) :
      m_histSize(histSize), m_min(min), m_max(max)
    {
      /* variables needed for recording function arguments
         init the histogram to 0 for easy incrementing*/
      raw_hist = new unsigned int[histSize];
      std::fill(mp_histogram, mp_histogram+histSize, 0);
      // make this pointer unique and threadsafe
      this->mp_histogram = std::unique_ptr<unsigned int[]>(raw_hist);
      this->mp_histogram_mutex.reset(new std::mutex[histSize]);
      this->m_peak=0;
    }

    // place x in the histogram
    void record_arg(IN_TYPE x)
    {
      unsigned int index = COMPUTE_INDEX(x);
      // lock only exists for the scope of this function
      std::lock_guard<std::mutex> lock1(mp_histogram_mutex[index]);
      mp_histogram[index]++;

      std::lock_guard<std::mutex> lock2(m_peak_mutex);
      unsigned int peak_index = COMPUTE_INDEX(m_peak);
      if(mp_histogram[index]>mp_histogram[peak_index])
        m_peak=x;
    }

    // make a string representation of the histogram
    std::string to_string();

    // print out the various fields in this class
    void print_details(std::ostream& out);
    void print_details_json(std::ostream& out);

    // free up the space used by the histogram
    ~ArgumentRecord(){}
};

template <typename IN_TYPE>
inline std::string ArgumentRecord<IN_TYPE>::to_string()
{
  // avoid division by zero by printing nothing if the histogram is empty
  unsigned int peak_index=COMPUTE_INDEX(m_peak);
  if(mp_histogram[peak_index]==0)
    return "";
  
  // print out the histogram horizontally such that the longest row is 15 stars wide
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

template <typename IN_TYPE>
inline void ArgumentRecord<IN_TYPE>::print_details_json(std::ostream& out)
{
  using nlohmann::json;
  json jsonStats;

  jsonStats["_comment"] = "Histogram of function evaluations.";
  jsonStats["histogramSize"] = m_histSize;
  for(unsigned int i=0; i<m_histSize; i++)
    jsonStats["histogram["+std::to_string(i)+"]"]=mp_histogram[i];
  jsonStats["peak_idx"] = COMPUTE_INDEX(m_peak);
  jsonStats["max"] = m_max;
  jsonStats["min"] = m_min;
  //other statistics
  out << jsonStats.dump(2) << std::endl;
}

template <typename IN_TYPE>
inline void ArgumentRecord<IN_TYPE>::print_details(std::ostream& out)
{
  out<<"histogram: "<<std::endl;
  out<<this->to_string()<<std::endl;
  out<<"The peak is located at x=" << m_peak << " with ";
  out<< mp_histogram[COMPUTE_INDEX(m_peak)] << " evaluations." << std::endl;
}

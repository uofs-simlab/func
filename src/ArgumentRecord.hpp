/*
  Helper class which acts as an extension to any existing
  EvaluationImplementation. Wraps a vector of unsigned
  int which acts as a histogram for recording the
  usage of a function's domain.

  Notes:
  - Currently, this class doesn't end up in the current version 
  of FunC unless the -DFUNC_RECORD flag is specified at compile time.
  See DirectEvaluation for example usage.
  - Designed assuming it will be a private member var so some functions
  won't do everything on their own
  *specifically print_details_json assumes the encapsulating 
  class will do most of the work for it*
*/
#pragma once
#include <string> // to_string()
#include <memory>
#include <mutex>
#include <fstream>
#include "json.hpp"

// macro used to decide where an arg should be placed in the histogram
#define COMPUTE_INDEX(X) \
  (((unsigned int) (m_histSize*(X-m_minArg)/(m_maxArg+1-m_minArg)))%m_histSize)

template <typename IN_TYPE>
class ArgumentRecord
{
  // Histogram used to record locations of function evaluations
  // and other helper helper vars
  std::vector<unsigned int>     mv_histogram;
  std::unique_ptr<std::mutex[]> mp_histogram_mutex;
  unsigned int                  m_histSize;

  // min and max should always be the same as the EvaluationImpl
  // containing the ArgumentRecord
  IN_TYPE m_minArg;
  IN_TYPE m_maxArg;

  // vars containing any statistics each with a corresponding
  // mutex to keep things threadsafe
  IN_TYPE m_peak_arg;
  std::mutex m_peak_mutex;

  // Record the extreme args so the user
  // knows what bounds to use for their tables
  IN_TYPE m_max_recorded;
  std::mutex m_max_recorded_mutex;

  IN_TYPE m_min_recorded;
  std::mutex m_min_recorded_mutex;

  public:
    ArgumentRecord(unsigned int histSize, IN_TYPE min, IN_TYPE max) :
      m_histSize(histSize), m_minArg(min), m_maxArg(max)
    {
      /* variables needed for recording function arguments
         init the histogram to 0 for easy incrementing*/
      mv_histogram = std::vector<unsigned int>(histSize, 0);
      mp_histogram_mutex = std::unique_ptr<std::mutex[]>(new std::mutex[histSize]);

      // naive initial member arg values
      m_peak_arg=0;
      m_max_recorded=std::numeric_limits<IN_TYPE>::min();
      m_min_recorded=std::numeric_limits<IN_TYPE>::max();
    }

    ~ArgumentRecord(){}

    // Rebuild our argument record
    // Assuming the encapsulating EvaluationImplementation
    // gave us a valid json file
    ArgumentRecord(nlohmann::json jsonStats)
    {
      m_minArg = jsonStats["ArgumentRecord"]["minArg"].get<IN_TYPE>();
      m_maxArg = jsonStats["ArgumentRecord"]["maxArg"].get<IN_TYPE>();

      m_histSize = jsonStats["ArgumentRecord"]["histogramSize"].get<unsigned int>();
      for(unsigned int i=0; i<m_histSize; i++)
        mv_histogram[i] = jsonStats["ArgumentRecord"]["histogram"][std::to_string(i)].get<unsigned int>();
      // restore our histogram's thread safety
      mp_histogram_mutex = std::vector<std::unique_ptr<std::mutex>>(m_histSize, std::unique_ptr<std::mutex>(new std::mutex));

      m_peak_arg = jsonStats["ArgumentRecord"]["peak_arg"].get<IN_TYPE>();
      m_max_recorded = jsonStats["ArgumentRecord"]["m_min_recorded"].get<IN_TYPE>();
      m_min_recorded = jsonStats["ArgumentRecord"]["m_max_recorded"].get<IN_TYPE>();
    }

    // place x in the histogram
    void record_arg(IN_TYPE x)
    {
      // Record x if it's within our table's limits
      // note: each std::lock_guard only exists for the duration of its scope
      if(m_minArg <= x && x <= m_maxArg){
        int index = COMPUTE_INDEX(x);
        {
          std::lock_guard<std::mutex> lock(mp_histogram_mutex[index]);
          mv_histogram[index]++;
        }

        {
          std::lock_guard<std::mutex> lock(m_peak_mutex);
          unsigned int peak_index = COMPUTE_INDEX(m_peak_arg);
          if(mv_histogram[index] > mv_histogram[peak_index])
            m_peak_arg = x;
        } 
      }

      {
        std::lock_guard<std::mutex> lock(m_max_recorded_mutex);
        if(m_max_recorded < x)
          m_max_recorded = x;
      }

      {
        std::lock_guard<std::mutex> lock(m_min_recorded_mutex);
        if(m_min_recorded < x)
          m_min_recorded = x;
      }
    }

    // make a string representation of the histogram
    std::string to_string();

    // print out the various fields in this class
    void print_details(std::ostream& out);

    // add to the provided json object
    void print_details_json(nlohmann::json& jsonStats);
};

/* Implementations of functions that print to the console */
template <typename IN_TYPE>
inline std::string ArgumentRecord<IN_TYPE>::to_string()
{
  // avoid division by zero by printing nothing if the histogram is empty
  unsigned int peak_index=COMPUTE_INDEX(m_peak_arg);
  if(mv_histogram[peak_index]==0)
    return "";
  
  // print out the histogram horizontally such that the longest row is 15 stars wide
  std::string hist_str;
  hist_str=std::to_string(m_minArg)+'\n';
  for(unsigned int i=0; i<m_histSize; i++){
    for(unsigned int j=0; j<(unsigned int)(15*mv_histogram[i]/mv_histogram[peak_index]); j++)
      hist_str=hist_str+'*';
    hist_str=hist_str+'\n';
  }
  hist_str=hist_str+std::to_string(m_maxArg);
  return hist_str;
}

template <typename IN_TYPE>
inline void ArgumentRecord<IN_TYPE>::print_details_json(nlohmann::json& jsonStats)
{
  // Assume the caller has the same max/min
  jsonStats["ArgumentRecord"]["_comment"] = "Histogram of function evaluations.";
  jsonStats["ArgumentRecord"]["histogramSize"] = m_histSize;
  for(unsigned int i=0; i<m_histSize; i++)
    jsonStats["ArgumentRecord"]["histogram"][std::to_string(i)]=mv_histogram[i];
  jsonStats["ArgumentRecord"]["peak_arg"] = m_peak_arg;
  jsonStats["ArgumentRecord"]["m_min_recorded"] = m_min_recorded;
  jsonStats["ArgumentRecord"]["m_max_recorded"] = m_max_recorded;
  
  // insert more statistics here
}

template <typename IN_TYPE>
inline void ArgumentRecord<IN_TYPE>::print_details(std::ostream& out)
{
  out << "histogram: \n";
  out << this->to_string() << "\n";
  out << "The peak_arg is located at x=" << m_peak_arg << " with ";
  out << mv_histogram[COMPUTE_INDEX(m_peak_arg)] << " evaluations.\n";
  out << "The largest argument recorded was x=" << m_max_recorded << "\n";
  out << "The most negative argument recorded was x=" << m_min_recorded << std::endl;
}

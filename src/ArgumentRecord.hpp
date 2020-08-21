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
#include <limits> // numeric_limits<IN_TYPE>::max() / min()
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

  // set min/max to the min/max of the encapsulating EvaluationImplementation.
  // Or use std::numeric_limits<IN_TYPE>::max()/ min()
  // to make the histogram account for all possible values of x.
  // At the price of making the histogram less readable
  IN_TYPE m_minArg;
  IN_TYPE m_maxArg;

  // vars containing any statistics each with a corresponding
  // mutex to keep things threadsafe
  IN_TYPE m_peak_index; // the index of the bucket with the largest count
  std::mutex m_peak_mutex;

  // the sum of histogram elements and args outside the histogram's range
  unsigned int m_out_of_bounds;
  std::mutex m_out_of_bounds_mutex;

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
      m_peak_index = 0.0;
      m_max_recorded = -std::numeric_limits<IN_TYPE>::max();
      m_min_recorded =  std::numeric_limits<IN_TYPE>::max();
    }

    ~ArgumentRecord(){}

    // Rebuild our argument record
    // Assuming the encapsulating EvaluationImplementation
    // gave us a valid json object
    ArgumentRecord(nlohmann::json jsonStats)
    {
      m_minArg = jsonStats["ArgumentRecord"]["minArg"].get<IN_TYPE>();
      m_maxArg = jsonStats["ArgumentRecord"]["maxArg"].get<IN_TYPE>();

      m_histSize = jsonStats["ArgumentRecord"]["histogramSize"].get<unsigned int>();
      for(unsigned int i=0; i<m_histSize; i++)
        mv_histogram[i] = jsonStats["ArgumentRecord"]["histogram"][std::to_string(i)].get<unsigned int>();
      // restore our histogram's thread safety
      mp_histogram_mutex = std::vector<std::unique_ptr<std::mutex>>(m_histSize, std::unique_ptr<std::mutex>(new std::mutex));

      m_peak_index = jsonStats["ArgumentRecord"]["peak_index"].get<unsigned int>();
      m_max_recorded = jsonStats["ArgumentRecord"]["m_min_recorded"].get<IN_TYPE>();
      m_min_recorded = jsonStats["ArgumentRecord"]["m_max_recorded"].get<IN_TYPE>();
    }

    /* place x in the histogram */
    void record_arg(IN_TYPE x)
    {
      // Record x if it's within our table's limits
      // note: each std::lock_guard only exists for the duration of its scope
      if(m_minArg <= x && x <= m_maxArg){
        int m_index = COMPUTE_INDEX(x);
        {
          std::lock_guard<std::mutex> lock(mp_histogram_mutex[x_index]);
          mv_histogram[x_index]++;
        }

        {
          std::lock_guard<std::mutex> lock(m_peak_mutex);
          if(mv_histogram[x_index] > mv_histogram[m_peak_index])
            m_peak_index = x_index;
        } 
      }else{
        // count x if it doesn't qualify
        std::lock_guard<std::mutex> lock(m_out_of_bounds_mutex);
        m_out_of_bounds++;
      }

      {
        std::lock_guard<std::mutex> lock(m_max_recorded_mutex);
        if(m_max_recorded < x)
          m_max_recorded = x;
      }

      {
        std::lock_guard<std::mutex> lock(m_min_recorded_mutex);
        if(x < m_min_recorded)
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
  if(mv_histogram[m_peak_index]==0)
    return "";
  
  // print out the histogram horizontally such that the longest row is 15 stars wide
  std::string hist_str;
  hist_str=std::to_string(m_minArg)+'\n';
  for(unsigned int i=0; i<m_histSize; i++){
    for(unsigned int j=0; j<(unsigned int)(15*mv_histogram[i]/mv_histogram[m_peak_index]); j++)
      hist_str=hist_str+'*';
    hist_str=hist_str+'\n';
  }
  hist_str=hist_str+std::to_string(m_maxArg);
  return hist_str;
}

template <typename IN_TYPE>
inline void ArgumentRecord<IN_TYPE>::print_details_json(nlohmann::json& jsonStats)
{
  jsonStats["ArgumentRecord"]["_comment"] = "Histogram of function evaluations.";
  jsonStats["ArgumentRecord"]["minArg"] = m_minArg;
  jsonStats["ArgumentRecord"]["maxArg"] = m_maxArg;
  jsonStats["ArgumentRecord"]["histogramSize"] = m_histSize;
  for(unsigned int i=0; i<m_histSize; i++)
    jsonStats["ArgumentRecord"]["histogram"][std::to_string(i)]=mv_histogram[i];
  jsonStats["ArgumentRecord"]["peak_index"] = m_peak_index;
  jsonStats["ArgumentRecord"]["m_min_recorded"] = m_min_recorded;
  jsonStats["ArgumentRecord"]["m_max_recorded"] = m_max_recorded;
  
  // insert more statistics here
}

template <typename IN_TYPE>
inline void ArgumentRecord<IN_TYPE>::print_details(std::ostream& out)
{
  unsigned int total_recorded = 0.0;
  for(unsigned int i=0; i<m_histSize; i++)
    total_recorded += mv_histogram[i];

  out << "histogram: \n";
  out << this->to_string() << "\n";
  out << m_out_of_bounds + total_recorded << " total args were sampled. Of those, "
      << total_recorded << " were recorded by the histogram.\n"
  out << "Recorded args were sampled the most often from the subinterval"
      << "["
      << m_minArg + (m_maxArg - m_minArg)*m_peak_index/m_histSize << ", "
      << m_minArg + (m_maxArg - m_minArg)*(m_peak_index+1)/m_histSize << ") "
      << "with " << mv_histogram[m_peak_index] << " evaluations.\n";
  out << "The largest argument recorded was x=" << m_max_recorded << "\n";
  out << "The least argument recorded was x=" << m_min_recorded << std::endl;
}

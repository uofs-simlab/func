/*
  Helper class which acts as an extension to any existing
  EvaluationImplementation. Wraps a vector of unsigned
  int which acts as a histogram for recording the
  usage of a function's domain.

  Notes:
  - Currently, this class is only used in the current version 
  of FunC if the -DFUNC_DEBUG flag is specified at compile time.
  - Argument recording is threadsafe
  - See DirectEvaluation.hpp and FailureProofTable.hpp for example usage.
  - This is designed to be a private member variable of some class
*/
#pragma once
#include <string> // to_string()
#include <memory>
#include <fstream>
#include <limits> // numeric_limits<IN_TYPE>::max() / min()
#include "json.hpp"

/* TODO 
   Since the histogram will never have that many buckets,
   I think we could get away with just making every one of
   this class's member variables threadprivate. */

// can change FuncMutex if we decide to use std::mutex in
// the case where openmp isn't available
#ifdef _OPENMP // defined when -fopenmp is used at compile time
#include <omp.h>
class FuncMutex
{
  omp_lock_t m_lock;
public:
  FuncMutex(){ omp_init_lock(&m_lock); }
  ~FuncMutex(){ omp_destroy_lock(&m_lock); }
  void lock(){ omp_set_lock(&m_lock); }
  void unlock(){ omp_unset_lock(&m_lock); }
};

class FuncScopedLock
{
  FuncMutex& m_mutex;
  // prevent copying
  void operator=(const FuncScopedLock&);
  FuncScopedLock(const FuncScopedLock&);
public:
  explicit FuncScopedLock(FuncMutex& mutex) :
    m_mutex(mutex) { m_mutex.lock(); }
  // release the lock when FuncScopedLock goes out of scope
  ~FuncScopedLock(){ m_mutex.unlock(); }
};
#else
// just use empty locks for now
class FuncMutex
{
public:
  FuncMutex(){}
  ~FuncMutex(){}
  void lock(){}
  void unlock(){}
};

class FuncScopedLock
{
  // prevent copying
  void operator=(const FuncScopedLock&);
  FuncScopedLock(const FuncScopedLock&);
public:
  explicit FuncScopedLock(FuncMutex& mutex){}
  ~FuncScopedLock(){}
};
#endif // _OPENMP

// macro used to decide where an arg should be placed in the histogram
#define COMPUTE_INDEX(X) \
  (((unsigned int) (m_histSize*(X-m_minArg)/(m_maxArg+1-m_minArg)))%m_histSize)

template <typename IN_TYPE>
class ArgumentRecord
{
  // Histogram used to record locations of function evaluations
  // and other helper helper vars
  std::vector<unsigned int> mv_histogram;
  std::vector<FuncMutex>    mv_histogram_mutex;
  unsigned int              m_histSize;

  // Set table bounds. Can be altered to result in nicer output
  IN_TYPE m_minArg;
  IN_TYPE m_maxArg;

  /* vars containing any statistics */
  // the index of the bucket with the largest count
  unsigned int m_peak_index;

  // the number of elements outside the histogram's range
  unsigned int m_num_out_of_bounds;

  // Record the extreme args to help the user
  // decide what bounds to use for their tables
  IN_TYPE m_max_recorded;
  IN_TYPE m_min_recorded;

  public:
    ArgumentRecord(IN_TYPE min, IN_TYPE max, unsigned int histSize) :
      m_minArg(min), m_maxArg(max), m_histSize(histSize)
    {
      /* variables needed for recording function arguments.
         Init the entire histogram to 0 for easy incrementing */
      mv_histogram = std::vector<unsigned int>(histSize, 0);
      mv_histogram_mutex = std::vector<FuncMutex>(histSize);

      // naive initial member arg values
      m_peak_index = 0;
      m_num_out_of_bounds = 0;
      m_max_recorded = std::numeric_limits<IN_TYPE>::lowest();
      m_min_recorded = std::numeric_limits<IN_TYPE>::max();
    }

#ifdef FUNC_DEBUG
    ~ArgumentRecord(){print_details(std::cerr);}
#else
    ~ArgumentRecord(){}
#endif

    /* Rebuild our argument record
       Note: Assuming the encapsulating
       EvaluationImplementation gave us a valid json object */
    ArgumentRecord(nlohmann::json jsonStats)
    {
      m_minArg = jsonStats["ArgumentRecord"]["minArg"].get<IN_TYPE>();
      m_maxArg = jsonStats["ArgumentRecord"]["maxArg"].get<IN_TYPE>();

      m_histSize = jsonStats["ArgumentRecord"]["histogramSize"].get<unsigned int>();
      for(unsigned int i=0; i<m_histSize; i++)
        mv_histogram[i] = jsonStats["ArgumentRecord"]["histogram"][std::to_string(i)].get<unsigned int>();
      // rebuild the histogram's thread safety
      mv_histogram_mutex = std::vector<FuncMutex>(m_histSize);

      m_peak_index   = jsonStats["ArgumentRecord"]["peak_index"].get<unsigned int>();
      m_max_recorded = jsonStats["ArgumentRecord"]["m_max_recorded"].get<IN_TYPE>();
      m_min_recorded = jsonStats["ArgumentRecord"]["m_min_recorded"].get<IN_TYPE>();
    }

    /* place x in the histogram */
    void record_arg(IN_TYPE x)
    {
      // Record x if it's within our histogram's limits
      if(m_minArg <= x && x <= m_maxArg){
        int x_index = COMPUTE_INDEX(x);

        // manually lock bucket x_index until we leave this scope
        FuncScopedLock lock(mv_histogram_mutex[x_index]);
        //#pragma omp critical
        mv_histogram[x_index]++;

        #pragma omp critical
        {
          if(mv_histogram[x_index] > mv_histogram[m_peak_index]){
            m_peak_index = x_index;
          }
        }

      }else{
        // x is out of bounds
        #pragma omp critical
        m_num_out_of_bounds++;
      }

      #pragma omp critical
      {
        if(m_max_recorded < x)
          m_max_recorded = x;
      }

      #pragma omp critical
      {
        if(x < m_min_recorded)
          m_min_recorded = x;
      }

      // record more statistics here
    }

    // make a string representation of the histogram
    std::string to_string()
    {
      // avoid division by zero by printing nothing if the histogram is empty
      if(mv_histogram[m_peak_index]==0)
        return "";
      
      // print out the histogram horizontally such that
      // 1. the longest row is 15 stars wide
      // 2. If a bucket recorded _any_ args then it gets at least 1 star
      std::string hist_str;
      hist_str=std::to_string(m_minArg)+'\n';
      for(unsigned int i=0; i<m_histSize; i++){
        for(unsigned int j=0; j<(unsigned int)ceil(15*mv_histogram[i]/mv_histogram[m_peak_index]); j++)
          hist_str=hist_str+'*';
        hist_str=hist_str+'\n';
      }
      hist_str=hist_str+std::to_string(m_maxArg);
      return hist_str;
    }

    // human readable output of the histogram and other statistics
    void print_details(std::ostream& out)
    {
      unsigned int total_recorded = 0;
      for(unsigned int i=0; i<m_histSize; i++)
        total_recorded += mv_histogram[i];

      if(m_num_out_of_bounds + total_recorded == 0){
        out << "No arguments were recorded by arg record" << "\n";
        return;
      }


      out << "histogram: \n";
      out << this->to_string() << "\n";
      out << m_num_out_of_bounds + total_recorded << " total args were sampled. Of those, "
          << total_recorded << " were recorded by the histogram.\n";
      out << "Recorded args were sampled the most often from the subinterval "
          << "["
          << m_minArg + (m_maxArg - m_minArg)*m_peak_index/(IN_TYPE)m_histSize << ", "
          << m_minArg + (m_maxArg - m_minArg)*(m_peak_index+1)/(IN_TYPE)m_histSize << ") "
          << "with " << mv_histogram[m_peak_index] << " evaluations.\n";
      out << "The largest argument seen was x=" << m_max_recorded << "\n";
      out << "The lowest argument seen was x=" << m_min_recorded << std::endl;
    }

    // print each field in this class to the given ostream
    void print_details_json(nlohmann::json& jsonStats)
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
};


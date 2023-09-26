/* TODO USE BOOST HISTOGRAM INSTEAD
 *
 *
 *
  Helper class which acts as an extension to any existing
  EvaluationImplementation. Wraps a vector of unsigned
  int which acts as a histogram for recording the
  usage of a function's domain.

  Notes:
  - This class is only used in the current version 
  of FunC if the -DFUNC_DEBUG flag is specified at compile time.
  - Argument recording is threadsafe
  - See DirectEvaluation.hpp and FailureProofTable.hpp for example usage.
  - This is designed to be a private member variable of some class

  TODO this class should support to_json & from_json
*/
#pragma once
#include <string> // to_string()
#include <memory>
#include <fstream>
#include <limits> // numeric_limits<TIN>::max() / lowest()
#include <cmath> // pow, floor, ceil
#include "json.hpp"

namespace func {
/* TODO 
   The histogram will never have that many buckets so 
   we could likely get away with just making every one of
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
  explicit FuncScopedLock(FuncMutex& mutex){(void)mutex;}
  ~FuncScopedLock(){}
};
#endif // _OPENMP

// macro used to decide where an arg should be placed in the histogram
#define COMPUTE_INDEX(X) \
  (static_cast<unsigned int>(m_histSize*(X-m_minArg)/(m_maxArg+1.0-m_minArg))%m_histSize)

template <typename TIN>
class ArgumentRecord
{
  // Histogram used to record locations of function evaluations
  // and other helper helper vars
  std::vector<unsigned int> mv_histogram;
  std::vector<FuncMutex>    mv_histogram_mutex;
  unsigned int              m_histSize;

  // Set table bounds. Can be altered to result in nicer output
  TIN m_minArg;
  TIN m_maxArg;

  /* vars containing any statistics */
  // the index of the bucket with the largest count
  unsigned int m_peak_index;

  // the number of elements outside the histogram's range
  unsigned int m_num_out_of_bounds;

  // Record the extreme args to help the user
  // decide what bounds to use for their tables
  TIN m_max_recorded;
  TIN m_min_recorded;

  public:
    ArgumentRecord(TIN min, TIN max, unsigned int histSize) :
      m_minArg(min), m_maxArg(max), m_histSize(histSize)
    {
      /* variables needed for recording function arguments.
         Init the entire histogram to 0 for easy incrementing */
      mv_histogram = std::vector<unsigned int>(histSize, 0);
      mv_histogram_mutex = std::vector<FuncMutex>(histSize);

      // naive initial member arg values
      m_peak_index = 0;
      m_num_out_of_bounds = 0;
      m_max_recorded = std::numeric_limits<TIN>::lowest();
      m_min_recorded = std::numeric_limits<TIN>::max();
    }

#ifdef FUNC_DEBUG
    /* make the specific stream printed to depend on a printing function */
    ~ArgumentRecord(){ std::cerr << *this;}
#else
    ~ArgumentRecord(){}
#endif

    /* Rebuild our argument record
       Note: Assuming the encapsulating LookupTable gave us a valid json object */
    ArgumentRecord(nlohmann::json jsonStats)
    {
      m_minArg = jsonStats["ArgumentRecord"]["minArg"].get<TIN>();
      m_maxArg = jsonStats["ArgumentRecord"]["maxArg"].get<TIN>();

      m_histSize = jsonStats["ArgumentRecord"]["histogramSize"].get<unsigned int>();
      for(unsigned int i=0; i<m_histSize; i++)
        mv_histogram[i] = jsonStats["ArgumentRecord"]["histogram"][std::to_string(i)].get<unsigned int>();
      // rebuild the histogram's thread safety
      mv_histogram_mutex = std::vector<FuncMutex>(m_histSize);

      m_peak_index   = jsonStats["ArgumentRecord"]["peak_index"].get<unsigned int>();
      m_max_recorded = jsonStats["ArgumentRecord"]["m_max_recorded"].get<TIN>();
      m_min_recorded = jsonStats["ArgumentRecord"]["m_min_recorded"].get<TIN>();
    }

    /* place x in the histogram. Mimic pipeline parallelism for any statistics with only one instance */
    void record_arg(TIN x)
    {
      // Record x if it's within our histogram's limits
      if(m_minArg <= x && x <= m_maxArg){
        int x_index = COMPUTE_INDEX(x);

        // manually lock bucket x_index until we leave this scope
        FuncScopedLock lock(mv_histogram_mutex[x_index]);
        mv_histogram[x_index]++;

        /* TODO using critical sections is problematic design because only one
         * ArgumentRecord instance can access some specific critical section code at a time.
         * Use a scoped lock instead */
        #pragma omp critical(func_argrecord1)
        {
          if(mv_histogram[x_index] > mv_histogram[m_peak_index]){
            m_peak_index = x_index;
          }
        }

      }else{
        // x is out of bounds
        #pragma omp critical(func_argrecord2)
        m_num_out_of_bounds++;
      }

      #pragma omp critical(func_argrecord3)
      {
        if(m_max_recorded < x)
          m_max_recorded = x;
      }

      #pragma omp critical(func_argrecord4)
      {
        if(x < m_min_recorded)
          m_min_recorded = x;
      }

      // record more statistics here
    }

    /* std::to_string(1e-7) == "0" which is unacceptable so we'll use this code from this SO post
     * https://stackoverflow.com/questions/16605967/set-precision-of-stdto-string-when-converting-floating-point-values
     * Default is the max possible precision by so users can choose how they'll round the answer on their own */
    template <typename T>
    std::string to_string_with_precision(const T val, const int n = std::numeric_limits<T>::max_digits10) const 
    {
        std::ostringstream out;
        out.precision(n);
        out << std::scientific << val;
        return out.str();
    }

    std::string ith_interval(unsigned int i, const int n = 3) const {
      return "[" + to_string_with_precision(m_minArg + (m_maxArg - m_minArg)*i/static_cast<TIN>(m_histSize), n) + ", "
          + to_string_with_precision(m_minArg + (m_maxArg - m_minArg)*(i+1)/static_cast<TIN>(m_histSize), n) + ")";
    }


    // make a string representation of the histogram
    std::string to_string() const
    {
      // avoid division by zero by printing nothing if the histogram is empty
      if(mv_histogram[m_peak_index]==0)
        return "";
      
      // print out the histogram horizontally such that
      // 1. the longest row is 15 stars wide
      // 2. If a bucket recorded _any_ args then it gets at least 1 asterisk
      std::string hist_str;
      hist_str=to_string_with_precision(m_minArg,3)+'\n';
      for(unsigned int i=0; i<m_histSize; i++){
        unsigned int row_length = ceil(15*mv_histogram[i]/mv_histogram[m_peak_index]);
        for(unsigned int j=0; j<row_length; j++)
          hist_str=hist_str+'*';
        for(unsigned int j=row_length; j<15; j++)
          hist_str=hist_str+' ';
        hist_str=hist_str + " " + ith_interval(i) + " with " + std::to_string(mv_histogram[i]) + " evaluations\n";
      }
      hist_str=hist_str+to_string_with_precision(m_maxArg,3);
      return hist_str;
    }

    // print each field in this class to the given ostream
    void print_json(nlohmann::json& jsonStats) const
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

  TIN min_arg() const { return m_minArg; }
  TIN max_arg() const { return m_maxArg; }
  /* TODO use a standard algorithm */
  unsigned int total_recorded() const {
    unsigned int t = 0;
    for(unsigned int i=0; i<m_histSize; i++)
      t += mv_histogram[i];
    return t;
  }
  
  /* the index of the bucket with the largest count */
  unsigned int index_of_peak() const { return m_peak_index; }
  unsigned int peak() const { return mv_histogram[m_peak_index]; }

  /* the number of elements outside the histogram's range */
  unsigned int num_out_of_bounds() const { return m_num_out_of_bounds; }

  /* Return the extreme args to help the user decide what bounds to use for their LUTs */
  TIN max_recorded() const { return m_max_recorded; }
  TIN min_recorded() const { return m_min_recorded; }
};


template <typename TIN>
std::ostream& operator<<(std::ostream& out, const ArgumentRecord<TIN>& arg_record){
  unsigned int complete_total = arg_record.num_out_of_bounds() + arg_record.total_recorded();

  if(complete_total == 0){
    out << "No arguments were recorded by arg record" << "\n";
    return out;
  }

  out << "histogram: \n";
  out << arg_record.to_string() << "\n";
  out << complete_total << " total args were sampled. Of those, "
      << arg_record.total_recorded() << " were recorded by the histogram.\n";
  out << "Recorded args were sampled the most often from the subinterval "
      << arg_record.ith_interval(arg_record.index_of_peak()) << " with " << arg_record.peak() << " evaluations ("
      << 100.0*arg_record.peak()/static_cast<double>(complete_total) << "% of the total evaluations).\n";
  /* iostream rounds to like 6 digits by default but the rounding can make the min/max args too large/small 
   * which is annoying so we'll just print every digit of the output and let users round on their own */
  out << "The largest argument recorded was x=" << arg_record.to_string_with_precision(arg_record.max_recorded()) << "\n";
  out << "The lowest argument recorded was x=" << arg_record.to_string_with_precision(arg_record.min_recorded()) << std::endl;
  return out;
}


//TODO make to/from_json()
} // namespace func

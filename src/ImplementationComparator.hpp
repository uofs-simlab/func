/*
  Class for comparing variations of EvaluationImplementations

  - takes ownership of the vector of implementations passed to it
*/
#pragma once
#include "Timer.hpp"
#include "RngInterface.hpp"
#include "StdRng.hpp"
#include "EvaluationImplementation.hpp"
#include "json.hpp"

#include <memory>
#include <vector>
#include <algorithm>
#include <limits> // epsilon()
#include <cassert>
#include <iostream>

namespace func {
/*
  Data types and containers used in FunC's ImplementationComparator
*/
typedef std::vector<double> TimeContainer;

template <typename IN_TYPE, typename OUT_TYPE>
using ImplType = EvaluationImplementation<IN_TYPE,OUT_TYPE>;

template <typename IN_TYPE, typename OUT_TYPE>
using ImplContainer = std::vector<std::unique_ptr<ImplType<IN_TYPE,OUT_TYPE>>>;

/*
  ImplTimer struct attaches additional data for timing an implementation
  to the implementation
*/
template <typename IN_TYPE, typename OUT_TYPE>
struct ImplTimer
{
  /* Note the the ImplType can NOT be a reference here, because we
     want to be able to sort a container of these ImplTimers. sort
     requires operator=, and classes that have non-static reference
     members cannot implement this */
  ImplType<IN_TYPE,OUT_TYPE> *impl;
  TimeContainer evaluationTimes;
  double maxTime, minTime, meanTime;
  ImplTimer(ImplType<IN_TYPE,OUT_TYPE> *inImpl) : impl(inImpl), maxTime(0), minTime(0), meanTime(0) {};
  void append_runtime(double time){ evaluationTimes.push_back(time); };
  void compute_timing_stats()
  {
    /* get min and max execution times */
    auto extremeTimes = minmax_element( evaluationTimes.begin(), evaluationTimes.end() );
    minTime = *(extremeTimes.first);
    maxTime = *(extremeTimes.second);
    /* mean execution time */
    auto sum = std::accumulate(evaluationTimes.begin(), evaluationTimes.end(), 0.0);
    meanTime =  sum / evaluationTimes.size();
  };
  void print_timing_stats(std::ostream& out)
  {
    out << "Min " << minTime << "s"
      << " Max " << maxTime << "s"
      << " Mean " << meanTime << "s"
      << "\n";
  }
};

/* ------------------------------------------------------------------------ */
template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE>
class ImplementationComparator
{
private:

  ImplContainer<IN_TYPE,OUT_TYPE>          m_implementations;
  std::vector<ImplTimer<IN_TYPE,OUT_TYPE>> m_implTimers;
  unsigned                                 m_numberOfImplementations;

  std::vector<TimeContainer>  m_evaluationTimers;

  IN_TYPE                     m_minArg,m_maxArg;

  std::unique_ptr<OUT_TYPE[]> m_evalHolder;

  /*
    RNG for evaluations
    - By default uses a std::uniform_real_distribution<IN_TYPE>
      with the std::mt19937 variant of the std::mersenne_twister_engine
  */
  std::unique_ptr<RngInterface<IN_TYPE>> mp_sampler;
  std::unique_ptr<IN_TYPE[]>             mp_randomEvaluations;
  int                                    m_nEvals;

  struct TimingStatistics
  {
    double minTime, maxTime, meanTime;
  } m_timingStatistics;

  /* Fill mp_randomEvaluations with random points to be evaluated */
  void draw_new_sample_points()
  {
    for (int ii=0;ii<m_nEvals;++ii)
      mp_randomEvaluations[ii] = mp_sampler->get_point();
  }

  /* Time implementation evaluations */
  void run_all_single()
  {
    for (auto &itImplTimer : m_implTimers) {
      Timer evalTimer;
      for (int ii=0;ii<m_nEvals;++ii)
         m_evalHolder[ii] = (*(itImplTimer.impl))(mp_randomEvaluations[ii]);

      evalTimer.stop();
      itImplTimer.append_runtime(evalTimer.duration());
    }
  }

public:

  ImplementationComparator(ImplContainer<IN_TYPE,OUT_TYPE> &inImpl, int nEvals = 100000,
      unsigned int seed = 2017, std::unique_ptr<RngInterface<IN_TYPE>> inRng = nullptr);
  ~ImplementationComparator(){}

  /* Run timings with different set of random arguments */
  void run_timings(int nRuns = 1)
  {
    for (int ii=0;ii<nRuns;++ii) {
      draw_new_sample_points();
      run_all_single();
    }
  }

  /* Obtain a const reference to the current array of evaluation points */
  const std::unique_ptr<IN_TYPE> sample_points(){ return mp_randomEvaluations; }

  /* Compute fastest and slowest times */
  void compute_timing_statistics()
  {
    for (auto &itImplTimer : m_implTimers)
      itImplTimer.compute_timing_stats();
  }

  void sort_timings(std::string type = "mean");
  void print_statistics_json(std::ostream&);
  void print_details(std::ostream&);
  void print_csv_header(std::ostream&);
  void print_details_csv(std::ostream&);
  void print_summary(std::ostream&);

  std::vector<double> fastest_times()
  {
    std::vector<double> minTimes;
    for (auto timer : m_implTimers) {
      minTimes.push_back(timer.minTime);
    }
    return minTimes;
  }

  std::vector<double> slowest_times()
  {
    std::vector<double> maxTimes;
    for (auto timer : m_implTimers) {
      maxTimes.push_back(timer.maxTime);
    }
    return maxTimes;
  }
};

/* Constructor's implementation */
template <typename IN_TYPE, typename OUT_TYPE>
inline ImplementationComparator<IN_TYPE,OUT_TYPE>::ImplementationComparator(
    ImplContainer<IN_TYPE,OUT_TYPE> &inImpl, int nEvals, unsigned int seed, std::unique_ptr<RngInterface<IN_TYPE>> inRng) :
  m_implementations(std::move(inImpl)), mp_sampler(std::move(inRng)), m_nEvals(nEvals)
{
  /*
     Allocate enough timer containers
  */
  m_numberOfImplementations = m_implementations.size();
  m_evaluationTimers.reserve(m_numberOfImplementations);
  m_evalHolder.reset(new OUT_TYPE[m_nEvals]);

  /*
    Initialize the timer structs with pointers to implementations
  */
  for (auto & itImpl : m_implementations) {
    m_implTimers.emplace_back(ImplTimer<IN_TYPE,OUT_TYPE>(itImpl.get()));
  }

  /*
    Ensure all implementations are using the same min/max range
  */
  auto itImpl = m_implementations.begin();
  m_minArg = (*itImpl)->min_arg(); // set to vals in first impl
  m_maxArg = (*itImpl)->max_arg();

  ++itImpl; // assert the rest of the impls are the same
  for (; (itImpl != m_implementations.end()); ++itImpl) {
    assert( abs( (**(itImpl)).min_arg() - m_minArg) <
  	    (std::numeric_limits<IN_TYPE>::epsilon())   );
    assert( abs( (*itImpl)->max_arg() - m_maxArg) <
  	    (std::numeric_limits<IN_TYPE>::epsilon())   );
  }

  /*
    Prepare to generate random points in the table interval to evaluate.
    If an RngInterface was not provided, set mp_sampler to a uniform real
    distribution on the table's endpoints
  */
  if(mp_sampler == nullptr)
    mp_sampler = std::unique_ptr<StdRng<IN_TYPE>>(new StdRng<IN_TYPE>(m_minArg, m_maxArg));
  mp_sampler->init(seed);
  mp_randomEvaluations = std::unique_ptr<IN_TYPE[]>(new IN_TYPE[m_nEvals]);
}

/* sort the vector of timings based on the min, max, or mean times */
template <typename IN_TYPE, typename OUT_TYPE>
inline void ImplementationComparator<IN_TYPE,OUT_TYPE>::sort_timings(std::string type)
{
  // default sort by mean time
  if (!type.compare("mean")) {
    sort( m_implTimers.begin(), m_implTimers.end(),
	       [](const ImplTimer<IN_TYPE,OUT_TYPE> &a, const ImplTimer<IN_TYPE,OUT_TYPE> &b)
	       { return (a.meanTime < b.meanTime); } );
  } else if (!type.compare("min")) {   // or sort by minimum (ie. best cases)
    sort( m_implTimers.begin(), m_implTimers.end(),
	       [](const ImplTimer<IN_TYPE,OUT_TYPE> &a, const ImplTimer<IN_TYPE,OUT_TYPE> &b)
	       { return (a.minTime < b.minTime); } );
  } else if (!type.compare("max")) {   // or sort by maximum time (ie. worst cases)
    sort( m_implTimers.begin(), m_implTimers.end(),
	       [](const ImplTimer<IN_TYPE,OUT_TYPE> &a, const ImplTimer<IN_TYPE,OUT_TYPE> &b)
	       { return (a.maxTime < b.maxTime); } );
  }
}

/* Implementation of functions that print to an ostream */
template <typename IN_TYPE, typename OUT_TYPE>
inline void ImplementationComparator<IN_TYPE,OUT_TYPE>::print_statistics_json(std::ostream &out)
{
  /* add implementations to vector */
  nlohmann::json jsonStats;

  jsonStats["_comment"] = "Timing data for implementations.";
  jsonStats["nEvals"] = m_nEvals;
  jsonStats["nTrials"] = (m_implTimers[0].evaluationTimes).size();

  for (auto itImplTimer : m_implTimers) {
    jsonStats[(itImplTimer.impl)->name()]["min"]  = itImplTimer.minTime;
    jsonStats[(itImplTimer.impl)->name()]["max"]  = itImplTimer.maxTime;
    jsonStats[(itImplTimer.impl)->name()]["mean"] = itImplTimer.meanTime;
    jsonStats[(itImplTimer.impl)->name()]["raw"]  = itImplTimer.evaluationTimes;
  }
  out << jsonStats.dump(2) << std::endl;
}

template <typename IN_TYPE, typename OUT_TYPE>
inline void ImplementationComparator<IN_TYPE,OUT_TYPE>::print_details(std::ostream &out)
{
  /*
     Print details of all timings
  */
  // temporary vars to hopefully get the template values
  IN_TYPE temp_in_type;
  OUT_TYPE temp_out_type;

  out << "----------------------------------------------------------------------------\n";
  out << "Table domain and range: "
      << typeid(temp_in_type).name()
      << " -> "
      << typeid(temp_out_type).name()
      << std::endl;
  out << "Number of trials performed: "
	    << m_implTimers[0].evaluationTimes.size()
	    << "\n";
  out << "Number of evaluations used: "
	    << m_nEvals
	    << std::endl;

  for (auto itImplTimer : m_implTimers) {
    out << "----------------------------------------------------------------------------\n";
    out << "| ";
    (itImplTimer.impl)->print_details(out);
    out << "| Memory usage (B): " << (itImplTimer.impl)->size() << "\n";
    out << "| Run times (seconds / " << m_nEvals << " evals):" << "\n";

    for ( auto itTimer : (itImplTimer.evaluationTimes) ) {
      out << "|    " << itTimer << "\n";
    }
    out << std::endl;

  }
  out << "----------------------------------------------------------------------------\n";
}

template <typename IN_TYPE, typename OUT_TYPE>
inline void ImplementationComparator<IN_TYPE,OUT_TYPE>::print_csv_header(std::ostream &out)
{
  /* Print header */
  for (auto itImplTimer : m_implTimers) {
    out << (itImplTimer.impl)->name() << " ";
  }
  out << "\n";
}

template <typename IN_TYPE, typename OUT_TYPE>
inline void ImplementationComparator<IN_TYPE,OUT_TYPE>::print_details_csv(std::ostream &out)
{
  print_csv_header(out);

  /* Print all times, row by row */
  int numTimes = (m_implTimers[0].evaluationTimes).size();
  for (int i = 0; i<numTimes; ++i) {
    for (auto itImplTimer : m_implTimers) {
      out << itImplTimer.evaluationTimes[i] << " ";
    }
    out << "\n";
  }
}

template <typename IN_TYPE, typename OUT_TYPE>
inline void ImplementationComparator<IN_TYPE,OUT_TYPE>::print_summary(std::ostream &out)
{
  /*
    Print summary of the timing statistics
  */
  out << "----------------------------------------------------------------------------\n";
  out << "Number of trials performed: "
	    << m_implTimers[0].evaluationTimes.size()
	    << std::endl;
  out << "Number of evaluations used: "
	    << m_nEvals
	    << std::endl;

  for (auto itImplTimer : m_implTimers) {
    out << "----------------------------------------------------------------------------\n";
    out << "| Implementation:   "; (itImplTimer.impl)->print_details(out);
    out << "\n| Memory usage (B): " << (itImplTimer.impl)->size() << "\n";
    out << "| Timings:          "; itImplTimer.print_timing_stats(out);
  }
  out << "----------------------------------------------------------------------------\n";
}

} // namespace func

#ifndef LUTDE_HPP
#define LUTDE_HPP
/*
  Class for comparing variations of EvaluationImplementations

  - takes ownership of the vector of implementations passed to it
*/
#pragma once
#include "Timer.hpp"
#include "RngInterface.hpp"
#include "EvaluationImplementation.hpp"

#include <memory>
#include <vector>
#include <random>
#include <algorithm>

/*
  Data types and containers used in Lutde
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
  ImplTimer(ImplType<IN_TYPE,OUT_TYPE> *inImpl) : impl(inImpl){};
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
    out << "Min " << minTime
      << " Max " << maxTime
      << " Mean " << meanTime
      << "\n";
  }
};

/* ------------------------------------------------------------------------ */
template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE>
class ImplementationComparator
{
private:

  std::vector<ImplTimer<IN_TYPE,OUT_TYPE>> m_implTimers;
  ImplContainer<IN_TYPE,OUT_TYPE>          m_implementations;
  unsigned                                 m_numberOfImplementations;

  std::vector<TimeContainer>      m_evaluationTimers;

  IN_TYPE                 m_minArg,m_maxArg;

  std::unique_ptr<OUT_TYPE[]> m_evalHolder;

  /*
    RNG for evaluations
    - By default uses a std::uniform_real_distribution<IN_TYPE>
      with the std::mt19937 variant of the std::mersenne_twister_engine
  */
  RngInterface<IN_TYPE> *mp_sampler;
  IN_TYPE               *mp_randomEvaluations;
  int                    m_nEvals;

  struct TimingStatistics
  {
    double minTime, maxTime, meanTime;
  } m_timingStatistics;

  void draw_new_sample_points();

  void run_all_single();

public:

  ImplementationComparator(ImplContainer<IN_TYPE,OUT_TYPE> &inImpl, int nEvals = 100000,
      unsigned int seed = 2017, RngInterface<IN_TYPE> *inRng = NULL);
  ~ImplementationComparator();

  void run_timings(int nRuns = 1);
  const IN_TYPE* sample_points();
  void compute_timing_statistics();
  void print_statistics_json(std::ostream&);
  void sort_timings(std::string type = "mean");
  void print_details(std::ostream&);
  void print_csv_header(std::ostream&);
  void print_details_csv(std::ostream&);
  void print_summary(std::ostream&);
  std::vector<double> fastest_times();
  std::vector<double> slowest_times();

};

#endif  // LUTDE_HPP

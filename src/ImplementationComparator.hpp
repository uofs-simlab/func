#ifndef LUTDE_HPP
#define LUTDE_HPP
/*
  Class for comparing variations of EvaluationImplementations

  - takes ownership of the vector of implementations passed to it
*/
#pragma once
#include "Timer.hpp"
#include "EvaluationImplementation.hpp"

#include <memory>
#include <vector>
#include <random>
#include <algorithm>

/*
  Data types and containers used in Lutde
*/
typedef std::vector<double> TimeContainer;

typedef double table_type;
typedef EvaluationImplementation ImplType;
typedef std::vector< std::unique_ptr<ImplType> > ImplContainer;

/*
  ImplTimer struct attaches additional data for timing an implementation
  to the implementation
*/
struct ImplTimer
{
  /* Note the the ImplType can NOT be a reference here, because we
     want to be able to sort a container of these ImplTimers. sort
     requires operator=, and classes that have non-static reference
     members cannot implement this */
  ImplType *impl;
  TimeContainer evaluationTimes;
  double maxTime, minTime, meanTime;
  ImplTimer(ImplType *inImpl) : impl(inImpl){};
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
class ImplementationComparator
{
private:

  ImplContainer              m_implementations;
  unsigned                   m_numberOfImplementations;
  std::vector<ImplTimer>     m_implTimers;

  std::vector<TimeContainer> m_evaluationTimers;

  table_type                 m_minArg,m_maxArg;

  std::unique_ptr<table_type[]> m_evalHolder;

  /*
    RNG for evaluations
    - Mersenne Twister
  */
  int                                     m_seed;
  std::mt19937                           *mp_generator;
  std::uniform_real_distribution<double> *mp_distribution;
  double                                 *mp_randomEvaluations;
  int                                     m_nEvals;

  struct TimingStatistics
  {
    double minTime, maxTime, meanTime;
  } m_timingStatistics;

  void draw_new_sample_points();

  void run_all_single();

public:

  ImplementationComparator(ImplContainer &inImpl, int nEvals = 100000, int seed = 2017);
  ~ImplementationComparator();

  void run_timings(int nRuns = 1);
  const double* sample_points();
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

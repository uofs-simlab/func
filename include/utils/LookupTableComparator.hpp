#pragma once
#include "Timer.hpp"
#include "RngInterface.hpp"
#include "StdRng.hpp"
#include "LookupTable.hpp"
#include "json.hpp"

#include <memory>
#include <vector>
#include <algorithm>
#include <limits> // epsilon()
#include <cassert>
#include <iostream>
#include <typeinfo> // typeid

namespace func {

template <typename TIN, typename TOUT>
using ImplContainer = std::vector<std::unique_ptr<LookupTable<TIN,TOUT>>>;

/**
 \brief specify how a LookupTableComparator should sort its results.
*/
enum class Sorter {NONE, BEST, MEAN, WORST};

/**
 \brief Helper struct: takes a LookupTable and attaches a set of timings to it
*/
template <typename TIN, typename TOUT>
struct ImplTimer
{
  /* Note that the LookupTable can NOT be a reference here, because we
     want to be able to sort a container of these ImplTimers. sort
     requires operator=, and classes that have non-static reference
     members cannot implement this */
  LookupTable<TIN,TOUT> *impl;
  std::vector<double> evaluationTimes;
  double maxTime, minTime, meanTime;

  ImplTimer(LookupTable<TIN,TOUT> *inImpl) : impl(inImpl), maxTime(0), minTime(0), meanTime(0) {};

  void append_runtime(double time){ evaluationTimes.push_back(time); };
  void compute_timing_stats() {
    /* get min and max execution times */
    auto extremeTimes = minmax_element( evaluationTimes.begin(), evaluationTimes.end() );
    minTime = *(extremeTimes.first);
    maxTime = *(extremeTimes.second);
    /* mean execution time */
    auto sum = std::accumulate(evaluationTimes.begin(), evaluationTimes.end(), 0.0);
    meanTime =  sum / evaluationTimes.size();
  }
  void print_timing_stats(std::ostream& out)
  {
    out << "Min " << minTime << "s"
        << " Max " << maxTime << "s"
        << " Mean " << meanTime << "s"
        << "\n";
  }
};

/**
  \brief Compare the average time taken to call the `operator()` of any
  LookupTable implementation.

  \ingroup Utils

  For example usage, see any file in the examples directory.

  \note This class takes ownership of the vector of LUT implementations it is
  constructed with
  \note Points are randomly sampled, and by default uses a
   `std::uniform_real_distribution<TIN>` with the `std::mt19937` variant of the
   `std::mersenne_twister_engine`. This can be changed by passing in a different `StdRng`

  \todo `LookupTableComparator`'s constructor should accept any callable type
*/
template <typename TIN, typename TOUT = TIN>
class LookupTableComparator
{
private:

  ImplContainer<TIN,TOUT>          m_implementations;
  std::vector<ImplTimer<TIN,TOUT>> m_implTimers;
  unsigned int                     m_numberOfImplementations;

  std::vector<std::vector<double>> m_evaluationTimers;

  TIN m_minArg,m_maxArg;

  std::unique_ptr<TOUT[]> m_evalHolder;

  /* RNG for evaluations */
  std::unique_ptr<RngInterface<TIN>> mp_sampler;
  std::unique_ptr<TIN[]>             mp_randomEvaluations;
  int                                m_nEvals;

  //struct TimingStatistics
  //{
  //  double minTime, maxTime, meanTime;
  //} m_timingStatistics;

  /** Fill mp_randomEvaluations with random points to be evaluated */
  void draw_new_sample_points() {
    for (int ii=0;ii<m_nEvals;++ii)
      mp_randomEvaluations[ii] = mp_sampler->get_point();
  }

  /** Time implementation evaluations */
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

  /** Prepare to run several timings for each LUT in the vector inImpl */
  LookupTableComparator(ImplContainer<TIN,TOUT> &inImpl, TIN minArg, TIN maxArg, unsigned int nEvals = 100000,
      unsigned int seed = 2017, std::unique_ptr<RngInterface<TIN>> inRng = nullptr);
  ~LookupTableComparator(){}

  /** Run timings with different set of random arguments */
  void run_timings(int nRuns = 1)
  {
    for (int ii=0;ii<nRuns;++ii) {
      draw_new_sample_points();
      run_all_single();
    }
  }

  /** Compute fastest and slowest times */
  void compute_statistics()
  {
    for (auto &itImplTimer : m_implTimers)
      itImplTimer.compute_timing_stats();
  }

  /** Sort the vector of implementations (m_implTimers) based on their max Sorter::WORST, mean Sorter::MEAN, or min Sorter::BEST times */
  void sort_timings(Sorter type = Sorter::MEAN);

  /** Print out the computed statistics for each LookupTable (no raw timings are displayed) */
  void print_summary(std::ostream&);

  /** Print out the raw timings for each LookupTable */
  void print_json(std::ostream&);
  void print_csv_header(std::ostream&);
  /** Output a space separated listing of timing results. Does not print a final newline if
   * type != Sorter::NONE which is helpful for plotting timing results with external programs (e.g. Python) */
  void print_csv(std::ostream&, Sorter type = Sorter::NONE);

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
template <typename TIN, typename TOUT>
inline LookupTableComparator<TIN,TOUT>::LookupTableComparator(
    ImplContainer<TIN,TOUT> &inImpl, TIN minArg, TIN maxArg, unsigned int nEvals, unsigned int seed, std::unique_ptr<RngInterface<TIN>> inRng) :
  m_implementations(std::move(inImpl)), mp_sampler(std::move(inRng)), m_nEvals(nEvals)
{
  /*
     Allocate enough timer containers
  */
  m_numberOfImplementations = m_implementations.size();
  m_evaluationTimers.reserve(m_numberOfImplementations);
  m_evalHolder.reset(new TOUT[m_nEvals]);

  /*
    Initialize the timer structs with pointers to implementations
  */
  for (auto & itImpl : m_implementations) {
    m_implTimers.emplace_back(ImplTimer<TIN,TOUT>(itImpl.get()));
  }

  // TODO alert the user if [minArg,maxArg] are outside the domain of any LUT
  m_minArg = minArg;
  m_maxArg = maxArg;

  /*
    Prepare to generate random points in the table interval to evaluate.
    If an RngInterface was not provided, set mp_sampler to a uniform real
    distribution on the table's endpoints
  */
  if(mp_sampler == nullptr)
    mp_sampler = std::unique_ptr<StdRng<TIN>>(new StdRng<TIN>(static_cast<double>(m_minArg), static_cast<double>(m_maxArg)));
  mp_sampler->init(seed);
  mp_randomEvaluations = std::unique_ptr<TIN[]>(new TIN[m_nEvals]);
}

template <typename TIN, typename TOUT>
inline void LookupTableComparator<TIN,TOUT>::sort_timings(Sorter type)
{
  switch(type){
  case Sorter::BEST: // sort by minimum (ie. best case performance)
  {
    sort(m_implTimers.begin(), m_implTimers.end(),
        [](const ImplTimer<TIN,TOUT> &a, const ImplTimer<TIN,TOUT> &b)
        { return (a.minTime < b.minTime); } );
    break;
  }
  case Sorter::MEAN: // default sort by mean time
  {
    sort(m_implTimers.begin(), m_implTimers.end(),
        [](const ImplTimer<TIN,TOUT> &a, const ImplTimer<TIN,TOUT> &b)
        { return (a.meanTime < b.meanTime); } );
    break;
  }
  case Sorter::WORST: // sort by maximum time (ie. worst case performance)
  {
    sort(m_implTimers.begin(), m_implTimers.end(),
        [](const ImplTimer<TIN,TOUT> &a, const ImplTimer<TIN,TOUT> &b)
        { return (a.maxTime < b.maxTime); } );
    break;
  }
  default: { throw std::logic_error("func::LookupTableComparator::sort_timings given an invalid Sorter."
               " Options are Sorter::BEST, Sorter::MEAN, or Sorter::WORST."); }
  }
}

/* -- Implementation of functions that print to an ostream -- */

template <typename TIN, typename TOUT>
inline void LookupTableComparator<TIN,TOUT>::print_json(std::ostream &out)
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

template <typename TIN, typename TOUT>
inline void LookupTableComparator<TIN,TOUT>::print_csv_header(std::ostream &out)
{
  /* Print header */
  for (auto itImplTimer : m_implTimers) {
    out << (itImplTimer.impl)->name() << " ";
  }
  out << "\n";
}

template <typename TIN, typename TOUT>
inline void LookupTableComparator<TIN,TOUT>::print_csv(std::ostream &out, Sorter type) {
  switch(type){
  case Sorter::BEST: {
    for (auto itImplTimer : m_implTimers)
      out << std::scientific << itImplTimer.minTime << "s ";
    break;
  }
  case Sorter::MEAN: {
    for (auto itImplTimer : m_implTimers)
      out << std::scientific << itImplTimer.meanTime << "s ";
    break;
  }
  case Sorter::WORST: {
    for (auto itImplTimer : m_implTimers)
      out << std::scientific << itImplTimer.maxTime << "s ";
    break;
  }
  default: {
    /* Print all times, row by row */
    int numTimes = (m_implTimers[0].evaluationTimes).size();
    for (int i = 0; i<numTimes; ++i) {
      for (auto itImplTimer : m_implTimers) {
        out << std::scientific << itImplTimer.evaluationTimes[i] << " ";
      }
      out << "\n";
    }
  }
  }
}

/* print basic info about a LookupTable */
template <typename TIN, typename TOUT = TIN>
std::ostream& operator<<(std::ostream& out, const LookupTableComparator<TIN,TOUT>& comparator){
  /*
    Print summary of the timing statistics
  */
  comparator.print_summary(out);
  return out;
}


template <typename TIN, typename TOUT>
inline void LookupTableComparator<TIN,TOUT>::print_summary(std::ostream &out)
{
  /* Print details of all timings
     temp_tin and temp_tout are needed to print template values (but the string is probably mangled).
     There are more legible solutions but this is good enough until anyone actually cares */
  TIN temp_tin;
  TOUT temp_tout;

  out << "----------------------------------------------------------------------------\n";
  out << "Table input and output types: "
      << typeid(temp_tin).name() << " -> " << typeid(temp_tout).name() << "\n";
  out << "Number of trials performed: "
	    << m_implTimers[0].evaluationTimes.size()
	    << "\n";
  out << "Number of evaluations used: "
	    << m_nEvals
	    << std::endl;

  for (auto itImplTimer : m_implTimers) {
    out << "----------------------------------------------------------------------------\n";
    out << "| LookupTable:      " << (*itImplTimer.impl) << "\n";
    out << "| Memory usage (B): " << (itImplTimer.impl)->size() << "\n";
    out << "| Timings:          "; itImplTimer.print_timing_stats(out);
  }
  out << "----------------------------------------------------------------------------" << std::endl;
}

} // namespace func

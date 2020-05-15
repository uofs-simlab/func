/*
  Implementation of ImplementationComparator.
*/
#include "ImplementationComparator.hpp"
#include "RngInterface.hpp"
#include "StdRng.hpp"
#include "json.hpp"

#include <iostream>
#include <limits>
#include <cassert>

ImplementationComparator::ImplementationComparator(ImplContainer &inImpl, int nEvals, unsigned int seed, RngInterface<> *inRng) : 
  m_implementations(std::move(inImpl)), m_nEvals(nEvals)
{
  /*
     Allocate enough timer containers
  */
  m_numberOfImplementations = m_implementations.size();
  m_evaluationTimers.reserve(m_numberOfImplementations);
  m_evalHolder.reset(new table_type[m_nEvals]);

  /*
    Initialize the timer structs with pointers to implementations
  */
  for (auto & itImpl : m_implementations) {
    m_implTimers.push_back(ImplTimer(itImpl.get()));
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
  	    (std::numeric_limits<table_type>::epsilon())   );
    assert( abs( (*itImpl)->max_arg() - m_maxArg) <
  	    (std::numeric_limits<table_type>::epsilon())   );
  }

  /*
    Generate random points in the table interval to evaluate
  */
  mp_sampler=inRng;
  // if inRng was not provided, set mp_sampler to a uniform distribution on the table's endpoints
  if(mp_sampler==nullptr)
    mp_sampler = new StdRng<std::uniform_real_distribution<double>>
      (new std::uniform_real_distribution<double>(m_minArg, m_maxArg));
  mp_sampler->init(seed);
  mp_randomEvaluations = new table_type[m_nEvals];
}

ImplementationComparator::~ImplementationComparator()
{
  delete[] mp_randomEvaluations;
}

void ImplementationComparator::draw_new_sample_points()
{
  /*
    Regenerate evaluation points.
  */
  for (int ii=0;ii<m_nEvals;++ii)
    mp_randomEvaluations[ii] = mp_sampler->getPt();
}

const double* ImplementationComparator:: sample_points()
{
  /*
     Obtain a const reference to the current array of evaluation points
  */
  return mp_randomEvaluations;
}
void ImplementationComparator::run_all_single()
{
  /*
    Time implementation evaluations
   */
  for (auto &itImplTimer : m_implTimers) {
    Timer evalTimer;
    for (int ii=0;ii<m_nEvals;++ii) {
       m_evalHolder[ii] = (*(itImplTimer.impl))(mp_randomEvaluations[ii]);
    }
    evalTimer.stop();
    itImplTimer.append_runtime(evalTimer.duration());
  }
}

void ImplementationComparator::run_timings(int nRuns)
{
  /*
     Run timings with different set of random arguments
  */
  for (int ii=0;ii<nRuns;++ii) {
    draw_new_sample_points();
    run_all_single();
  }
}

void ImplementationComparator::compute_timing_statistics()
{
  /*
    Compute fastest and slowest
  */
  for (auto &itImplTimer : m_implTimers) {
    itImplTimer.compute_timing_stats();
  }
}
void ImplementationComparator::print_statistics_json(std::ostream &out)
{
  /* add implementations to vector */
  using nlohmann::json;
  json jsonStats;

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


void ImplementationComparator::sort_timings(std::string type)
{
  // default sort by mean time
  if (!type.compare("mean")) {
    sort( m_implTimers.begin(), m_implTimers.end(),
	       [](const ImplTimer &a, const ImplTimer &b)
	       { return (a.meanTime < b.meanTime); } );
  } else if (!type.compare("min")) {   // or sort by minimum (ie. best cases)
    sort( m_implTimers.begin(), m_implTimers.end(),
	       [](const ImplTimer &a, const ImplTimer &b)
	       { return (a.minTime < b.minTime); } );
  } else if (!type.compare("max")) {   // or sort by maximum time (ie. worst cases)
    sort( m_implTimers.begin(), m_implTimers.end(),
	       [](const ImplTimer &a, const ImplTimer &b)
	       { return (a.maxTime < b.maxTime); } );
  }
}

void ImplementationComparator::print_details(std::ostream &out)
{
  /*
     Print details of all timings
  */
  out << "----------------------------------------------------------------------------\n";
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
void ImplementationComparator::print_csv_header(std::ostream &out)
{
  /* Print header */
  for (auto itImplTimer : m_implTimers) {
    out << (itImplTimer.impl)->name() << " ";
  }
  out << "\n";
}
void ImplementationComparator::print_details_csv(std::ostream &out)
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

void ImplementationComparator::print_summary(std::ostream &out)
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

std::vector<double> ImplementationComparator::fastest_times()
{
  std::vector<double> minTimes;
  for (auto timer : m_implTimers) {
    minTimes.push_back(timer.minTime);
  }
  return minTimes;
}
std::vector<double> ImplementationComparator::slowest_times()
{
  std::vector<double> maxTimes;
  for (auto timer : m_implTimers) {
    maxTimes.push_back(timer.maxTime);
  }
  return maxTimes;
}

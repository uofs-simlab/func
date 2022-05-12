/*
  Creating implementations big enough to guarantee they aren't stored
  in cache.

  Usage:
      ./experiment_best_worst <tableMin> <tableMax> <tableTol> <nExperiments> <nEvals> <seed>
 */

// #include "sin.hpp"
// #include "func.hpp"
#include "ZeroFunction.hpp"
// #include "chaste_log_function.hpp"

#include "ImplementationComparator.hpp"
#include "UniformLookupTableGenerator.hpp"
#include "UniformTables.hpp"
#include "DirectEvaluation.hpp"

#include <iostream>
#include <fstream>
#include <cassert>

void print_usage()
{
  std::cout << "Usage:" << std::endl
	    << "    ./experiment_best_worst <tableSizeFactor> <nExperiments> <nEvals> <seed>"
	    << std::endl;
}

//---------------------------------------------------------------------------->>
int main(int argc, char* argv[])
{

  using namespace std;

  if (argc < 5) {
      print_usage();
      exit(0);
  }

  double    tableSizeFactor = std::stof(argv[1]);
  int    nExperiments = std::stod(argv[2]);
  int    nEvals       = std::stoi(argv[3]);
  int    seed         = std::stoi(argv[4]);

  FunctionContainer func_container;
  func_container.double_func = new ZeroFunction<double>;
  func_container.fvar1_func  = new ZeroFunction<fvar1>;
  func_container.fvar2_func  = new ZeroFunction<fvar2>;
  func_container.fvar3_func  = new ZeroFunction<fvar3>;

  double stepSize;

  /* Which LUT implementations to use */
  std::vector<std::string> implNames {"UniformLinearInterpolationTable",
      "UniformLinearPrecomputedInterpolationTable",
      "UniformQuadraticPrecomputedInterpolationTable",
      "UniformCubicPrecomputedInterpolationTable",
      "UniformLinearTaylorTable",
      "UniformQuadraticTaylorTable",
      "UniformCubicTaylorTable"};

  /* get cache sizes running something like `lscpu | grep cache` */
  const unsigned long cacheSize = (3072u+256u+32u)*1024u;
  const unsigned long ramSize = 8192lu*1024lu*1024lu;
  const unsigned desiredTableSize = (unsigned)(tableSizeFactor*ramSize/implNames.size());
  // const unsigned desiredTableSize = (unsigned)tableSizeFactor*cacheSize;

  double percentRam = 100*implNames.size()*(double)desiredTableSize/(double)ramSize;

  assert(percentRam < 75);

  std::cout << "\n# impls using ~ " << percentRam <<"% of RAM\n";
  cout << "# Function:  " << FUNCNAME << endl << endl;

  UniformLookupTableGenerator gen(&func_container, 0, 1);

  /* Fill in the implementations */
  std::vector<unique_ptr<EvaluationImplementation>> impls;


  /* Run the best case */
  std::cout << "# Running Best case (small tables)" << std::endl;

  for (auto itName : implNames) {
    impls.emplace_back(gen.generate_by_step(itName,1.0));
  }

  /* Run comparator */
  ImplementationComparator implCompare_best(impls, nEvals, seed);
  implCompare_best.run_timings(nExperiments);

  /* Summarize the results */
  implCompare_best.compute_timing_statistics();
  std::ofstream jsonfs;
  jsonfs.open("best_case.json");
  implCompare_best.print_statistics_json(jsonfs);
  jsonfs.close();

  implCompare_best.sort_timings("min");
  implCompare_best.print_summary(std::cout);


  impls.clear();


  /* Run the worst case */
  std::cout << "# Running Worst case (large tables)" << std::endl;

  /* add implementations to vector */
  for (auto itName : implNames) {
    impls.emplace_back(gen.generate_by_impl_size(itName,desiredTableSize));
  }

  /* Run comparator */
  ImplementationComparator implCompare_worst(impls, nEvals, seed);
  implCompare_worst.run_timings(nExperiments);

  /* Summarize the results */
  implCompare_worst.compute_timing_statistics();
  jsonfs.open("worst_case.json");
  implCompare_worst.print_statistics_json(jsonfs);
  jsonfs.close();

  implCompare_worst.sort_timings("min");
  implCompare_worst.print_summary(std::cout);

  return 0;
}

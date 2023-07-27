/*
  Creating implementations big enough to guarantee they aren't stored
  in cache.

  Usage:
      ./experiment_best_worst <tableMin> <tableMax> <tableTol> <nExperiments> <nEvals> <seed>
 */

#include "func.hpp"
#include "ZeroFunction.hpp"

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
  using namespace func;

  if (argc < 5) {
      print_usage();
      exit(0);
  }

  double    tableSizeFactor = std::stof(argv[1]);
  int    nExperiments = std::stod(argv[2]);
  int    nEvals       = std::stoi(argv[3]);
  int    seed         = std::stoi(argv[4]);

  FunctionContainer<double> func_container {FUNC_SET_F(ZeroFunction,double)};

  // double stepSize;

  /* Which LUT implementations to use */
  std::vector<std::string> implNames {
    "UniformLinearRawInterpTable",
    "UniformInterpTable<1>",
    "UniformInterpTable<2>",
    "UniformInterpTable<3>",
    "UniformTaylorTable<1>",
    "UniformTaylorTable<2>",
    "UniformTaylorTable<3>",
  };

  /* get cache sizes running something like `lscpu | grep cache` */
  // const unsigned long cacheSize = (3072u+256u+32u)*1024u;
  const unsigned long ramSize = 8192lu*1024lu*1024lu;
  const unsigned desiredTableSize = (unsigned)(tableSizeFactor*ramSize/implNames.size());
  // const unsigned desiredTableSize = (unsigned)tableSizeFactor*cacheSize;

  double percentRam = 100*implNames.size()*(double)desiredTableSize/(double)ramSize;

  assert(percentRam < 75);

  std::cout << "\n# impls using ~ " << percentRam <<"% of RAM\n";
  cout << "# Function:  " << FUNCNAME << endl << endl;

  LookupTableGenerator<double> gen(func_container, 0, 1);

  /* Fill in the implementations */
  std::vector<unique_ptr<LookupTable<double>>> impls;


  /* Run the best case */
  std::cout << "# Running Best case (small tables)" << std::endl;

  for (auto itName : implNames) {
    impls.emplace_back(gen.generate_by_step(itName,1.0));
  }

  /* Run comparator */
  LookupTableComparator<double> implCompare_best(impls, 0.0, 1.0, nEvals, seed);
  implCompare_best.run_timings(nExperiments);

  /* Summarize the results */
  implCompare_best.compute_statistics();
  std::ofstream jsonfs;
  jsonfs.open("best_case.json");
  implCompare_best.print_json(jsonfs);
  jsonfs.close();

  implCompare_best.sort_timings(Sorter::BEST);
  implCompare_best.print_summary(std::cout);


  impls.clear();


  /* Run the worst case */
  std::cout << "# Running Worst case (large tables)" << std::endl;

  /* add implementations to vector */
  for (auto itName : implNames) {
    impls.emplace_back(gen.generate_by_impl_size(itName,desiredTableSize));
  }

  /* Run comparator */
  LookupTableComparator<double> implCompare_worst(impls, 0.0, 1.0, nEvals, seed);
  implCompare_worst.run_timings(nExperiments);

  /* Summarize the results */
  implCompare_worst.compute_statistics();
  jsonfs.open("worst_case.json");
  implCompare_worst.print_json(jsonfs);
  jsonfs.close();

  implCompare_worst.sort_timings(Sorter::BEST);
  implCompare_worst.print_summary(std::cout);

  return 0;
}

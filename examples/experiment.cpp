/*
  Main program for using ImplementationComparator to compare LookUpTable and Direct
  Evaluation performance.

  Usage:
      ./experiment <tableMin> <tableMax> <tableTol> <nExperiments> <nEvals> <seed>
 */

#include "func.hpp"

#include "chaste_log_function.hpp"

/*
  func.hpp contains the following:
  - func is used to pass into LookupTable for initialization (generality)
  - FUNC(x) macro can be used to avoid as much overhead as possible in direct
    evaluations, but loses the direct
*/

#include <iostream>
#define TYPE double

void print_usage()
{
  std::cout << "Usage:" << std::endl
	    << "    ./experiment <tableMin> <tableMax> <tableTol> <nExperiments> <nEvals> <seed>"
	    << std::endl;
}

//---------------------------------------------------------------------------->>
int main(int argc, char* argv[])
{

  using namespace std;
  using namespace func;

  if (argc < 7) {
      print_usage();
      exit(0);
  }

  double tableMin     = std::stod(argv[1]);
  double tableMax     = std::stod(argv[2]);
  double tableTol     = std::stod(argv[3]);
  int    nExperiments = std::stod(argv[4]);
  int    nEvals       = std::stoi(argv[5]);
  unsigned int seed   = std::stoi(argv[6]);

  FunctionContainer<TYPE> func_container{FUNC_SET_F(MyFunction,TYPE)};

  LookupTableFactory<TYPE> factory;

  /* Check which implementations are available */
  std::cout << "# Registered tables: \n#  ";
  for (auto it : factory.get_registered_keys() ) {
    std::cout << it << "\n#  ";
  }
  std::cout << "\n";

  /* Fill in the implementations */
  std::vector<unique_ptr<LookupTable<TYPE>>> impls;

  /* Which LUT implementations to use */
  std::vector<std::string> implNames {
    //"UniformLinearRawInterpTable",
    //"UniformEqSpaceInterpTable<1>",
    //"UniformEqSpaceInterpTable<2>",
    //"UniformEqSpaceInterpTable<3>",
    //"UniformTaylorTable<1>",
    //"UniformTaylorTable<2>",
    //"UniformTaylorTable<3>",
    "UniformChebyInterpTable<1>",
    "UniformChebyInterpTable<2>",
    "UniformChebyInterpTable<3>",
    "NonUniformChebyInterpTable<1>",
    "NonUniformChebyInterpTable<2>",
    "NonUniformChebyInterpTable<3>",
  };

  //std::vector<std::string> padeNames {
  //  "UniformPadeTable<1,1>",
  //  "UniformPadeTable<2,1>",
  //  "UniformPadeTable<3,1>",
  //  "UniformPadeTable<4,1>",
  //  "UniformPadeTable<5,1>",
  //  "UniformPadeTable<6,1>",
  //  "UniformPadeTable<2,2>",
  //  "UniformPadeTable<3,2>",
  //  "UniformPadeTable<4,2>",
  //  "UniformPadeTable<5,2>",
  //  "UniformPadeTable<3,3>",
  //  "UniformPadeTable<4,3>",
  //};

  LookupTableGenerator<TYPE> gen(func_container, tableMin, tableMax);

  impls.emplace_back(unique_ptr<LookupTable<TYPE>>(new DirectEvaluation<TYPE>(func_container,tableMin,tableMax)));
  for (auto itName : implNames) {
    cerr << "Building " << itName << " ..." << endl;
    impls.emplace_back(gen.generate_by_tol(itName,tableTol));
  }

  std::cout << "Running timings ..." << std::endl;

  ImplementationComparator<TYPE> implCompare(impls, tableMin, tableMax, nEvals, seed);
  implCompare.run_timings(nExperiments);

  /* Summarize the results */
  cout << "# Function:  " << FUNCNAME << std::endl;
  cout << "# Range:      (" << tableMin << "," << tableMax << ")" << std::endl;

  implCompare.compute_timing_statistics();
  implCompare.sort_timings(SortType::worst);
  implCompare.print_summary(std::cout);
}

/*
  Main program for using LookupTableComparator to compare LookUpTable and Direct
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

void print_usage(){
  std::cout << "Usage:" << std::endl
	    << "    ./experiment <tableMin> <tableMax> <tableTol> <nExperiments> <nEvals> <seed>"
	    << std::endl;
}

//---------------------------------------------------------------------------->>
int main(int argc, char* argv[]){
  using namespace std;
  using namespace func;

  if(argc < 7){
      print_usage();
      return 1;
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
  std::vector<std::string> uniformNames {
    "UniformChebyInterpTable<1>",
    "UniformChebyInterpTable<2>",
    "UniformChebyInterpTable<3>",
    "UniformChebyInterpTable<4>",
    "UniformChebyInterpTable<5>",
    "UniformChebyInterpTable<6>",
    "UniformChebyInterpTable<7>",
    //"UniformCubicHermiteTable",
    //"UniformEqSpaceInterpTable<1>",
    //"UniformEqSpaceInterpTable<2>",
    //"UniformEqSpaceInterpTable<3>",
    //"UniformLinearRawInterpTable",
    //"UniformTaylorTable<1>",
    //"UniformTaylorTable<2>",
    //"UniformTaylorTable<3>",
    //"UniformTaylorTable<4>",
    //"UniformTaylorTable<5>",
    //"UniformTaylorTable<6>",
    //"UniformTaylorTable<7>",
  };

  std::vector<std::string> nonuniformNames {
    "NonUniformChebyInterpTable<1>",
    "NonUniformChebyInterpTable<2>",
    "NonUniformChebyInterpTable<3>",
    "NonUniformChebyInterpTable<4>",
    "NonUniformChebyInterpTable<5>",
    //"NonUniformChebyInterpTable<6>",
    //"NonUniformChebyInterpTable<7>",
    //"NonUniformCubicHermiteTable",
    //"NonUniformEqSpaceInterpTable<1>",
    //"NonUniformEqSpaceInterpTable<2>",
    //"NonUniformEqSpaceInterpTable<3>",
    //"NonUniformTaylorTable<1>",
    //"NonUniformTaylorTable<2>",
    //"NonUniformTaylorTable<3>",
    //"NonUniformTaylorTable<4>",
    //"NonUniformTaylorTable<5>",
    //"NonUniformTaylorTable<6>",
    //"NonUniformTaylorTable<7>",
  };


  std::vector<std::string> padeNames {
    //"UniformPadeTable<1,1>",
    //"UniformPadeTable<2,1>",
    //"UniformPadeTable<3,1>",
    //"UniformPadeTable<4,1>",
    //"UniformPadeTable<5,1>",
    //"UniformPadeTable<6,1>",
    //"UniformPadeTable<2,2>",
    //"UniformPadeTable<3,2>",
    //"UniformPadeTable<4,2>",
    //"UniformPadeTable<5,2>",
    //"UniformPadeTable<3,3>",
    //"UniformPadeTable<4,3>",
  };

  LookupTableGenerator<TYPE> gen(func_container, tableMin, tableMax);

  impls.emplace_back(unique_ptr<LookupTable<TYPE>>(new DirectEvaluation<TYPE>(func_container,tableMin,tableMax)));
  for (auto itName : uniformNames) {
    std::cerr << "Building " << itName << " ..." << std::endl;
    impls.emplace_back(gen.generate_by_tol(itName,tableTol));
  }
  for (auto itName : nonuniformNames) {
    std::cerr << "Building " << itName << " ..." << std::endl;
    impls.emplace_back(gen.generate_by_tol(itName,tableTol));
  }
  for (auto itName : padeNames) {
    std::cerr << "Building " << itName << " ..." << std::endl;
    impls.emplace_back(gen.generate_by_tol(itName,tableTol));
  }


  std::cout << "Running timings ..." << std::endl;

  LookupTableComparator<TYPE> implCompare(impls, tableMin, tableMax, nEvals, seed);
  implCompare.run_timings(nExperiments);

  /* Summarize the results */
  cout << "# Function:  " << FUNCNAME << std::endl;
  cout << "# Domain:      (" << tableMin << "," << tableMax << ")" << std::endl;

  implCompare.compute_statistics();
  implCompare.sort_timings(Sorter::WORST);
  implCompare.print_summary(std::cout);
  return 0;
}

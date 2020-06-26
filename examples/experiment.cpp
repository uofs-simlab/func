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

  MyFunction<boost_fvar> func;
  MyFunction<double>     double_func;
  double stepSize;

  /* Check which implementations are available */
  std::cout << "# Registered uniform tables: \n#  ";
  for (auto it : UniformLookupTableFactory::get_registry_keys() ) {
    std::cout << it << "\n#  ";
  }
  std::cout << "\n";

  /* Fill in the implementations */
  std::vector<unique_ptr<EvaluationImplementation>> impls;

  /* Which LUT implementations to use */
  std::vector<std::string> implNames {
      //"UniformCubicHermiteTable",
      //"UniformLinearInterpolationTable",
      //"UniformLinearPrecomputedInterpolationTable",
      //"UniformQuadraticPrecomputedInterpolationTable",
      "UniformCubicPrecomputedInterpolationTable",
      //"UniformArmadilloPrecomputedInterpolationTable<4>",
      //"UniformArmadilloPrecomputedInterpolationTable<5>",
      //"UniformArmadilloPrecomputedInterpolationTable<6>",
      //"UniformArmadilloPrecomputedInterpolationTable<7>",
      //"UniformPadeTable<1,1>",
      //"UniformPadeTable<2,2>",
      //"UniformPadeTable<3,3>",
      //"UniformPadeTable<4,3>",
      "UniformLinearTaylorTable"
      //"UniformQuadraticTaylorTable",
      //"UniformCubicTaylorTable"
      };

  UniformLookupTableGenerator gen(&func, tableMin, tableMax);
  cout << "foo1" << endl;

  /* add implementations to vector */
  // unique_ptr<EvaluationImplementation> test = make_unique<DirectEvaluation>(&func,tableMin,tableMax);

  impls.emplace_back(unique_ptr<EvaluationImplementation>(new DirectEvaluation(&double_func,tableMin,tableMax)));
  cout << "foo4" << endl;
  for (auto itName : implNames) {
    impls.emplace_back(gen.generate_by_tol(itName,tableTol));
    cout << "foo6" << endl;
  }
  cout << "foo2" << endl;

  /* Run comparator */
  ImplementationComparator implCompare(impls, nEvals, seed);
  cout << "foo3" << endl;
  implCompare.run_timings(nExperiments);
  cout << "foo3" << endl;

  /* Summarize the results */
  cout << "# Function:  " << FUNCNAME << endl;
  cout << "# Range:      (" << tableMin << "," << tableMax << ")" << endl;

  implCompare.compute_timing_statistics();
  implCompare.sort_timings("max");
  implCompare.print_summary(std::cout);

  return 0;
}

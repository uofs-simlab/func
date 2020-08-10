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

  FunctionContainer<double> func_container{SET_F(MyFunction,double)};
  FunctionContainer<float> func_container_f{SET_F(MyFunction,float)};

  double stepSize;

  /* Check which implementations are available */
  std::cout << "# Registered uniform tables: \n#  ";
  for (auto it : UniformLookupTableFactory<double>::get_registry_keys() ) {
    std::cout << it << "\n#  ";
  }
  std::cout << "\n";

  /* Fill in the implementations */
  std::vector<unique_ptr<EvaluationImplementation<double>>> impls;
  std::vector<unique_ptr<EvaluationImplementation<float>>> implsf;

  /* Which LUT implementations to use */
  std::vector<std::string> implNames {
    //"UniformArmadilloPrecomputedInterpolationTable<4>",
    //"UniformArmadilloPrecomputedInterpolationTable<5>",
    //"UniformArmadilloPrecomputedInterpolationTable<6>",
    //"UniformArmadilloPrecomputedInterpolationTable<7>",
    //"UniformCubicHermiteTable",
    "UniformCubicPrecomputedInterpolationTable",
    //"UniformCubicTaylorTable",
    "UniformLinearInterpolationTable",
    //"UniformLinearPrecomputedInterpolationTable",
    //"UniformLinearTaylorTable",
    //"UniformQuadraticPrecomputedInterpolationTable",
    //"UniformQuadraticTaylorTable",
    "NonUniformLinearInterpolationTable",
    "NonUniformCubicPrecomputedInterpolationTable"
  };

  std::vector<std::string> padeNames {
    "UniformPadeTable<1,1>",
    "UniformPadeTable<2,1>",
    "UniformPadeTable<3,1>",
    "UniformPadeTable<4,1>",
    "UniformPadeTable<5,1>",
    "UniformPadeTable<6,1>",
    "UniformPadeTable<2,2>",
    "UniformPadeTable<3,2>",
    "UniformPadeTable<4,2>",
    "UniformPadeTable<5,2>",
    "UniformPadeTable<3,3>",
    "UniformPadeTable<4,3>"
  };

  UniformLookupTableGenerator<double> gen(&func_container, tableMin, tableMax);
  UniformLookupTableGenerator<float> genf(&func_container_f, tableMin, tableMax);

  impls.emplace_back(unique_ptr<EvaluationImplementation<double>>(new DirectEvaluation<double>(&func_container,tableMin,tableMax)));
  for (auto itName : implNames) {
    cout << "about to build " << itName << endl;
    impls.emplace_back(gen.generate_by_tol(itName,tableTol));
  }

  //implsf.emplace_back(unique_ptr<EvaluationImplementation<float>>(new DirectEvaluation<float>(&func_container_f,tableMin,tableMax)));
  //for(auto itName : implNames){
  //  implsf.emplace_back(genf.generate_by_tol(itName,tableTol));
  //}
  //for (auto itName : padeNames) {
  //  impls.emplace_back(gen.generate_by_tol(itName,tableTol));
  //}
  
  //add a composite table. The generator doesn't converge with mid = the root 
  //and also this is far too specific for this experiment.
  //double mid = exp(7.7/13.0287)+1;
  //UniformLookupTableGenerator gen1(&func_container, tableMin, mid);
  //UniformLookupTableGenerator gen2(&func_container, mid, tableMax);

  //impls.emplace_back(unique_ptr<EvaluationImplementation>(
  //      new CompositeLookupTable({gen1.generate_by_tol(implNames[3],tableTol*(mid-tableMin)/(tableMax-tableMin)),
  //        gen2.generate_by_tol(implNames[7],tableTol*(tableMax-mid)/(tableMax-tableMin))})));
 
  ImplementationComparator<double> implCompare(impls, nEvals, seed);
  implCompare.run_timings(nExperiments);

  //ImplementationComparator<float> implComparef(implsf, nEvals, seed);
  //implComparef.run_timings(nExperiments);

  /* Summarize the results */
  cout << "# Function:  " << FUNCNAME << endl;
  cout << "# Range:      (" << tableMin << "," << tableMax << ")" << endl;

  implCompare.compute_timing_statistics();
  implCompare.sort_timings("max");
  implCompare.print_summary(std::cout);

  //implComparef.compute_timing_statistics();
  //implComparef.sort_timings("max");
  //implComparef.print_summary(std::cout);

  return 0;
}

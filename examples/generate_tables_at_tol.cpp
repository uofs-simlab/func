#include "func.hpp"

#include "chaste_log_function.hpp"

#define MIN_ARG 0.001
#define MAX_ARG 30.0
#define TOL 1e-4

#include <iostream>
#include <vector>
#include <iomanip>

/*
  Simple program that uses the UniformLookupTableGenerator class to
  generate various table implementations at a desired tolerance
*/
int main()
{
  using namespace std;

  FunctionContainer<double> func_container{SET_F(MyFunction,double)};

  cout << "# Function: " << FUNCNAME << "\n";
  cout << "# Tol:      " << TOL << "\n";

  vector<unique_ptr<EvaluationImplementation<double>>> impls;

  /* Which LUT implementations to use */
  vector<string> implNames {"UniformLinearInterpolationTable",
      "UniformQuadraticPrecomputedInterpolationTable",
      "UniformCubicPrecomputedInterpolationTable",
      "UniformLinearTaylorTable",
      "UniformQuadraticTaylorTable",
      "UniformCubicTaylorTable"};


  UniformLookupTableGenerator<double> gen(&func_container, MIN_ARG, MAX_ARG);

  for (auto itName : implNames) {
    std::cout << "\nGenerating " << itName << ":" << std::endl;
    impls.emplace_back(gen.generate_by_tol(itName,TOL));
  }

  cout << "# Type, min_arg, max_arg, step_size, num_intervals\n";
  for (auto& impl : impls) {
    impl->print_details(cout);
    cout << "\n";
  }

  return 0;
}  // int main()

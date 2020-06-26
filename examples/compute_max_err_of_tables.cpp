#include "func.hpp"
#include "chaste_log_function.hpp"

#define MIN_ARG 0.001
#define MAX_ARG 30.0

#include <iostream>

/*
  Simple program that uses the UniformLookupTableGenerator class to
  compute errors for varying table step sizes
*/
int main()
{
  using namespace std;

  MyFunction<boost_fvar> func;

  /* Which implementations to use */
  std::vector<std::string> implNames {"UniformLinearInterpolationTable",
      "UniformQuadraticPrecomputedInterpolationTable",
      "UniformCubicPrecomputedInterpolationTable",
      "UniformLinearTaylorTable",
      "UniformQuadraticTaylorTable",
      "UniformCubicTaylorTable"};

  UniformLookupTableGenerator gen(&func, MIN_ARG, MAX_ARG);

  cout << "# Function: " << FUNCNAME << endl;
  cout << "# h ";
  for (auto itName : implNames) {
    std::cout << itName << " ";
    }

  cout << endl;
  for (double h = 1.0; h>0.01; h /=2) {

    std::cout << h;

    for (auto itName : implNames) {
      double err = gen.error_at_step_size(itName,h);
      std::cout << " " << err;
    }

    std::cout << std::endl;
  }

  return 0;
}  // int main()

#include "func.hpp"
#include "chaste_log_function.hpp"

#define MIN_ARG 0.001
#define MAX_ARG 30.0

#include <iostream>

/*
  Simple program that uses the LookupTableGenerator class to
  compute errors for varying table step sizes
*/
int main()
{
  using namespace std;
  using namespace func;

  FunctionContainer<double> func_container{FUNC_SET_F(MyFunction,double)};

  /* Which implementations to use */
  std::vector<std::string> implNames {"UniformLinearInterpolationTable",
      "UniformQuadraticPrecomputedInterpolationTable",
      "UniformCubicPrecomputedInterpolationTable",
      "UniformLinearTaylorTable",
      "UniformQuadraticTaylorTable",
      "UniformCubicTaylorTable"};

  LookupTableGenerator<double> gen(&func_container, MIN_ARG, MAX_ARG);

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

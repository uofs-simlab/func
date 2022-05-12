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

  FunctionContainer func_container;
  func_container.double_func = new MyFunction<double>;
  func_container.fvar1_func  = new MyFunction<fvar1>;
  func_container.fvar2_func  = new MyFunction<fvar2>;
  func_container.fvar3_func  = new MyFunction<fvar3>;

  /* Which implementations to use */
  std::vector<std::string> implNames {"UniformLinearInterpolationTable",
      "UniformQuadraticPrecomputedInterpolationTable",
      "UniformCubicPrecomputedInterpolationTable",
      "UniformLinearTaylorTable",
      "UniformQuadraticTaylorTable",
      "UniformCubicTaylorTable"};

  UniformLookupTableGenerator gen(&func_container, MIN_ARG, MAX_ARG);

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

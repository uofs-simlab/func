#include "func.hpp"

#include "chaste_log_function.hpp"
// #include "test_function.hpp"

#define MIN_ARG 0.5
#define MAX_ARG 4
#define STEP (MAX_ARG-MIN_ARG)/2

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

  cout << "# Function: " << FUNCNAME << endl;
  cout << "# h";
  cout << endl;

  LookupTableGenerator<double> gen(&func_container,MIN_ARG,MAX_ARG);

  // gen.plot_implementation_at_step_size("UniformLinearTaylorTable",STEP);
  // gen.plot_implementation_at_step_size("UniformQuadraticTaylorTable",STEP);
  gen.plot_implementation_at_step_size("UniformCubicTaylorTable",STEP);

  return 0;
}  // int main()

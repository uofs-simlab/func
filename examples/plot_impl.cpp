#include "func.hpp"

#include "test_function.hpp"

#define MIN_ARG 1.0
#define MAX_ARG 3.0
#define STEP (MAX_ARG-MIN_ARG)/2

#include <iostream>

/*
  Simple program that uses the UniformLookupTableGenerator class to
  compute errors for varying table step sizes
*/
int main()
{
  using namespace std;

  MyFunction func;

  cout << "# Function: " << FUNCNAME << endl;
  cout << "# h";
  cout << endl;

  // UniformLookupTableGenerator<UniformLinearInterpolationTable>
  //   gen(&func,MIN_ARG,MAX_ARG,0.0);

  UniformLookupTableGenerator gen(&func,MIN_ARG,MAX_ARG);

  gen.plot_implementation_at_step_size("UniformLinearInterpolationTable",STEP);

  return 0;
}  // int main()

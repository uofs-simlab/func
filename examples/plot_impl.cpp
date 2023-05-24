#include "func.hpp"

#include "chaste_log_function.hpp"
// #include "test_function.hpp"


#include <iostream>

void print_usage()
{
  std::cout << "Usage:\n"
	    << "    ./experiment <tableMin> <tableMax> <tableStep>"
	    << std::endl;
}

/*
  Simple program that uses the LookupTableGenerator class to
  compute errors for varying table step sizes
*/
int main(int argc, char* argv[])
{
  using namespace std;
  using namespace func;

  if (argc != 4) {
      print_usage();
      exit(0);
  }

  double tableMin  = std::stod(argv[1]);
  double tableMax  = std::stod(argv[2]);
  double tableStep = std::stod(argv[3]);

  FunctionContainer<double> func_container{FUNC_SET_F(MyFunction,double)};

  //std::cout << "# Function: " << FUNCNAME << "\n";
  //std::cout << "# h " << tableStep << "\n";
  //std::cout << std::endl;

  LookupTableGenerator<double> gen(func_container,tableMin,tableMax);
  gen.plot_implementation_at_step_size("UniformTaylorTable<1>",tableStep);
}

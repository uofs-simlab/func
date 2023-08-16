#include <iostream>
#include "func.hpp"
#include "chaste_log_function.hpp"

void print_usage(){
  std::cout << "Usage:\n"
    << "    ./experiment <tableKey> <tableMin> <tableMax> <tableStep>"
    << std::endl;
  std::cout << "Acceptable values of tableKey are any of the following:"
    << std::endl;

  /* Print every available implementation */
  func::LookupTableFactory<double> factory;
  std::cout << "# Registered tables: \n#  ";
  for (auto it : factory.get_registered_keys() ) {
    std::cout << it << "\n#  ";
  }
  std::cout << std::endl;
}

/*
  Simple program that uses the LookupTableGenerator class to
  print x y values to std::cout
*/
int main(int argc, char* argv[]){
  if(argc != 5){
    print_usage();
    exit(1);
  }

  std::string tableKey = argv[1];
  double tableMin  = std::stod(argv[2]);
  double tableMax  = std::stod(argv[3]);
  double tableStep = std::stod(argv[4]);

  func::FunctionContainer<double> func_container{FUNC_SET_F(MyFunction,double)};
  func::LookupTableGenerator<double> gen(func_container,tableMin,tableMax);
  gen.plot_implementation_at_step_size(tableKey,tableStep);
}

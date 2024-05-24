#include <iostream>
#include "func.hpp"
#include "chaste_log_function.hpp"

void print_usage(){
  std::cout << "Usage:\n"
    << "    ./experiment <tableKey> <tableMin> <tableMax> <tableStep> <plotRefinement>"
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
  Use LookupTableGenerator to print x y values to std::cout
*/
int main(int argc, char* argv[]){
  if(argc != 5 && argc != 6){
    print_usage();
    exit(1);
  }

  std::string tableKey = argv[1];
  double tableMin  = std::stod(argv[2]);
  double tableMax  = std::stod(argv[3]);
  double tableStep = std::stod(argv[4]);
  double plotRefinement; // how many points sampled per subinterval
  if (argc == 6)
    plotRefinement = std::stod(argv[5]);
  else
    plotRefinement = 100;

  func::FunctionContainer<double> func_container{FUNC_SET_F(MyFunction,double)};
  func::LookupTableParameters<double> par {tableMin, tableMax, 0.0};
  par.special_points = {{std::exp(7.7/13.0287), 0, 0.0}};
  func::LookupTableGenerator<double> gen(func_container, par);
  gen.plot_implementation_at_step_size(tableKey, tableStep, plotRefinement);
}

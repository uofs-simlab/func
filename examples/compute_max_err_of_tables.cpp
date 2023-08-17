#include <iostream>
#include "func.hpp"
#include "chaste_log_function.hpp"

#define TYPE double
#define MIN_ARG 0.001
#define MAX_ARG 30.0

/*
  Simple program that uses the LookupTableGenerator class to
  compute errors for varying table step sizes
*/
int main(){
  using namespace func;
  FunctionContainer<TYPE> func_container{FUNC_SET_F(MyFunction,TYPE)};

  /* Use every registered LookupTable! */
  std::cout << "Function: " << FUNCNAME << std::endl;
  LookupTableFactory<TYPE> factory;
  LookupTableGenerator<TYPE> gen(func_container, MIN_ARG, MAX_ARG);
  for(auto lutName : factory.get_registered_keys()){
    std::cout << "h vs E(h) for " << lutName << " with rtol=atol=1." << std::endl;
    for(TYPE h = 1.0; h>0.01; h /= 2.0){
      TYPE err = gen.error_at_step_size(lutName,h);
      std::cout << h << " " << err << std::endl;
    }
    std::cout << std::endl;
  }
}

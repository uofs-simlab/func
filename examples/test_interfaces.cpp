/**
 * \file Test if FunC can build and evaluate LUT containers (FailureProofTable, and CompositeLookupTable)
 */

#include "func.hpp"
#include "chaste_log_function.hpp"
#include <iostream>
#define TYPE double
#define MIN_ARG 0.1
#define MAX_ARG 30
#define STEP 0.1
#define NEVALS 10000
#define SEED 5
#define NEXPERIMENTS 10

int main()
{
  using namespace func;

  FunctionContainer<TYPE> func_container{FUNC_SET_F(MyFunction,TYPE)};

  /* Build a FailureProofTable */
  // (std::exp(7.7/13.0287)-MIN_ARG)/17 \approx 0.1
  FailureProofTable<UniformEqSpaceInterpTable<3,double>,double> F(func_container, {MIN_ARG,MAX_ARG,(std::exp(7.7/13.0287)-MIN_ARG)/17});
  std::cout << "F(0.01) = " << F(0.01) << std::endl;
  std::cout << "F(1)  = " << F(1)  << std::endl;
  std::cout << "F(std::exp(7.7/13.0287))  = " << F(std::exp(7.7/13.0287))  << std::endl;
  std::cout << "F(50) = " << F(50) <<  std::endl;

  std::vector<std::tuple<std::string,TYPE,TYPE,TYPE>> v{
    /*{tableKey, left, right, step},*/
    {"UniformEqSpaceInterpTable<3>",MIN_ARG,std::exp(7.7/13.0287),STEP},
    {"UniformEqSpaceInterpTable<3>",std::exp(7.7/13.0287),MAX_ARG,STEP}
  };

  CompositeLookupTable<TYPE> C(func_container, v);
  std::cout << "C(0.01) = " << C(0.01) << std::endl;
  std::cout << "C(1)  = " << C(1)  << std::endl;
  std::cout << "C(std::exp(7.7/13.0287))  = " << C(std::exp(7.7/13.0287))  << std::endl;
  std::cout << "C(29) = " << C(29) << std::endl;
  std::cout << "C(50) = " << C(50) << std::endl;

  LookupTableGenerator<TYPE> gen(func_container, MIN_ARG, MAX_ARG);
  std::cout << gen.error_of_table(F) << std::endl;
  std::cout << gen.error_of_table(C) << std::endl;

  // TODO 
  ///* copy the above objects into unique_ptrs so we can use a LookupTableComparator */
  //std::vector<std::unique_ptr<LookupTable<TYPE>>> impls;
  //impls.emplace_back(std::unique_ptr<LookupTable<TYPE>>(new FailureProofTable<UniformEqSpaceInterpTable<3,double>,double>(F)));
  //impls.emplace_back(std::unique_ptr<LookupTable<TYPE>>(new CompositeLookupTable<TYPE>(C)));

  ///* Time how long it takes to apply the LUTs to large random vectors a couple times and summarize the statistics */
  //LookupTableComparator<TYPE> implCompare(impls, MIN_ARG, MAX_ARG, NEVALS, SEED);
  //implCompare.run_timings(NEXPERIMENTS);
  //implCompare.compute_statistics();
  //implCompare.sort_timings(Sorter::WORST);
  //implCompare.print_summary(std::cout);
}

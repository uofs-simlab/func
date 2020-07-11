/*
  Direct evaluation of a function

  Usage example:
    DirectEvaluation de(&function,0,10);
    double f = de(0.87354);

  Notes:
  - basically just a wrapper around a function
  - requires a max and min value (for consistency with EvaluationImplementation)
  - Record where the function is being evaluated by specifying the -DFUNC_RECORD flag at compile time
*/
#pragma once

#include "EvaluationImplementation.hpp"
#include "FunctionContainer.hpp"
#include <cmath>

#ifdef FUNC_RECORD
#include "ArgumentRecord.hpp"
#endif

class DirectEvaluation final : public EvaluationImplementation
{  
private:
  #ifdef FUNC_RECORD
  std::unique_ptr<ArgumentRecord> mp_recorder;
  #endif
public:
  DirectEvaluation(FunctionContainer *func_container, double min = 0, double max = 1, unsigned int histSize = 10);
  double operator()(double x) override;
  void print_details(std::ostream& out) override;
  void print_details_json(std::ostream& out);
  ~DirectEvaluation();
};

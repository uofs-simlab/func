/*
  A wrapper for a standard func table class. If an argument is outside a
  table's range then the original function is used and the arg is recorded.

  Usage example:
    FailureProofTable<double> failsafe(unique_ptr<UniformLookupTable<double>>(
      new UniformCubicPrecomputedInterpolationTable<double>(&function,0,10,0.0001))
    );
    double val = failsafe(0.87354);
    // or
    FailureProofTable<double,double,UniformCubicPrecomputedInterpolationTable<double>> failsafe(
      &function,{0,10,0.0001}
    );
    double val = failsafe(0.87354);

  Notes:
  - ownership of the original table is moved to this class upon construction
  (not a problem since tables don't have move constructors for unique_ptr)
  - specify the LUT_TYPE in the third template arg is recommended because it
  avoids a virtual operator() call, and doing so lets you built this table like
  any other FunC lookup table
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
  - specify the NDEBUG flag to turn off argument recording for args outside
  the table's range.
*/
#pragma once
#include "EvaluationImplementation.hpp"
#include "FunctionContainer.hpp"
#include "UniformLookupTable.hpp"
#include "json.hpp"
#include <fstream> //ifstream
#include <limits> // std::numeric_limits<IN_TYPE>::max() / min()

#ifdef FUNC_RECORD
  #include "ArgumentRecord.hpp"
#endif

template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE, class LUT_TYPE = UniformLookupTable<IN_TYPE,OUT_TYPE>>
class FailureProofTable final : public EvaluationImplementation<IN_TYPE,OUT_TYPE> {
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  #ifdef FUNC_RECORD
    std::unique_ptr<ArgumentRecord<IN_TYPE>> mp_recorder;
  #endif

  std::unique_ptr<LUT_TYPE> mp_LUT;
public:
  /* Steal the given LUTs identity */
  FailureProofTable(std::unique_ptr<LUT_TYPE> LUT, unsigned int histSize = 10) :
    mp_LUT(std::move(LUT)),
    EvaluationImplementation<IN_TYPE,OUT_TYPE>(mp_LUT->function(), mp_LUT->name())
  {
    m_minArg = mp_LUT->min_arg();
    m_maxArg = mp_LUT->max_arg();
    m_order  = mp_LUT->order();
    m_dataSize = mp_LUT->size();
    #ifdef FUNC_RECORD
      mp_recorder.reset(new ArgumentRecord<IN_TYPE>(histSize,
            std::numeric_limits<IN_TYPE>::min(), std::numeric_limits::max()));
    #endif
    (void) histSize; // ignore histSize
  }

  /* Build our own LUT_TYPE */
  FailureProofTable(FunctionContainer<IN_TYPE,OUT_TYPE> *fc, UniformLookupTableParameters<IN_TYPE> par, unsigned int histSize = 10) :
    FailureProofTable(std::unique_ptr<LUT_TYPE>(new LUT_TYPE(fc,par)), histSize) {}

  /* Build our own LUT_TYPE from a file */
  FailureProofTable(FunctionContainer<IN_TYPE,OUT_TYPE> *fc, std::string filename) :
    FailureProofTable(std::unique_ptr<LUT_TYPE>(new LUT_TYPE(fc,filename)))
  {
    std::ifstream file_reader(filename);
    using nlohmann::json;
    json jsonStats;
    file_reader >> jsonStats;
    #ifdef FUNC_RECORD
      // reconstruct our arg record
      mp_recorder = std::unique_ptr<ArgumentRecord<IN_TYPE>>(new ArgumentRecord<IN_TYPE>(jsonStats));
    #endif
  }

  /* if x isn't contained within the tables bounds,
   * then resort to evaluating the original function */
  OUT_TYPE operator()(IN_TYPE x) override
  {
    // check if x is in the range of the table
    if(x < this->m_minArg || this->m_maxArg < x){
      #ifdef FUNC_RECORD
        mp_recorder->record_arg(x);
      #endif
      return this->mp_func(x);
    }
    return (*(this->mp_LUT))(x);
  }

  void print_details(std::ostream &out) override 
  {
    out << "FailureProof" << this->m_name << " " << this->m_minArg << " " << this->m_maxArg << " "
        << mp_LUT->step_size() << " " << mp_LUT->num_intervals() << " ";
  #ifdef FUNC_RECORD
    out << std::endl;
    mp_recorder->print_details(out);
  #endif
  }

  void print_details_json(std::ostream &out) override
  {
    using nlohmann::json;
    json jsonStats;

    jsonStats["_comment"] = "FunC FailureProofTable data";
    jsonStats["name"] = m_name;
    jsonStats["minArg"] = m_minArg;
    jsonStats["maxArg"] = m_maxArg;
    #ifdef FUNC_RECORD
      // have our ArgumentRecord add it's own data
      mp_recorder->print_details_json(jsonStats);
    #endif
    
    out << jsonStats.dump(2) << std::endl;
  }  
};

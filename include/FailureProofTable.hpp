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
  - specify the FUNC_DEBUG flag to turn on argument recording for args outside
  the table's range.
  - optional ArgumentRecord args available if you want nicer looking output
*/
#pragma once
#include "EvaluationImplementation.hpp"
#include "FunctionContainer.hpp"
#include "LookupTable.hpp"
#include "LookupTableGenerator.hpp" // generate_by_filename
#include "json.hpp"
#include <memory> // unique_ptr
#include <fstream> //ifstream
#include <limits> // std::numeric_limits<IN_TYPE>::max() / lowest()

#ifdef FUNC_DEBUG
  #include "ArgumentRecord.hpp"
#endif

namespace func {

template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE, class LUT_TYPE = LookupTable<IN_TYPE,OUT_TYPE>>
class FailureProofTable final : public EvaluationImplementation<IN_TYPE,OUT_TYPE> {
  std::unique_ptr<LUT_TYPE> mp_LUT;
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  #ifdef FUNC_DEBUG
    std::unique_ptr<ArgumentRecord<IN_TYPE>> mp_recorder;
  #endif
public:
  /* Steal the given LUTs identity */
  FailureProofTable(std::unique_ptr<LUT_TYPE> LUT,
      IN_TYPE histMin = 1,
      IN_TYPE histMax = 0,
      unsigned int histSize = 10
      ) :
    mp_LUT(std::move(LUT))
  {
    // m_func and m_name can't be set in the super class b/c the
    // base class constructor would be evaluated before mp_LUT is set
    m_func   = mp_LUT->function();
    m_name   = mp_LUT->name();
    m_minArg = mp_LUT->min_arg();
    m_maxArg = mp_LUT->max_arg();
    m_order  = mp_LUT->order();
    m_dataSize = mp_LUT->size();
    #ifdef FUNC_DEBUG
      // check if we're using default/bad histogram arguments
      if(histMin >= histMax){
        histMin = m_minArg;
        histMax = m_maxArg;
      }
      mp_recorder.reset(new ArgumentRecord<IN_TYPE>(
            histMin, histMax, histSize
            ));
    #endif
    // ignore hist parameters if they're unused
    (void) histMin;
    (void) histMax;
    (void) histSize;
  }

  /* Build our own LUT_TYPE. Only works if the template is specific enough */
  FailureProofTable(FunctionContainer<IN_TYPE,OUT_TYPE> *fc,
      LookupTableParameters<IN_TYPE> par,
      IN_TYPE histMin = 1, IN_TYPE histMax = 0, unsigned int histSize = 10) :
    FailureProofTable(std::unique_ptr<LUT_TYPE>(new LUT_TYPE(fc,par)), histMin, histMax, histSize) {}

  /* Build our own LUT_TYPE from a file */
  FailureProofTable(FunctionContainer<IN_TYPE,OUT_TYPE> *fc, std::string filename) :
    FailureProofTable(LookupTableGenerator<IN_TYPE,OUT_TYPE>(fc,1,0).generate_by_file(filename))
  {
    nlohmann::json jsonStats;
    std::ifstream(filename) >> jsonStats;

    #ifdef FUNC_DEBUG
      // reconstruct our arg record if it was ever saved to jsonStats
      try{
        mp_recorder =
          std::unique_ptr<ArgumentRecord<IN_TYPE>>(new ArgumentRecord<IN_TYPE>(jsonStats));
      } catch (const nlohmann::detail::type_error&) {}
    #endif
  }

  /* if x isn't contained within the tables bounds,
   * then resort to evaluating the original function */
  OUT_TYPE operator()(IN_TYPE x) override
  {
    // check if x is in the range of the table
    if(x < m_minArg || m_maxArg < x){
      #ifdef FUNC_DEBUG
        mp_recorder->record_arg(x);
      #endif
      return m_func(x);
    }
    return (*mp_LUT)(x);
  }

  void print_details(std::ostream &out) override 
  {
    out << "FailureProof" << m_name << " " << m_minArg << " " << m_maxArg << " "
        << mp_LUT->step_size() << " " << mp_LUT->num_intervals() << " ";
  #ifdef FUNC_DEBUG
    out << std::endl;
    mp_recorder->print_details(out);
  #endif
  }

  // TODO test
  void print_details_json(std::ostream &out) override
  {
    nlohmann::json jsonStats;

    jsonStats["_comment"] = "FunC FailureProofTable data";
    jsonStats["name"] = m_name;
    jsonStats["minArg"] = m_minArg;
    jsonStats["maxArg"] = m_maxArg;
    #ifdef FUNC_DEBUG
      // have our ArgumentRecord add it's own data
      mp_recorder->print_details_json(jsonStats);
    #endif
    mp_LUT->print_details_json(out);
    
    out << jsonStats.dump(2) << std::endl;
  }  
};
} // namespace func
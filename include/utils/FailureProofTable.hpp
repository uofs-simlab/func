/* TODO don't wrap a smart pointer. Just store a LUT directly!
 *
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

  TODO this class should support to_json but not from_json! Add another constructor
  to build this class from a FunctionContainer and a filename
*/
#pragma once
#include "FunctionContainer.hpp"
#include "LookupTableGenerator.hpp"
#include "LookupTable.hpp"
#include "json.hpp"
#include <memory> // unique_ptr
#include <fstream> //ifstream
#include <limits> // std::numeric_limits<TIN>::max() / lowest()

#ifdef FUNC_DEBUG
  #include "ArgumentRecord.hpp"
#endif

namespace func {

template <class LUT_TYPE, typename TIN, typename TOUT = TIN>
class FailureProofTable final : public LookupTable<TIN,TOUT> {
  std::function<TOUT(TIN)> m_fun;
  LUT_TYPE m_LUT;
  #ifdef FUNC_DEBUG
  std::unique_ptr<mutable ArgumentRecord<TIN>> mp_recorder;
  #endif
public:
  /* Deep copy the given LUT */
  //FailureProofTable(const FunctionContainer<TIN,TOUT>& fc, const LUT_TYPE& LUT,
  //    TIN histMin = 1, TIN histMax = 0, unsigned int histSize = 10) : m_LUT(LUT)
  //{
  //  m_fun = fc.standard_fun;
  //  #ifdef FUNC_DEBUG
  //    // check if we're using the default (aka bad) histogram arguments
  //    if(histMin >= histMax){
  //      histMin = m_LUT.min_arg();
  //      histMax = m_LUT.max_arg();
  //    }
  //    mp_recorder.reset(new ArgumentRecord<TIN>(histMin, histMax, histSize));
  //  #endif
  //  // ignore histogram parameters if they're unused
  //  (void) histMin; (void) histMax; (void) histSize;
  //}

  /* Build our own LUT_TYPE. Only works if the template is specific enough */
  FailureProofTable(const FunctionContainer<TIN,TOUT>& fc, const LookupTableParameters<TIN>& par,
      TIN histMin = 1, TIN histMax = 0, unsigned int histSize = 10) : m_LUT(fc, par)
  {
  m_fun = fc.standard_fun;
  #ifdef FUNC_DEBUG
  // check if we're using the default (aka bad) histogram arguments
  if(histMin >= histMax){
    histMin = m_LUT.min_arg();
    histMax = m_LUT.max_arg();
  }
  mp_recorder.reset(new ArgumentRecord<TIN>(histMin, histMax, histSize));
  #endif
  // ignore histogram parameters if they're unused
  (void) histMin; (void) histMax; (void) histSize;
  }

  /* if x isn't in the LUT's bounds, then return f(x) */
  TOUT operator()(TIN x) const final
  {
    // check if x is in the range of the table
    //if((x - m_minArg)*(m_maxArg - x) > 0){ // is this a possible micro-optimization?
    if((m_LUT.min_arg() < x) && (x < m_LUT.max_arg())){
      return m_LUT(x);
    }else{
      #ifdef FUNC_DEBUG
        mp_recorder->record_arg(x);
      #endif
      return m_fun(x);
    }
  }

  std::string name() const final { return std::string("FailureProof") + m_LUT.name(); }
  TIN min_arg() const final { return m_LUT.min_arg(); }
  TIN max_arg() const final { return m_LUT.max_arg();}
  unsigned int order() const final { return m_LUT.order();}
  unsigned int size() const final { return m_LUT.size();}
  unsigned int num_subintervals() const final { return m_LUT.num_subintervals();}
  TIN step_size() const final { return m_LUT.step_size(); }
  std::pair<TIN,TIN> bounds_of_subinterval(unsigned int intervalNumber) const final { return m_LUT.bounds_of_subinterval(intervalNumber);}
  void print_json(std::ostream& out) const final { (void) out; /* TODO call to_json() */ };
};
// TODO make to_json()
} // namespace func

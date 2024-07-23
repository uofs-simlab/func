/** \brief A wrapper for any implementation of LookupTable L. The
   operator()(x) ensures x is within the bounds of L before returning L(x).
   Returns f(x) for out of bounds arguments. If FUNC_DEBUG is defined then
   out of bounds arguments are recorded in a histogram.


   \tparam LUT_TYPE is a specific implementation of LookupTable (eg. ChebyInterpTable<3,double>)

  TODO broken
  Usage example:
    FailureProofTable<UniformChebyInterpTable<3,double>> failsafe(
      &function,{0,10,0.0001}
    );
    double val = failsafe(0.87354);
    double val = failsafe(100);

  Notes:
  - data is static after constructor call
  - optional ArgumentRecord args available to improve binning (max & min is accuratly tracked)

  TODO this class will support to_json but not from_json. Add another constructor
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
#include <boost/math/special_functions/sign.hpp>

#ifdef FUNC_DEBUG
  #include "ArgumentRecord.hpp"
#endif

namespace func {

template <class LUT_TYPE>
class FailureProofTable final : public LookupTable<typename LUT_TYPE::input_type, typename LUT_TYPE::output_type> {
  using TIN = typename LUT_TYPE::input_type;
  using TOUT = typename LUT_TYPE::output_type;

  LUT_TYPE m_LUT;
  std::function<TIN(TOUT)> m_fun;
  
  #ifdef FUNC_DEBUG
    mutable std::unique_ptr<ArgumentRecord<TIN>> mp_recorder;
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
      TIN histMin = 1, TIN histMax = 0, unsigned int histSize = 10, std::ostream* streamer = nullptr) :
    m_LUT(fc, par)
  {
    m_fun = fc.standard_fun;
    #ifdef FUNC_DEBUG
    // check if we're using the default (aka bad) histogram arguments
    if(histMin >= histMax){
      histMin = m_LUT.min_arg();
      histMax = m_LUT.max_arg();
    }
    mp_recorder.reset(new ArgumentRecord<TIN>(histMin, histMax, histSize, streamer));
    #endif
    // ignore histogram parameters if they're unused
    (void) histMin; (void) histMax; (void) histSize; (void) streamer;
  }

  /* if x isn't in the LUT's bounds, then return f(x) */
  TOUT operator()(TIN x) const final
  {
    // check if x is in the range of the LUT
    // try boost sign????
    //if(boost::math::sign(x - m_LUT.min_arg()) * boost::math::sign(m_LUT.max_arg() - x) > 0){
    // if((x - m_LUT.min_arg())*(m_LUT.max_arg() - x) > 0){
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
  std::size_t size() const final { return m_LUT.size();}
  unsigned int num_subintervals() const final { return m_LUT.num_subintervals();}
  TIN step_size() const final { return m_LUT.step_size(); }
  std::pair<TIN,TIN> bounds_of_subinterval(unsigned int intervalNumber) const final { return m_LUT.bounds_of_subinterval(intervalNumber);}
  
  void print_json(std::ostream& out) const final
  { 
    (void) out; 
    #ifdef FUNC_DEBUG
      out << mp_recorder->print_json(out);
    #endif
  };
};

} // namespace func

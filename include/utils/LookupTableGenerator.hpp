/*
  Generate a FunC LookupTable when given that table's name and one of the following:
  - stepsize
  - tolerance
  - memory size limit
  - filename
  If gen_by_XXX is given a filename then it will generate a table and save that output to
  to filename. Future runs will use the table in filename instead of generating that table
  from scratch again. filename will be with respect to the cwd unless users provide an
  absolute path.

  If Boost is not available then users can only build tables by file.

  Provided as a header only class because it is templated on the error precision TERR.
  We MUST be able to cast TERR to TIN and vice versa.
  Ideally TERR satisfies: sqrt(epsilon_TERR) <= epsilon_TOUT.

  This class is also equipped to
  - compute table error estimates at a given stepsize
  - plot a table implementation against the exact function

  TODO:
  - Newton's iterate is currently unused because sometimes it'll try
    building a LUT with a stepsize that will mortal computers.
  - It might behave better with the new error estimate? or if we start approximating function's roots
*/
#pragma once
#include "LookupTable.hpp"
#include "LookupTableFactory.hpp"
#include "config.hpp" // FUNC_USE_QUADMATH, FUNC_USE_BOOST

#include <string>
#include <memory>
#include <limits>
#include <stdexcept>
#include <cmath> // std::abs, std::min
//#include <algorithm> // std::min

#ifdef FUNC_USE_BOOST
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#endif

namespace func {


#if defined(FUNC_USE_BOOST)
template <typename TIN, typename TOUT = TIN, typename TERR = boost::multiprecision::cpp_bin_float_quad>
#else
template <typename TIN, typename TOUT = TIN, typename TERR = long double> // TERR is unused if Boost is not available
#endif
class LookupTableGenerator
{
private:
  FunctionContainer<TIN,TOUT> m_fc;
  LookupTableParameters<TIN> m_par;
  TIN m_min, m_max; // min/max member variables are convenient but they should be removed in favour of m_par...

  LookupTableFactory<TIN,TOUT> factory;

  /* Nested functor for error evaluation */
  struct LookupTableErrorFunctor;

  /* Nested functor for optimal grid spacing determination */
  struct OptimalStepSizeFunctor;


  /* check if filename exists or if the user can access it */
  bool file_exists(std::string filename){
    return static_cast<bool>(std::ifstream(filename));
  }

  inline void save_lut(LookupTable<TIN,TOUT>* lut, std::string filename){
    if(filename == "") return;

    std::ofstream out_file(filename);
    lut->print_json(out_file);
    return;
  }

public:
  LookupTableGenerator(const FunctionContainer<TIN,TOUT>& fc, const LookupTableParameters<TIN>& par) :
    m_fc(fc), m_par(par), m_min(par.minArg), m_max(par.maxArg) {}

  LookupTableGenerator(const FunctionContainer<TIN,TOUT>& fc, TIN minArg, TIN maxArg) :
    LookupTableGenerator(fc, {minArg, maxArg, static_cast<TIN>(0.0)}) {}

  ~LookupTableGenerator(){}

  /* A wrapper for the LookupTableFactory<std::string> which builds tables from a file
   * tableKey arg only exists as a sanity check (it's pointless otherwise) */
  std::unique_ptr<LookupTable<TIN,TOUT>> generate_by_file(std::string filename, std::string tableKey = "")
  {
    if(filename.find(".json") == std::string::npos) // TODO are there any other json filename extensions?
      throw std::invalid_argument("Error in func::LookupTableGenerator.generate_by_file: filename is not a valid json file.");

    nlohmann::json jsonStats;
    std::ifstream(filename) >> jsonStats;
    // get the tableKey from filename
    if(tableKey == ""){
      tableKey = jsonStats["name"].get<std::string>();
    }
    // MetaTable will check that tableKey actually matches the name in filename
    return factory.create(tableKey, m_fc, LookupTableParameters<TIN>{0,0,0}, jsonStats);
  }

  /* A wrapper for the LookupTableFactory */
  std::unique_ptr<LookupTable<TIN,TOUT>> generate_by_step(std::string tableKey, TIN stepSize, std::string filename = "")
  {
    if(filename != "" && file_exists(filename))
      return generate_by_file(filename, tableKey);

    // LookupTable will make sure the stepsize is positive
    LookupTableParameters<TIN> par = m_par;
    par.stepSize = stepSize;

    auto lut = factory.create(tableKey, m_fc, par);
    save_lut(lut.get(), filename);
    return lut;
  }

  /* Generate a table that has the largest possible stepsize such that the error is less than desiredErr */
  std::unique_ptr<LookupTable<TIN,TOUT>> generate_by_tol(std::string tableKey, TIN a_tol, TIN r_tol, std::string filename = "");
  std::unique_ptr<LookupTable<TIN,TOUT>> generate_by_tol(std::string tableKey, TIN desiredErr, std::string filename = ""){
    return generate_by_tol(tableKey,desiredErr,desiredErr,filename);
  }

  /* Generate a table takes up desiredSize bytes */
  std::unique_ptr<LookupTable<TIN,TOUT>> generate_by_impl_size(std::string tableKey, unsigned long desiredSize, std::string filename = "");

  /* Return the approx error in tableKey at stepSize
   * - relTol is a parameter which determines how much effect small f(x) values have on the error calculation */
  long double error_at_step_size(std::string tableKey, TIN stepSize, TIN relTol = static_cast<TIN>(1.0));

  /* compare tableKey to the original function at stepSize */
  void plot_implementation_at_step_size(std::string tableKey, TIN stepSize);

  TIN min_arg(){ return m_par.minArg; }
  TIN max_arg(){ return m_par.maxArg; }
};

/*----------------------------------------------------------------------------*/
/* Non-trivial member definitions */
/*----------------------------------------------------------------------------*/

/* lambda used for computing error in a given lookup table */
template <typename TIN, typename TOUT, typename TERR>
struct LookupTableGenerator<TIN,TOUT,TERR>::LookupTableErrorFunctor
{
  LookupTableErrorFunctor(LookupTable<TIN,TOUT>* impl, std::function<TOUT(TIN)> fun, TERR relTol) :
    m_fun(fun), m_impl(impl), m_relTol(relTol) {}

  /* Compute -|f(x) - L(x)| / (1 + r_tol/a_tol*|f(x)|). Notes:
   * - We're maximizing this function with brent_find_minima so operator() must always return a negative value
   * - Only parameterized on the tolerance relative to a_tol (needed for error_at_step_size()) */
  TERR operator()(TERR const& x)
  {
    TERR f_value = static_cast<TERR>(m_fun(static_cast<TIN>(x)));
    TERR lut_value = static_cast<TERR>((*m_impl)(static_cast<TIN>(x)));
    return -fabs(f_value - lut_value) / (static_cast<TERR>(1.0) + m_relTol*fabs(f_value));
  }

private:

  std::function<TOUT(TIN)> m_fun;
  LookupTable<TIN,TOUT> *m_impl;
  TERR m_relTol;
};

/* lambda used for finding error over table bounds */
template <typename TIN, typename TOUT, typename TERR>
struct LookupTableGenerator<TIN,TOUT,TERR>::OptimalStepSizeFunctor
{
  // small relTol => don't worry about small f(x) in error metric. large relTol => small |f(x)| must be approximated very well!
  // Set relTol = r_tol/a_tol.
  OptimalStepSizeFunctor(LookupTableGenerator<TIN,TOUT,TERR> &parent, std::string tableKey, TERR relTol, TERR desiredErr) :
    m_parent(parent), m_tableKey(tableKey), m_relTol(relTol), m_desiredErr(desiredErr) {}

  // returns something close to 0 when err is close to the desiredErr
  TIN operator()(TIN const& stepSize)
  {
    using namespace boost::math::tools;

    LookupTableParameters<TIN> par = m_parent.m_par; // deep copy
		par.stepSize = stepSize;
    auto impl = m_parent.factory.create(m_tableKey, m_parent.m_fc, par);

    /* get number of binary bits in mantissa */
    int bits = std::numeric_limits<TERR>::digits/2; // effective maximum
    boost::uintmax_t max_it = 20;

    TERR max_err = 0;

    /* Want a small bracket for brent's method so for each interval in the table,
     * compute the maximum error.
     * - Be careful about the top most interval b/c tableMaxArg can be greater than max
     *   (and we don't care about error outside of table bounds)
     * - TODO Parallelizing this for loop is worthwhile because software implementations of
     *   high precision floats are quite slow. Is this the best pragma possible?
     * - TODO can be too slow for high order tables with very few subintervals
     *   */
    #pragma omp parallel for
    for(unsigned ii=0; ii<impl->num_subintervals(); ii++){
      std::pair<TIN,TIN> intEndPoints = impl->bounds_of_subinterval(ii);
      TERR x = static_cast<TERR>(boost::math::float_next(intEndPoints.first));
      TERR xtop = static_cast<TERR>(boost::math::float_prior(intEndPoints.second));

      /* tableMaxArg >= max. We don't care about any error outside [min,max]
       * so ignore any points between (max,tableMax] */
      if(static_cast<TIN>(xtop) > m_parent.max_arg())
        xtop = static_cast<TERR>(m_parent.max_arg());

      std::pair<TERR, TERR> r = brent_find_minima(LookupTableErrorFunctor(impl.get(),m_parent.m_fc.standard_fun,m_relTol),x,xtop,bits,max_it);

      #pragma omp critical
      {
        max_err = std::min(max_err, r.second);
      }
      //std::cerr << -err << " error at x=" << r.first << "\n";
    }

    /* want return to be 0 if the same, +/- on either side */
    max_err = -max_err;
    //std::cerr << m_tableKey << " with stepSize=" << stepSize << " has approx error=" << double(max_err) << "\n";
    return static_cast<TIN>(max_err-m_desiredErr);
  }

private:

  LookupTableGenerator<TIN,TOUT,TERR>& m_parent;
  std::string m_tableKey;
  TERR m_relTol;
  TERR m_desiredErr;

};

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*
  LookupTableGenerator functions
*/

template <typename TIN, typename TOUT, typename TERR>
std::unique_ptr<LookupTable<TIN,TOUT>> LookupTableGenerator<TIN,TOUT,TERR>::generate_by_impl_size(
    std::string tableKey, unsigned long desiredSize, std::string filename)
{
  if(filename != "" && file_exists(filename))
    return generate_by_file(filename, tableKey);

  /* Use 2 query points to get relationship */
  const unsigned long N1 = 2;
  const TIN step1 = (m_max-m_min)/N1;
  const unsigned long N2 = 10;
  const TIN step2 = (m_max-m_min)/N2;

  LookupTableParameters<TIN> par1 = m_par;
  par1.stepSize = step1;

  LookupTableParameters<TIN> par2 = m_par;
  par2.stepSize = step2;

  auto impl1 = factory.create(tableKey, m_fc, par1);
  auto impl2 = factory.create(tableKey, m_fc, par2);

  unsigned long size1 = impl1->size();
  unsigned long size2 = impl2->size();

  // TODO logically, this can't ever be a problem... but I guess it's nice to check?
  if(size2 == size1)
    throw std::logic_error("Error in func::LookupTableGenerator.generate_by_impl_size: Query tables have same size");

  /* approximate step size for for desired impl size
   * (assuming linear relationship of num_subintervals to size */
  const unsigned long N3 = (N2-N1)*(desiredSize-size1)/(size2-size1) + N1;
  par1.stepSize = (m_max-m_min)/static_cast<TIN>(N3);

  if(par1.stepSize <= 0)
    throw std::invalid_argument("Error in func::LookupTableGenerator.generate_by_impl_size: Requested memory size is too small");

  auto lut = factory.create(tableKey, m_fc, par1);
  save_lut(lut.get(), filename);
  return lut;
}

template <typename TIN, typename TOUT, typename TERR>
std::unique_ptr<LookupTable<TIN,TOUT>> LookupTableGenerator<TIN,TOUT,TERR>::generate_by_tol(std::string tableKey, TIN a_tol, TIN r_tol, std::string filename)
{
#ifndef FUNC_USE_BOOST
    static_assert(sizeof(TIN)!=sizeof(TIN), "Cannot generate any LUT by tol without Boost");
#else
  if(filename != "" && file_exists(filename))
    return generate_by_file(filename, tableKey);

  LookupTableParameters<TIN> par = m_par;
  par.stepSize = m_max-m_min; // max reasonable stepsize

  /* generate a first approximation for the implementation */
  auto impl = factory.create(tableKey, m_fc, par);
  /* And initialize the functor used for refinement:
   * Note relTol is considered the "relative tolerance relative to the absolute tolerance" */
  TERR relTol = r_tol/a_tol;
  OptimalStepSizeFunctor f(*this,tableKey,relTol,0.0);

  /* quit now if this table is amazing already. Useful for high order tables on small intervals
   * where bracket_and_solve tries to use a stepsize larger than the table range */
  auto fmax_step = f(m_max-m_min);
  if(fmax_step <= a_tol){
    save_lut(impl.get(), filename);
    //std::cerr << "estimated max error of " << fmax_step << " with stepsize of " << m_max-m_min << std::endl;
    return impl;
  }

  TIN stepSize = (m_max-m_min)/1000.0; // initial guess stepsize. Probably too small for most high order tables
  const int order  = impl->order();

  const TIN logTol = log(a_tol);

  /*
    APPLY NEWTON'S METHOD IN log-log SPACE
    ('known' slope = order)

     Approximate the solution that satisfies tolerance based on the
     order of the implementation
     - assumes that the initial guess is in the asymptotic regime
       of the table's convergence
     - even if not in asymptotic regime, this thing has some pretty
       robust convergence!
     - Currently unused. TODO test with new error metric!
  */
  const int  N_NEWTON_MAX_IT  = 0;     // max log-Newton-iterations
  const double NEWTON_IT_RTOL = 1e-5;
  const double NEWTON_IT_ATOL = 1e-10;

  int NEWTON_SUCCESS_FLAG = 0;
  std::vector<std::pair<double,double>> iterates;
  for (int iNewton = 0; iNewton < N_NEWTON_MAX_IT; iNewton++) {
    double err = f(stepSize);
    if (fabs(err-a_tol) <  fabs(a_tol)*NEWTON_IT_RTOL+NEWTON_IT_ATOL) {
      std::cout << "Newton iter: " << iNewton << "\n";
      NEWTON_SUCCESS_FLAG = 1;
      break;
    }
    double logStepSize = log(stepSize);
    double logError    = log(err);
    logStepSize += (logTol-logError)/order;
    stepSize = exp(logStepSize);
    iterates.push_back(std::make_pair(stepSize,err));
  }

  /*
    Output Newton iterates:
  */
  // double err = f(stepSize);
  // std::cout << "Iterates (normal space): step, error" << std::endl;
  // for (auto it : iterates) {
  //   std::cout << "    " << it.first << ", " << it.second << std::endl;
  // }

  /*
    APPLY A BRACKETING ALGORITHM IN log-log SPACE

    If a suitable bracket is found, this guarantees a solution below a_tol.
    Otherwise, bracket_and_solve throws an error.
  */
  const boost::uintmax_t BRACKET_MAX_IT = 50;  // Limit to maximum iterations.
  boost::uintmax_t it = BRACKET_MAX_IT;        // Initally our chosen max iterations, but updated with actual.

  /*
    Throw when the log-Newton method did not converge AND there are no
    bracketing iterations performed.
   */
  if(!NEWTON_SUCCESS_FLAG && !BRACKET_MAX_IT)
    throw std::logic_error("Error in func::LookupTableGenerator: No bracketing iterations specified. log-Newton method did not converge in " + std::to_string(N_NEWTON_MAX_IT) + " steps.");
  
  /*
    Use the guess step size as an initialization to bracket_and_solve
  */
  using namespace boost::math::tools;

  int digits = std::numeric_limits<TIN>::digits;  // Maximum possible binary digits accuracy for type T.
  // Some fraction of digits is used to control how accurate to try to make the result.
  int get_digits = digits-30;  // We have to have a non-zero interval at each step, so
                               // doesn't have to be so accurate maximum accuracy is digits - 1.
                               // But we also have to allow for inaccuracy in f(x), otherwise the last
                               // few iterations just thrash around.
  eps_tolerance<TIN> tol(get_digits); // Set the tolerance.
  OptimalStepSizeFunctor g(*this,tableKey,relTol,a_tol); // functor for solving log(E(h)/TOL) = 0. g = f-a_tol

  /*
    Run bisection
    - answer is taken to be the lower point of the bracketing interval
    - this guarantees a tolerance lower than desired
    - we have to watch out for stepsizes larger than the table range
    so we can't just use bracket_and_solve_root
  */
  //TIN factor = 2; // Mult/divide factor for bracket expansion when searching.
  //bool is_rising = 1; // The error curve should always be rising
  //std::pair<TIN, TIN> r = bracket_and_solve_root(g, stepSize, factor, is_rising, tol, it);

  /* f(m_max-m_min) - a_tol > 0 (otherwise we would've quit already) */
  std::pair<TIN, TIN> r = toms748_solve(g, static_cast<TIN>(0.0), m_max-m_min, -static_cast<TIN>(a_tol), static_cast<TIN>(fmax_step-a_tol), tol, it);

  /* Save and return the implementation with the desired stepSize */
  par.stepSize = r.first;
  auto lut = factory.create(tableKey,m_fc,par);
  save_lut(lut.get(), filename);
  return lut;
#endif
}

template <typename TIN, typename TOUT, typename TERR>
long double LookupTableGenerator<TIN,TOUT,TERR>::error_at_step_size(
    std::string tableKey, TIN stepSize, TIN relTol)
{
#ifndef FUNC_USE_BOOST
    static_assert(sizeof(TIN)!=sizeof(TIN), "Cannot compute error at step without Boost");
    return 0; // just so the compiler doesn't complain
#else
  /* Use brent_find_minima to compute the max error over [min,max] */
  OptimalStepSizeFunctor f(*this,tableKey,static_cast<TERR>(relTol),0.0);
  TERR err = f(stepSize);
  return static_cast<long double>(err);
#endif
}

template <typename TIN, typename TOUT, typename TERR>
void LookupTableGenerator<TIN,TOUT,TERR>::plot_implementation_at_step_size(std::string tableKey, TIN stepSize)
{
  /*
    Can be implemented in terms of the Functor used in solving for a specific
    tolerance, so that is reused.
  */
  LookupTableParameters<TIN> par = m_par;
  par.stepSize = stepSize;
  auto impl = factory.create(tableKey, m_fc, par);

  std::cout << "# x func impl" << std::endl;
  for (TIN x = impl->min_arg();
           x < impl->max_arg();
           x+= impl->step_size()/static_cast<double>(100)){
    std::cout << x << " " <<
      (m_fc.standard_fun)(x) << " " <<
      (*impl)(x) << std::endl;
  }
}

} // namespace func

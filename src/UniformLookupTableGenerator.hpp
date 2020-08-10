/*
  Generate a FunC UniformLookupTable by name
  and with a given stepsize, tolerance, or memory size limit

  NOTES:
    Also equipped to
    - compute table error estimates at a given stepsize
    - plot a table implementation against the exact function
*/
#pragma once
#include "UniformLookupTable.hpp"
#include "TransferFunction.hpp"
#include <iostream>
#include <memory>
#include <limits>

#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/special_functions/next.hpp>

// If quadmath is used, work in the boost::multiprecision namespace
#ifdef USE_QUADMATH
#include <boost/multiprecision/float128.hpp>
using namespace boost::multiprecision;
using errprecision = float128;
#else
using errprecision = double;
#endif

// If Armadillo is used, add functionality for generating nonuniform lookup tables
#ifdef USE_ARMADILLO
  #include "TransferFunctionSinh.hpp"
#endif


template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE>
class UniformLookupTableGenerator
{
private:
  FunctionContainer<IN_TYPE,OUT_TYPE> *mp_func_container;

  std::shared_ptr<TransferFunction<IN_TYPE>> mp_transfer_function;

  IN_TYPE m_min;
  IN_TYPE m_max;

  /* Nested functor for error evaluation */
  struct LookupTableErrorFunctor;

  /* Nested functor for optimal grid spacing determination */
  struct OptimalStepSizeFunctor;

public:
  /* set member variables */
  UniformLookupTableGenerator(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container,
      IN_TYPE minArg, IN_TYPE maxArg,
      std::shared_ptr<TransferFunction<IN_TYPE>> transferFunction = nullptr) :
    mp_func_container(func_container), m_min(minArg), m_max(maxArg), mp_transfer_function(transferFunction)
  {
  #ifdef USE_ARMADILLO
    if(mp_transfer_function == nullptr)
      mp_transfer_function.reset(new TransferFunctionSinh<IN_TYPE>(mp_func_container, m_min, m_max));
  #endif
  }

  ~UniformLookupTableGenerator(){}

  /* essentially a wrapper function for the UniformLookupTableFactory */
  std::unique_ptr<UniformLookupTable<IN_TYPE,OUT_TYPE>> generate_by_step(std::string tableKey, IN_TYPE stepSize)
  {
    UniformLookupTableParameters<IN_TYPE> par;
    par.minArg = m_min; 
    par.maxArg = m_max; 
    par.stepSize = stepSize;
    par.transferFunction = mp_transfer_function;
    return UniformLookupTableFactory<IN_TYPE,OUT_TYPE>::Create(tableKey, mp_func_container, par);
  }

  /* Generate a table that is accurate to desiredTolerance */
  std::unique_ptr<UniformLookupTable<IN_TYPE,OUT_TYPE>> generate_by_tol(std::string tableKey, double desiredTolerance);

  /* Generate a table takes up desiredSize bytes */
  std::unique_ptr<UniformLookupTable<IN_TYPE,OUT_TYPE>> generate_by_impl_size(std::string tableKey, unsigned long desiredSize);

  /* Return the approx error in tableKey at stepSize */
  double error_at_step_size(std::string tableKey, IN_TYPE stepSize);

  /* compare tableKey to the original function at stepSize */
  void plot_implementation_at_step_size(std::string tableKey, IN_TYPE stepSize);
};

/*----------------------------------------------------------------------------*/
/* Non-trivial member definitions */
/*----------------------------------------------------------------------------*/
/*
   Nested Functor used for computing error in a given lookup table
*/
template <typename IN_TYPE, typename OUT_TYPE>
struct UniformLookupTableGenerator<IN_TYPE,OUT_TYPE>::LookupTableErrorFunctor
{
  LookupTableErrorFunctor(UniformLookupTable<IN_TYPE,OUT_TYPE>* impl) : m_impl(impl) {}
  /* operator() always returns a negative value */
  errprecision operator()(errprecision const& x)
  {
    errprecision f_value = static_cast<errprecision>((m_impl->function())(double(x)));
    errprecision lut_value = static_cast<errprecision>((*m_impl)(double(x)));
    return -static_cast<errprecision>(2.0) * fabs( (f_value - lut_value) ) /
      (fabs(f_value)+fabs(lut_value));
  }

private:

  UniformLookupTable<IN_TYPE,OUT_TYPE> *m_impl;
};

/* Nested Functor used for finding optimal stepsize that satisfies TOL */
template <typename IN_TYPE, typename OUT_TYPE>
struct UniformLookupTableGenerator<IN_TYPE,OUT_TYPE>::OptimalStepSizeFunctor
{
  OptimalStepSizeFunctor(UniformLookupTableGenerator<IN_TYPE,OUT_TYPE> &parent, std::string tableKey, double tol) :
    m_parent(parent), m_tableKey(tableKey), m_tol(tol) {}

  double operator()(IN_TYPE const& stepSize)
  {
    using namespace boost::math::tools;

    UniformLookupTableParameters<IN_TYPE> par;
    par.minArg = m_parent.m_min;
		par.maxArg = m_parent.m_max;
		par.stepSize = stepSize;
    par.transferFunction = m_parent.mp_transfer_function;
    auto impl = UniformLookupTableFactory<IN_TYPE,OUT_TYPE>::Create(m_tableKey, m_parent.mp_func_container, par);

    boost::uintmax_t max_it = 20;

    errprecision max_err = 0;
    errprecision xstar, err;

    /* get number of binary bits in mantissa */
    int bits = std::numeric_limits<errprecision>::digits;

    double eps = std::numeric_limits<double>::epsilon();

    /*
      for each interval in the uniform table, compute the maximum error
      - be careful about the top most interval, it may reach beyond the
      table range due to rounding errors
    */
    unsigned index = 0;
    for(unsigned ii=0; ii<impl->num_intervals()-1; ii++){

      std::pair<IN_TYPE,IN_TYPE> intEndPoints = impl->arg_bounds_of_interval(ii);
      errprecision x = static_cast<errprecision>(boost::math::float_next(intEndPoints.first));
      errprecision xtop = static_cast<errprecision>(boost::math::float_prior(intEndPoints.second));
      if ( double(xtop) > m_parent.m_max )
        break;
      std::pair<errprecision, errprecision> r =
        brent_find_minima(LookupTableErrorFunctor(impl.get()),x,xtop,bits,max_it);
      xstar = r.first; err = r.second;
      if( err < max_err ) {
        index = ii+1;
        max_err = err;
      }
    }

    /* want return to be 0 if the same, +/- on either side */
    max_err = -max_err;
    // std::cout << "stepSize: " << stepSize << " max_err: " << max_err << std::endl;
    return double(max_err-m_tol);
  }

private:

  UniformLookupTableGenerator<IN_TYPE,OUT_TYPE> &m_parent;
  std::string m_tableKey;
  double m_tol;

};

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*
  UniformLookupTableGenerator functions
*/

template <typename IN_TYPE, typename OUT_TYPE>
inline std::unique_ptr<UniformLookupTable<IN_TYPE,OUT_TYPE>> UniformLookupTableGenerator<IN_TYPE,OUT_TYPE>::generate_by_impl_size(
    std::string tableKey, unsigned long desiredSize)
{
  /* Use 2 query points to get relationship */
  const unsigned long N1  = 2;
  const double step1 = (m_max-m_min)/N1;
  const unsigned long N2  = 10;
  const double step2 = (m_max-m_min)/N2;

  UniformLookupTableParameters<IN_TYPE> par1;
  par1.minArg = m_min;
  par1.maxArg = m_max;
  par1.stepSize = step1;
  par1.transferFunction = mp_transfer_function;

  UniformLookupTableParameters<IN_TYPE> par2;
  par2.minArg = m_min;
  par2.maxArg = m_max;
  par2.stepSize = step2;
  par2.transferFunction = mp_transfer_function;

  std::unique_ptr<EvaluationImplementation<IN_TYPE,OUT_TYPE>> impl1 =
   UniformLookupTableFactory<IN_TYPE,OUT_TYPE>::Create(tableKey, mp_func_container, par1);
  std::unique_ptr<EvaluationImplementation<IN_TYPE,OUT_TYPE>> impl2 =
    UniformLookupTableFactory<IN_TYPE,OUT_TYPE>::Create(tableKey, mp_func_container, par2);

  unsigned long size1 = impl1->size();
  unsigned long size2 = impl2->size();

  if (size2 == size1) {
    throw "Query tables have same size.";
  }

  /* approximate step size for for desired impl size
     (assuming linear relationship of num_intervals to size */
  par1.stepSize = 1.0/((double)((N2-N1)*(desiredSize-size1)/(size2-size1) + N1));

  return UniformLookupTableFactory<IN_TYPE,OUT_TYPE>::Create(tableKey, mp_func_container, par1);
}

template <typename IN_TYPE, typename OUT_TYPE>
inline std::unique_ptr<UniformLookupTable<IN_TYPE,OUT_TYPE>> UniformLookupTableGenerator<IN_TYPE,OUT_TYPE>::generate_by_tol(std::string tableKey, double desiredTolerance)
{
  UniformLookupTableParameters<IN_TYPE> par;
  par.minArg = m_min;
  par.maxArg = m_max;
  par.stepSize =  (m_max-m_min)/1000.0;
  par.transferFunction = mp_transfer_function;
  /* generate a first approximation for the implementation */
  auto impl =
    UniformLookupTableFactory<IN_TYPE,OUT_TYPE>::Create(tableKey, mp_func_container, par);
  /* And initialize the functor used for refinement */
  OptimalStepSizeFunctor f(*this,tableKey,0);
  IN_TYPE stepSize  = impl->step_size();
  const int order  = impl->order();

  const double logTol = log(desiredTolerance);

  /*
    APPLY NEWTON'S METHOD IN log-log SPACE
    ('known' slope = order)

     Approximate the solution that satisfies tolerance based on the
     order of the implementation
     - assumes that the initial guess is in the asymptotic regime
       of the table's convergence
     - even if not in asymptotic regime, this thing has some pretty
       robust convergence!
  */
  const int  N_NEWTON_MAX_IT  = 0;     // max log-Newton-iterations
  const double NEWTON_IT_RTOL = 1e-5;
  const double NEWTON_IT_ATOL = 1e-10;

  int NEWTON_SUCCESS_FLAG = 0;
  std::vector<std::pair<double,double>> iterates;
  for (int iNewton = 0; iNewton < N_NEWTON_MAX_IT; iNewton++) {
    double err = f(stepSize);
    if ( fabs(err-desiredTolerance) <  fabs(desiredTolerance)*NEWTON_IT_RTOL+NEWTON_IT_ATOL ) {
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

    If a suitable bracket is found, this guarantees a solution below
    desiredTolerance. If not, bracket_and_solve throws.
  */
  const boost::uintmax_t BRACKET_MAX_IT = 50;  // Limit to maximum iterations.
  boost::uintmax_t it = BRACKET_MAX_IT;        // Initally our chosen max iterations, but updated with actual.

  /*
    Throw when the log-Newton method did not converge AND there are no
    bracketing iterations performed.
   */
  if ( !NEWTON_SUCCESS_FLAG && !BRACKET_MAX_IT) {
    std::cerr << "WARNING: No bracketing iterations specified." << std::endl;
    std::stringstream throwMessage;
    throwMessage << "log-Newton method did not converge in " << N_NEWTON_MAX_IT << " steps.";
    throw throwMessage;
  }
  /*
    Use the guess step size as an initialization to bracket_and_solve
  */
  using namespace boost::math::tools;

  IN_TYPE factor = 2;   // Mult/divide factor for bracket expansion when searching.
  bool is_rising = 1;  // The error curve should always be rising
  int digits = std::numeric_limits<double>::digits;  // Maximum possible binary digits accuracy for type T.
  // Some fraction of digits is used to control how accurate to try to make the result.
  int get_digits = digits-30;  // We have to have a non-zero interval at each step, so
                               // doesn't have to be so accurate
                               // maximum accuracy is digits - 1.  But we also have to
                               // allow for inaccuracy in f(x), otherwise the last few iterations
                               // just thrash around.
  eps_tolerance<double> tol(get_digits); // Set the tolerance.
  OptimalStepSizeFunctor g(*this,tableKey,desiredTolerance); // functor for solving log(E(h)/TOL)-1 = 0

  /*
    Run the bracket and solve algorithm
    - answer is taken to be the lower point of the bracketing interval
    - this guarantees a tolerance lower than desired
  */
  std::pair<IN_TYPE, IN_TYPE> r = bracket_and_solve_root(g, stepSize, factor, is_rising, tol, it);

  /* Finally, return the implementation with the desired stepSize*/
  par.stepSize = r.first;
  return UniformLookupTableFactory<IN_TYPE,OUT_TYPE>::Create(tableKey,mp_func_container,par);
}

template <typename IN_TYPE, typename OUT_TYPE>
inline double UniformLookupTableGenerator<IN_TYPE,OUT_TYPE>::error_at_step_size(
    std::string tableKey, IN_TYPE stepSize)
{
  /*
    Can be implemented in terms of the Functor used in solving for a specific
    tolerance, so that is reused.
  */
  OptimalStepSizeFunctor f(*this,tableKey,0);
  double err = f(stepSize);
  return err;
}

template <typename IN_TYPE, typename OUT_TYPE>
inline void UniformLookupTableGenerator<IN_TYPE,OUT_TYPE>::plot_implementation_at_step_size(
    std::string tableKey, IN_TYPE stepSize)
{
  /*
    Can be implemented in terms of the Functor used in solving for a specific
    tolerance, so that is reused.
  */
  UniformLookupTableParameters<IN_TYPE> par;
  par.minArg = m_min;
  par.maxArg = m_max; 
  par.stepSize = stepSize;
  par.transferFunction = mp_transfer_function;
  auto impl =
    UniformLookupTableFactory<IN_TYPE,OUT_TYPE>::Create(tableKey, mp_func_container, par);

  std::cout << "# x func impl" << std::endl;
  for (double x=impl->min_arg();
       x < impl->max_arg();
       x+=impl->step_size()/10 ) {
    std::cout << x << " " <<
      (impl->function())(x) << " " <<
      (*impl)(x) << std::endl;
  }
}

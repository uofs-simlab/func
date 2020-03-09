/* Implementation of UniformLookupTableGenerator */
#include "UniformLookupTableGenerator.hpp"
// #include "UniformTables.hpp"

#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/multiprecision/float128.hpp>
#include <boost/math/special_functions/next.hpp>

#include <limits>
#include <iostream>

using namespace boost::multiprecision;

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*
   Nested Functor used for computing error in a given lookup table
*/

struct UniformLookupTableGenerator::LookupTableErrorFunctor
{
  LookupTableErrorFunctor(UniformLookupTable* impl) : m_impl(impl) {}
  /* operator() always returns a negative value */
  float128 operator()(float128 const& x)
  {
    float128 f_value = float128((*(m_impl->function()))(double(x)));
    float128 lut_value = float128((*m_impl)(double(x)));
    return -float128(2.0) * fabs( (f_value - lut_value) ) /
      (fabs(f_value)+fabs(lut_value));
  }

private:

  UniformLookupTable *m_impl;
};

/*
  Nested Functor used for finding optimal stepsize that satisfies TOL
*/

struct UniformLookupTableGenerator::OptimalStepSizeFunctor
{
  OptimalStepSizeFunctor(UniformLookupTableGenerator &parent, std::string tableKey, double tol) :  m_parent(parent), m_tableKey(tableKey) , m_tol(tol)
  {
  }

  double operator()(double const& stepSize)
  {
    using namespace boost::math::tools;

    UniformLookupTableParameters par;
    par.minArg = m_parent.m_min;
		par.maxArg = m_parent.m_max;
		par.stepSize = stepSize;
    auto impl = UniformLookupTableFactory::Create(m_tableKey, m_parent.mp_func, par);

    boost::uintmax_t max_it = 20;

    float128 max_err = 0;
    float128 xstar, err;

    /* get number of binary bits in mantissa */
    int bits = std::numeric_limits<float128>::digits;

    double eps = std::numeric_limits<double>::epsilon();

    /*
      for each interval in the uniform table, compute the maximum error
      - be careful about the top most interval, it may reach beyond the
      table range due to rounding errors
    */
    unsigned index = 0;
    for(unsigned ii=0; ii<impl->num_intervals()-1; ii++){

      std::pair<double,double> intEndPoints = impl->arg_bounds_of_interval(ii);
      float128 x = float128(boost::math::float_next(intEndPoints.first));
      float128 xtop = float128(boost::math::float_prior(intEndPoints.second));
      if ( double(xtop) > m_parent.m_max )
	break;
      std::pair<float128, float128> r =
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

  UniformLookupTableGenerator &m_parent;
  std::string m_tableKey;
  double m_tol;

};

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*
  UniformLookupTableGenerator functions
*/
UniformLookupTableGenerator::UniformLookupTableGenerator(EvaluationFunctor<double,double> *func, double minArg, double maxArg) : mp_func(func), m_min(minArg), m_max(maxArg) {}


UniformLookupTableGenerator::~UniformLookupTableGenerator()
{
}

std::unique_ptr<UniformLookupTable> UniformLookupTableGenerator::generate_by_step(std::string tableKey, double stepSize)
{
  UniformLookupTableParameters par;
  par.minArg = m_min; 
  par.maxArg = m_max; 
  par.stepSize = stepSize;
  return UniformLookupTableFactory::Create(tableKey, mp_func, par);
}

std::unique_ptr<UniformLookupTable> UniformLookupTableGenerator::generate_by_impl_size(std::string tableKey, unsigned long desiredSize)
{
  /* Use 2 query points to get relationship */
  const unsigned long N1  = 2;
  const double step1 = (m_max-m_min)/N1;
  const unsigned long N2  = 10;
  const double step2 = (m_max-m_min)/N2;
  UniformLookupTableParameters par1;
  par1.minArg = m_min;
  par1.maxArg = m_max;
  par1.stepSize = step1;
  UniformLookupTableParameters par2;
  par2.minArg = m_min;
  par2.maxArg = m_max;
  par2.stepSize = step2;

  std::unique_ptr<EvaluationImplementation> impl1 =
   UniformLookupTableFactory::Create(tableKey, mp_func, par1);
  std::unique_ptr<EvaluationImplementation> impl2 =
    UniformLookupTableFactory::Create(tableKey, mp_func, par2);

  unsigned long size1 = impl1->size();
  unsigned long size2 = impl2->size();

  if (size2 == size1) {
    throw "Query tables have same size.";
  }

  /* approximate step size for for desired impl size
     (assuming linear relationship of num_intervals to size */
  par1.stepSize = 1.0/((double)((N2-N1)*(desiredSize-size1)/(size2-size1) + N1));

  return UniformLookupTableFactory::Create(tableKey, mp_func, par1);
}

std::unique_ptr<UniformLookupTable> UniformLookupTableGenerator::generate_by_tol(std::string tableKey, double desiredTolerance)
{
  UniformLookupTableParameters par;
  par.minArg = m_min;
  par.maxArg = m_max;
  par.stepSize =  (m_max-m_min)/1000.0;
  /* generate a first approximation for the implementation */
  auto impl =
    UniformLookupTableFactory::Create(tableKey, mp_func, par);
  /* And initialize the functor used for refinement */
  OptimalStepSizeFunctor f(*this,tableKey,0);
  double stepSize  = impl->step_size();
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

  double factor  = 2;   // Mult/divide factor for bracket expansion when searching.
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
  std::pair<double, double> r = bracket_and_solve_root(g, stepSize, factor, is_rising, tol, it);

  /* Finally, return the implementation with the desired stepSize*/
  par.stepSize = r.first;
  return UniformLookupTableFactory::Create(tableKey,mp_func,par);
}


double UniformLookupTableGenerator::error_at_step_size(std::string tableKey,
						       double stepSize)
{
  /*
    Can be implemented in terms of the Functor used in solving for a specific
    tolerance, so that is reused.
  */
  OptimalStepSizeFunctor f(*this,tableKey,0);
  double err = f(stepSize);
  return err;
}

void UniformLookupTableGenerator::plot_implementation_at_step_size(std::string tableKey,
								   double stepSize)
{
  /*
    Can be implemented in terms of the Functor used in solving for a specific
    tolerance, so that is reused.
  */
  UniformLookupTableParameters par;
  par.minArg = m_min;
  par.maxArg = m_max; 
  par.stepSize = stepSize;
  auto impl =
    UniformLookupTableFactory::Create(tableKey, mp_func, par);

  std::cout << "# x func impl" << std::endl;
  for (double x=impl->min_arg();
       x < impl->max_arg();
       x+=impl->step_size()/10 ) {
    std::cout << x << " " <<
      (*(impl->function()))(x) << " " <<
      (*impl)(x) << std::endl;
  }
}

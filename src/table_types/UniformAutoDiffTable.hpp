/*
  Intermediate abstract class for LUTs with uniformly spaced grid points
*/
#pragma once
#include "EvaluationImplementation.hpp"

#include <memory>
#include <map>
#include <vector>
#include <functional>
#include <boost/math/differentiation/autodiff.hpp>

using boost::math::differentiation::autodiff_fvar;

template <unsigned int N>
class UniformAutoDiffTable : public UniformLookupTable
{
  friend class func_wrapper;
protected:

  EvaluationFunctor<autodiff_fvar<double,N>, autodiff_fvar<double,N>> *mp_boost_func;

  std::unique_ptr<double[]> m_grid;  // pointers to grid and evaluation data
  // a LUT array needs to be provided by each implementation

  unsigned                  m_numIntervals;   // sizes of grid and evaluation data
  unsigned                  m_numTableEntries;

  double                    m_stepSize;    // fixed grid spacing (and its inverse)
  double                    m_stepSize_inv;

  double                    m_tableMaxArg; // > m_maxArg if (m_maxArg-m_minArg)/m_stepSize is non-integer

public:

  UniformAutoDiffTable(EvaluationFunctor<autodiff_fvar<double,N>,autodiff_fvar<double,N>> *func, UniformLookupTableParameters par);
  virtual ~UniformAutoDiffTable(){};

  /* public access of protected data */
  double step_size(){ return m_stepSize; };
  unsigned num_table_entries(){ return m_numTableEntries; };
  unsigned num_intervals(){ return m_numIntervals; };
  void print_details(std::ostream&) override;
  std::pair<double,double> arg_bounds_of_interval(unsigned);
  EvaluationFunctor<autodiff_fvar<double,N>, autodiff_fvar<double,N>> *fvar_function();
};

/* This can probably be removed when we move to std::function */
class FuncWrapper : public EvaluationFunctor<double,double>
{
  public:
    EvaluationFunctor<autodiff_fvar<double,N>,autodiff_fvar<double,N>> *fvar_EvalFunctor;
    FuncWrapper( *new_EvalFuntor){ fvar_EvalFunctor = new_EvalFunctor; }
    double operator()(double x) override { return (double) (*fvar_EvalFunctor)(x); }
};

/* ////////////////////////////////////////////////////////////////////////////
   Factory and registration management

   UniformLookupTableFactory: singleton class responsible for
   - creating Derived classes of UniformAutoDiffTable
   - keeping a registry of derived classes
   I ASSUME THERE SHOULD ONLY BE ONE OF THESE
//////////////////////////////////////////////////////////////////////////// */
class UniformLookupTableFactory
{
public:
  // Only ever hand out unique pointers
  static std::unique_ptr<UniformAutoDiffTable> Create(std::string name,
					       EvaluationFunctor<autodiff_fvar<double,N>,autodiff_fvar<double,N>> *f,
					       UniformLookupTableParameters par);
  // Actual registration function
  static void RegisterFactoryFunction(std::string name,
        std::function<UniformAutoDiffTable*(EvaluationFunctor<autodiff_fvar<double,N>,autodiff_fvar<double,N>>*,UniformLookupTableParameters)> classFactoryFunction);
  // Get all keys from the registry
  static std::vector<std::string> get_registry_keys(void);

private:
  // the actual registry is private to this class
  static std::map<std::string, std::function<UniformAutoDiffTable*(
			       EvaluationFunctor<autodiff_fvar<double,N>,autodiff_fvar<double,N>>*,UniformLookupTableParameters)>>& get_registry();
  // Do NOT implement copy methods
  UniformLookupTableFactory(){};
  UniformLookupTableFactory(UniformLookupTableFactory const& copy);
  UniformLookupTableFactory& operator=(UniformLookupTableFactory const& copy);
};

/*
  Helper class for registering Uniform LUT implementations

  NOTE: implementation defined in this header so that templates get
  instantiated in derived LUT class files
*/
template<class T>
class UniformLookupTableRegistrar {
public:
  UniformLookupTableRegistrar(std::string className)
  {
    UniformLookupTableFactory::RegisterFactoryFunction(className,
	    [](EvaluationFunctor<autodiff_fvar<double,N>,autodiff_fvar<double,N>> *f, UniformLookupTableParameters par) -> UniformAutoDiffTable * { return new T(f, par);});
  }
};

/*
   Macros for class registration:
   - REGISTER_ULUT goes inside class declaration in .hpp
   - REGISTER_ULUT_IMPL goes inside .cpp implementation
*/
#define REGISTER_ULUT(classname) \
private: \
   static const UniformLookupTableRegistrar<classname> registrar;
#define STR_EXPAND(x...) #x
#define STR(x...) STR_EXPAND(x)
#define REGISTER_ULUT_IMPL(classname...) \
   const UniformLookupTableRegistrar<classname> classname::registrar(STR(classname));

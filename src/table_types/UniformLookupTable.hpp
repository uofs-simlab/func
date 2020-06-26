/*
  Intermediate abstract class for LUTs with uniformly spaced grid points
*/
#pragma once
#include "EvaluationImplementation.hpp"
#include <boost/math/differentiation/autodiff.hpp>

#include <memory>
#include <map>
#include <vector>
#include <functional>

struct UniformLookupTableParameters
{
  // UniformLookupTableParameters(): minArg(0), maxArg(1), stepSize(1) {}
  double minArg = 0;
  double maxArg = 1;
  double stepSize = 1;
};

template <unsigned int NUM_COEFS, std::size_t ALIGN=64>
struct alignas(ALIGN) polynomial{
  double coefs[NUM_COEFS];
};

class UniformLookupTable : public EvaluationImplementation
{
protected:

  std::unique_ptr<double[]> m_grid;  // pointers to grid and evaluation data
  // a LUT array needs to be provided by each implementation

  unsigned                  m_numIntervals;   // sizes of grid and evaluation data
  unsigned                  m_numTableEntries;

  double                    m_stepSize;    // fixed grid spacing (and its inverse)
  double                    m_stepSize_inv;

  double                    m_tableMaxArg; // > m_maxArg if (m_maxArg-m_minArg)/m_stepSize is non-integer

public:

  UniformLookupTable(EvaluationFunctor<double,double> *func, UniformLookupTableParameters par);
  virtual ~UniformLookupTable(){};

  /* public access of protected data */
  double step_size(){ return m_stepSize; };
  unsigned num_table_entries(){ return m_numTableEntries; };
  unsigned num_intervals(){ return m_numIntervals; };
  void print_details(std::ostream&) override;
  std::pair<double,double> arg_bounds_of_interval(unsigned);
};

/* ////////////////////////////////////////////////////////////////////////////
   Factory and registration management

   UniformLookupTableFactory: singleton class responsible for
   - creating Derived classes of UniformLookupTable
   - keeping a registry of derived classes
//////////////////////////////////////////////////////////////////////////// */
class UniformLookupTableFactory
{
public:
  // Only ever hand out unique pointers
  static std::unique_ptr<UniformLookupTable> Create(std::string name,
                          EvaluationFunctor<func_type,func_type> *f, // make this a data type with 8 functions in it
                          UniformLookupTableParameters par)
  {
    // Create a UniformLookupTable
    UniformLookupTable * instance = nullptr;

    // find the name in the registry and call factory method.
    auto it = get_registry().find(name);
    if(it != get_registry().end()){
      // find the index of the function in the function map
      instance = it->second(f,par);
    }

    // wrap instance in a unique ptr and return (if created)
    if(instance == nullptr)
      throw "Table name not found in registry."; // TODO better exception
    return std::unique_ptr<UniformLookupTable>(instance);

  }

  // Actual registration function
  static void RegisterFactoryFunction(std::string name,
        std::function<UniformLookupTable*(EvaluationFunctor<func_type,func_type>*,UniformLookupTableParameters)> classFactoryFunction)
  {
    // register a derived class factory function
    get_registry()[name] = classFactoryFunction;
  }

  // Get all keys from the registry
  static std::vector<std::string> get_registry_keys(void)
  {
    // copy all keys from the registry map into a vector
    std::vector<std::string> keys;
    for (auto const& elem : get_registry() ) {
      keys.push_back(elem.first);
    }
    return keys;
  }

private:
  // the actual registry is private to this class
  static std::map<std::string, std::function<UniformLookupTable*(
			       EvaluationFunctor<func_type,func_type>*,UniformLookupTableParameters)>>& get_registry()
  {
    // Get the singleton instance of the registry map
    static std::map<std::string, std::function<UniformLookupTable*(
                     EvaluationFunctor<func_type,func_type>*,UniformLookupTableParameters)>> registry;
    return registry;
  }

  // Do NOT implement copy methods
  UniformLookupTableFactory(){};
  UniformLookupTableFactory(UniformLookupTableFactory const& copy);
  UniformLookupTableFactory& operator=(UniformLookupTableFactory const& copy);
};

/*
  Helper class for registering Uniform LUT implementations.
  T is a class name and N is the number of differentiations T needs
  in order to be built.

  NOTE: implementation defined in this header so that templates get
  instantiated in derived LUT class files
*/
template <class T, unsigned int N = 0>
class UniformLookupTableRegistrar {
public:
  UniformLookupTableRegistrar(std::string className)
  {
    if(N==0){
      UniformLookupTableFactory::RegisterFactoryFunction(className,
          [](EvaluationFunctor<double,double> *f, UniformLookupTableParameters par) -> UniformLookupTable * { return new T(f, par);});
    }else{
      UniformLookupTableFactory::RegisterFactoryFunction(className,
          [](EvaluationFunctor<autodiff_fvar<double,N>,autodiff_fvar<double,N>> *f, UniformLookupTableParameters par) 
          -> UniformLookupTable * { return new T(f, par);});
    }
    // set num_derivs map
    UniformLookupTableFactory::RegisterFactoryFunctionDerivs(className,...);
  }
};

/*
   Macros for class registration:
   - REGISTER_ULUT goes inside class declaration in .hpp
   - REGISTER_ULUT_IMPL goes inside .cpp implementation
   - *_DIFF are variants of these macros for registering the
   differentiation tables (mind the order of args)
*/
#define REGISTER_ULUT(classname) \
private: \
   static const UniformLookupTableRegistrar<classname> registrar;
#define STR_EXPAND(x...) #x
#define STR(x...) STR_EXPAND(x)
#define REGISTER_ULUT_IMPL(classname...) \
   const UniformLookupTableRegistrar<classname> classname::registrar(STR(classname));

#define REGISTER_ULUT_DIFF(order,classname) \
private: \
   static const UniformLookupTableRegistrar<classname,autodiff_fvar<double,order>> registrar;
#define STR_EXPAND(x...) #x
#define STR(x...) STR_EXPAND(x)
#define REGISTER_ULUT_IMPL_DIFF(order,classname...) \
   const UniformLookupTableRegistrar<classname,autodiff_fvar<double,order>> classname::registrar(STR(classname));

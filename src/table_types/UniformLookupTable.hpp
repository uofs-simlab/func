/*
  Intermediate abstract class for LUTs with uniformly spaced grid points
  TODO use "using" statements instead of putting "this->" everywhere
*/
#pragma once
#include "FunctionContainer.hpp"
#include "EvaluationImplementation.hpp"
#include "TransferFunction.hpp"

#include <memory>
#include <map>
#include <vector>
#include <functional>

template <typename IN_TYPE>
struct UniformLookupTableParameters
{
  IN_TYPE minArg = 0;
  IN_TYPE maxArg = 1;
  IN_TYPE stepSize = 1;
  std::shared_ptr<TransferFunction<IN_TYPE>> transferFunction;
};

template <typename OUT_TYPE, unsigned int NUM_COEFS, std::size_t ALIGN=64>
struct alignas(ALIGN) polynomial{
  OUT_TYPE coefs[NUM_COEFS];
};

template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE>
class UniformLookupTable : public EvaluationImplementation<IN_TYPE,OUT_TYPE>
{
protected:

  std::unique_ptr<IN_TYPE[]> m_grid;  // pointers to grid and evaluation data
  // a LUT array needs to be provided by each implementation

  unsigned int m_numIntervals;   // sizes of grid and evaluation data
  unsigned int m_numTableEntries;

  IN_TYPE  m_stepSize;    // fixed grid spacing (and its inverse)
  IN_TYPE  m_stepSize_inv;

  IN_TYPE  m_tableMaxArg; // > m_maxArg if (m_maxArg-m_minArg)/m_stepSize is non-integer

public:

  UniformLookupTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container,
      UniformLookupTableParameters<IN_TYPE> par) :
    EvaluationImplementation<IN_TYPE,OUT_TYPE>(func_container->standard_func, "uniform_lookup_table")
  {
    /* Base class variables */
    this->m_minArg = par.minArg; this->m_maxArg = par.maxArg;

    /*
       Uniform Lookup Table variables
       - corresponds to the grid
    */
    m_stepSize     = par.stepSize;
    m_stepSize_inv = 1.0/m_stepSize;
    m_numIntervals = (unsigned) ceil(m_stepSize_inv*(this->m_maxArg-this->m_minArg))+1;
    /*
      If the step size does not exactly divide the arg domain, the max
      arg of the table is set to the nearest value above such that it
      does.
     */
    m_tableMaxArg = this->m_maxArg;
    if ( m_tableMaxArg < this->m_minArg+this->m_stepSize*(m_numIntervals-1) )
      m_tableMaxArg = this->m_minArg+m_stepSize*(m_numIntervals-1);
    m_grid.reset(new IN_TYPE[m_numIntervals]);
  }

  virtual ~UniformLookupTable(){};

  /* public access of protected data */
  IN_TYPE step_size(){ return m_stepSize; };
  unsigned num_table_entries(){ return m_numTableEntries; };
  unsigned num_intervals(){ return m_numIntervals; };

  void print_details(std::ostream &out)
  {
    out << this->m_name << " " << this->m_minArg << " " << this->m_maxArg << " "
        << m_stepSize << " " << m_numIntervals << " ";
  }

  std::pair<IN_TYPE,IN_TYPE> arg_bounds_of_interval(unsigned intervalNumber)
  {
    return std::make_pair(m_grid[intervalNumber],m_grid[intervalNumber+1]);
  }
};

/* ////////////////////////////////////////////////////////////////////////////
   Factory and registration management

   UniformLookupTableFactory: singleton class responsible for
   - creating Derived classes of UniformLookupTable
   - keeping a registry of derived classes
//////////////////////////////////////////////////////////////////////////// */
template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE>
class UniformLookupTableFactory
{
public:
  // Only ever hand out unique pointers
  static std::unique_ptr<UniformLookupTable<IN_TYPE,OUT_TYPE>> Create(
      std::string name, FunctionContainer<IN_TYPE,OUT_TYPE> *fc,
      UniformLookupTableParameters<IN_TYPE> par)
  {
    // Create a UniformLookupTable
    UniformLookupTable<IN_TYPE,OUT_TYPE> * instance = nullptr;

    // find the name in the registry and call factory method.
    auto it = get_registry().find(name);
    if(it != get_registry().end())
      instance = it->second(fc,par);

    // wrap instance in a unique ptr and return (if created)
    if(instance == nullptr)
      throw std::invalid_argument(name + " not found in registry.");
    return std::unique_ptr<UniformLookupTable<IN_TYPE,OUT_TYPE>>(instance);
  }

  // Actual registration function
  static void RegisterFactoryFunction(std::string name,
      std::function<UniformLookupTable<IN_TYPE,OUT_TYPE>*(
        FunctionContainer<IN_TYPE,OUT_TYPE>*, 
        UniformLookupTableParameters<IN_TYPE>
      )> classFactoryFunction)
  {
    // register a derived class factory function
    get_registry()[name] = classFactoryFunction;
  }

  // Get all keys from the registry
  static std::vector<std::string> get_registry_keys(void)
  {
    // copy all keys from the registry map into a vector
    std::vector<std::string> keys;
    for(auto const& elem : get_registry())
      keys.push_back(elem.first);
    return keys;
  }

private:
  // the actual registry is private to this class
  static std::map<std::string, 
    std::function<UniformLookupTable<IN_TYPE,OUT_TYPE>*(
			FunctionContainer<IN_TYPE,OUT_TYPE>*,
      UniformLookupTableParameters<IN_TYPE>
    )>>& get_registry()
  {
    // Get the singleton instance of the registry map
    static std::map<std::string, std::function<UniformLookupTable<IN_TYPE,OUT_TYPE>*(
               FunctionContainer<IN_TYPE,OUT_TYPE>*,UniformLookupTableParameters<IN_TYPE>)>> registry;
    return registry;
  }

  // Do NOT implement copy methods
  UniformLookupTableFactory(){};
  UniformLookupTableFactory(UniformLookupTableFactory<IN_TYPE,OUT_TYPE> const& copy);
  UniformLookupTableFactory& operator=(UniformLookupTableFactory<IN_TYPE,OUT_TYPE> const& copy);
};

/*
  Helper class for registering Uniform LUT implementations.

  NOTE: implementation defined in this header so that templates get
  instantiated in derived LUT class files
*/
template <typename IN_TYPE, typename OUT_TYPE, class T>
class UniformLookupTableRegistrar {
public:
  UniformLookupTableRegistrar(std::string className)
  {
    UniformLookupTableFactory<IN_TYPE,OUT_TYPE>::RegisterFactoryFunction(className,
      [](FunctionContainer<IN_TYPE,OUT_TYPE> *fc,
         UniformLookupTableParameters<IN_TYPE> par
      ) -> UniformLookupTable<IN_TYPE,OUT_TYPE>* { return new T(fc, par); }
    );
  }
};

/*
   Macros for class registration:
   - REGISTER_LUT goes inside class definition
   - TODO Currently assumes template names are IN_TYPE and OUT_TYPE
   which isn't great but it's a pain to generalize
*/
#define REGISTER_LUT(classname) \
private: \
   static const UniformLookupTableRegistrar<IN_TYPE,OUT_TYPE,classname> registrar;
#define STR_EXPAND(x...) #x
#define STR(x...) STR_EXPAND(x)
#define REGISTER_LUT_IMPL(classname) \
   const UniformLookupTableRegistrar<IN_TYPE,OUT_TYPE,classname> classname::registrar(STR(classname));
//
//#define STR_EXPAND(x...) #x
//#define STR(x...) STR_EXPAND(x)
//#define REGISTER_LUT(classname) \
//private: \
//   static const UniformLookupTableRegistrar<IN_TYPE,OUT_TYPE,classname> registrar(STR(classname));

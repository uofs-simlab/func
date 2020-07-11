/*
  Intermediate abstract class for LUTs with uniformly spaced grid points
*/
#pragma once
#include "FunctionContainer.hpp"
#include "EvaluationImplementation.hpp"

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

  UniformLookupTable(FunctionContainer *func_container, UniformLookupTableParameters par);
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
                          FunctionContainer *fc,
                          UniformLookupTableParameters par);
  // Actual registration function
  static void RegisterFactoryFunction(std::string name,
        std::function<UniformLookupTable*(FunctionContainer*,UniformLookupTableParameters)> classFactoryFunction);
  // Get all keys from the registry
  static std::vector<std::string> get_registry_keys(void);

private:
  // the actual registry is private to this class
  static std::map<std::string, std::function<UniformLookupTable*(
			       FunctionContainer*,UniformLookupTableParameters)>>& get_registry();
  // Do NOT implement copy methods
  UniformLookupTableFactory(){};
  UniformLookupTableFactory(UniformLookupTableFactory const& copy);
  UniformLookupTableFactory& operator=(UniformLookupTableFactory const& copy);
};

/*
  Helper class for registering Uniform LUT implementations.

  NOTE: implementation defined in this header so that templates get
  instantiated in derived LUT class files
*/
template <class T>
class UniformLookupTableRegistrar {
public:
  UniformLookupTableRegistrar(std::string className)
  {
    UniformLookupTableFactory::RegisterFactoryFunction(className,
        [](FunctionContainer *fc, UniformLookupTableParameters par) -> UniformLookupTable * { return new T(fc, par);});
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

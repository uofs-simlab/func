/*
  TODO explain how this links up with MetaTable & each other implementation
  Intermediate abstract class for LUTs with uniformly spaced grid points.
  Outfits tables with enough tools to sample from a nonuniform
  grid efficiently.

  Notes:
  - In the case where (max-min)/stepsize is not an integer then
  the real table max is greater than the user supplied max
*/
#pragma once
#include "FunctionContainer.hpp"
#include "EvaluationImplementation.hpp"
#include "json.hpp"
#include "config.hpp" // FUNC_USE_SMALL_REGISTRY

#include <map>
#include <memory>
#include <vector>
#include <functional>
#include <fstream>

template <typename TIN>
struct UniformLookupTableParameters
{
  TIN minArg;
  TIN maxArg;
  TIN stepSize;

  // support initializer lists
  // still some unfortunate redundancy b/c std::string also takes a templated
  // initializer list...
  UniformLookupTableParameters(TIN min, TIN max, TIN step) :
    minArg(min), maxArg(max), stepSize(step) {}
  UniformLookupTableParameters(){}
};

static constexpr unsigned int alignments[] = {0,1,2,4,4,8,8,8,8};

template <typename TOUT, unsigned int N>
struct alignas(sizeof(TOUT)*alignments[N]) polynomial{
  static const unsigned int num_coefs = N;
  TOUT coefs[N];
};

/* Use to inherit UniformLookupTable's member variables 
   without having to put "this->" everywhere. */
#define INHERIT_UNIFORM_LUT(TIN,TOUT) \
  using UniformLookupTable<TIN,TOUT>::m_grid; \
  using UniformLookupTable<TIN,TOUT>::m_numIntervals; \
  using UniformLookupTable<TIN,TOUT>::m_numTableEntries; \
  using UniformLookupTable<TIN,TOUT>::m_stepSize; \
  using UniformLookupTable<TIN,TOUT>::m_stepSize_inv; \
  using UniformLookupTable<TIN,TOUT>::m_tableMaxArg

template <typename TIN, typename TOUT = TIN>
class UniformLookupTable : public EvaluationImplementation<TIN,TOUT>
{
protected:
  INHERIT_EVALUATION_IMPL(TIN,TOUT);

  // TODO discuss having m_grid exclusive to NonUniformTables
  std::unique_ptr<TIN[]> m_grid;  // pointers to grid and evaluation data
  // a polynomial (above) array needs to be provided by each implementation

  // get the ith polynomial's jth coefficient
  virtual TOUT get_table_entry(unsigned int i, unsigned int j)=0;
  // get the number of coefficients used
  virtual unsigned int get_num_coefs()=0;

  unsigned int m_numIntervals;   // sizes of grid and evaluation data
  unsigned int m_numTableEntries;

  TIN  m_stepSize;    // fixed grid spacing (and its inverse)
  TIN  m_stepSize_inv;

  TIN  m_tableMaxArg; // > m_maxArg if (m_maxArg-m_minArg)/m_stepSize is non-integer

public:

  UniformLookupTable(FunctionContainer<TIN,TOUT> *func_container,
      UniformLookupTableParameters<TIN> par) :
    EvaluationImplementation<TIN,TOUT>(func_container->standard_func, "uniform_lookup_table")
  {
    /* Base class variables */
    m_minArg = par.minArg; m_maxArg = par.maxArg;

    /*
       Uniform Lookup Table variables
       - corresponds to the grid
    */
    m_stepSize     = par.stepSize;
    m_stepSize_inv = 1.0/m_stepSize;
    m_numIntervals = (unsigned) ceil(m_stepSize_inv*(m_maxArg-m_minArg))+1;
    /*
      If the step size does not exactly divide the arg domain, the max
      arg of the table is set to the nearest value above such that it
      does.
     */
    m_tableMaxArg = m_maxArg;
    if ( m_tableMaxArg < m_minArg+m_stepSize*(m_numIntervals-1) )
      m_tableMaxArg = m_minArg+m_stepSize*(m_numIntervals-1);
    m_grid.reset(new TIN[m_numIntervals]);
  }

  /* Set every generic member variable from a json file */
  UniformLookupTable(FunctionContainer<TIN,TOUT> *func_container,
      std::string filename) :
    EvaluationImplementation<TIN,TOUT>(func_container->standard_func, "uniform_lookup_table")
  {
    // Assuming we're still using the same function container as before
    // We'll end up making 2 jsonStats objects because the MetaTable needs
    // to check the table name and construct m_table. Not a big deal
    std::ifstream file_reader(filename);
    using nlohmann::json;
    json jsonStats;
    file_reader >> jsonStats;

    // name checking happens in MetaTable
    m_minArg = jsonStats["minArg"].get<TIN>();
    m_maxArg = jsonStats["maxArg"].get<TIN>();
    m_stepSize = jsonStats["stepSize"].get<TIN>();
    m_stepSize_inv = 1.0/m_stepSize;
    m_order = jsonStats["order"].get<unsigned int>();
    m_dataSize = jsonStats["dataSize"].get<unsigned int>();
    m_numIntervals = jsonStats["numIntervals"].get<unsigned int>();
    m_numTableEntries = jsonStats["numTableEntries"].get<unsigned int>();

    // recompute m_tableMaxArg
    m_tableMaxArg = m_maxArg;
    if ( m_tableMaxArg < m_minArg+m_stepSize*(m_numIntervals-1) )
      m_tableMaxArg = m_minArg+m_stepSize*(m_numIntervals-1);

    // Not recomputing m_grid so the NonuniformLUTs can use this code.
    // the array of polynomials will now be built by the MetaTable
  }

  virtual ~UniformLookupTable(){};

  /* public access of protected data */
  TIN step_size(){ return m_stepSize; };
  unsigned num_table_entries(){ return m_numTableEntries; };
  unsigned num_intervals(){ return m_numIntervals; };

  void print_details(std::ostream &out)
  {
    out << m_name << " " << m_minArg << " " << m_maxArg << " "
        << m_stepSize << " " << m_numIntervals << " ";
  }

  virtual std::pair<TIN,TIN> arg_bounds_of_interval(unsigned intervalNumber)
  {
    return std::make_pair(m_minArg + intervalNumber*m_stepSize,m_minArg + (intervalNumber+1)*m_stepSize);
  }

  /* Write table data to the provided ostream in the form of json. Virtual in case derived
     classes need to save their own member vars. eg. NonuniformLUTs will store their
     grid and rebuild their transfer function */
  virtual void print_details_json(std::ostream& out)
  {
    // TODO if FunC gets a namespace we should consider renaming our json
    // functions to to_json() and from_json() functions as seen here
    // https://github.com/nlohmann/json
    using nlohmann::json;
    json jsonStats;

    jsonStats["_comment"] = "FunC lookup table data";
    jsonStats["name"] = m_name;
    jsonStats["minArg"] = m_minArg;
    jsonStats["maxArg"] = m_maxArg;
    jsonStats["stepSize"] = m_stepSize;
    jsonStats["order"] = m_order;
    jsonStats["dataSize"] = m_dataSize;
    jsonStats["numIntervals"] = m_numIntervals;
    jsonStats["numTableEntries"] = m_numTableEntries;

    // save the polynomial coefs of each lookup table
    // Note: m_order is used as the number of polynomial coefs
    for(unsigned int i=0; i<m_numTableEntries; i++)
      for(unsigned int j=0; j<get_num_coefs(); j++)
        jsonStats["table"][std::to_string(i)]["coefs"][std::to_string(j)] = get_table_entry(i,j);

    out << jsonStats.dump(2) << std::endl;
  }
};

/* ////////////////////////////////////////////////////////////////////////////
  Factory and registration management

  UniformLookupTableFactory: singleton class responsible for
  - creating Derived classes of UniformLookupTable
  - keeping a registry of derived classes
  - It's usage throughout FunC implies 8 namespaces are made.
  UniformLookupTableFactory<double>::
  UniformLookupTableFactory<float>::
  UniformLookupTableFactory<double,float>::
  UniformLookupTableFactory<float,double>::
  UniformLookupTableFactory<double,double,std::string>::
  UniformLookupTableFactory<float,float,std::string>::
  UniformLookupTableFactory<double,float,std::string>::
  UniformLookupTableFactory<float,double,std::string>::
  and NonUniformTables just add to the pile
  TODO desperately need a dev flag for this

  Usage example:
  // given the fc is an initialized function container
  // and par is an initialized UniformLookupTableParameters
  std::unique_ptr<EvaluationImplementation<double>> p_table = 
    UniformLookupTableFactory<double>::Create("UniformCubicPrecomputedInterpolationTable",fc, par);

//////////////////////////////////////////////////////////////////////////// */
template <typename TIN, typename TOUT = TIN, class OTHER = UniformLookupTableParameters<TIN>>
class UniformLookupTableFactory
{
public:
  // Only ever hand out unique pointers
  static std::unique_ptr<UniformLookupTable<TIN,TOUT>> Create(
      std::string name, FunctionContainer<TIN,TOUT> *fc,
      OTHER args)
  {
    // Create a UniformLookupTable
    UniformLookupTable<TIN,TOUT> *instance = nullptr;

    // find the name in the registry and call factory method.
    auto it = get_registry().find(name);
    if(it != get_registry().end())
      instance = it->second(fc, args);

    // wrap instance in a unique ptr and return (if created)
    if(instance == nullptr)
      throw std::invalid_argument(name + " not found in registry.");
    return std::unique_ptr<UniformLookupTable<TIN,TOUT>>(instance);
  }

  // Actual registration function
  static void RegisterFactoryFunction(std::string name,
      std::function<UniformLookupTable<TIN,TOUT>*(
        FunctionContainer<TIN,TOUT>*, 
        OTHER
      )> classFactoryFunction)
  {
    // register a derived class factory function
    get_registry()[name] = classFactoryFunction;
  }

  // Get all keys from the registry
  static std::vector<std::string> get_registry_keys()
  {
    // copy all keys from the registry map into a vector
    std::vector<std::string> keys;
    for(auto const& elem : get_registry())
      keys.push_back(elem.first);
    return keys;
  }

private:
  // the actual registry is private to this class
  static std::map<std::string, std::function<UniformLookupTable<TIN,TOUT>*(
			FunctionContainer<TIN,TOUT>*, OTHER)>>& get_registry()
  {
    // Get the singleton instance of the registry map
    static std::map<std::string, std::function<UniformLookupTable<TIN,TOUT>*(
               FunctionContainer<TIN,TOUT>*,OTHER)>> registry;
    return registry;
  }

  // Do NOT implement copy methods
  UniformLookupTableFactory(){};
  UniformLookupTableFactory(UniformLookupTableFactory<TIN,TOUT> const& copy);
  UniformLookupTableFactory& operator=(UniformLookupTableFactory<TIN,TOUT> const& copy);
};

/*
  Helper class for registering Uniform LUT implementations.

  NOTE: implementation defined in this header so that templates get
  instantiated in derived LUT class files
*/
template <class T, typename TIN, typename TOUT, class OTHER = UniformLookupTableParameters<TIN>>
class UniformLookupTableRegistrar {
public:
  UniformLookupTableRegistrar(std::string className)
  {
    UniformLookupTableFactory<TIN,TOUT>::RegisterFactoryFunction(className,
      [](FunctionContainer<TIN,TOUT> *fc,
         OTHER args
      ) -> UniformLookupTable<TIN,TOUT>* { return new T(fc, args); }
    );
  }
};

/*
   Macros for class registration:
   - FUNC_REGISTER_LUT goes inside class definition
   - FUNC_REGISTER_EACH_ULUT_IMPL goes underneath the class definition.
   Several different versions of this macro exist for registering templated classes
   - other... is for template parameters unrelated to the tables TIN and TOUT
*/
// this macro is used for nonuniform tables
#define FUNC_REGISTER_LUT(classname) \
private: \
  static const UniformLookupTableRegistrar<classname,TIN,TOUT> registrar; \
  static const UniformLookupTableRegistrar<classname,TIN,TOUT,std::string> str_registrar

#define FUNC_STR_EXPAND(x...) #x
#define FUNC_STR(x...) FUNC_STR_EXPAND(x)

// everything after this point is specialized to uniform LUTs.
#define FUNC_REGISTER_ULUT_IMPL(classname,TIN,TOUT) \
  template<> const \
  UniformLookupTableRegistrar<classname<TIN,TOUT>,TIN,TOUT> \
    classname<TIN,TOUT>::registrar(FUNC_STR(classname))

#define FUNC_REGISTER_TEMPLATED_ULUT_IMPL(classname,TIN,TOUT,other...) \
  template<> const \
    UniformLookupTableRegistrar<classname<TIN,TOUT,other>,TIN,TOUT> \
    classname<TIN,TOUT,other>::registrar(FUNC_STR(classname<other>))

#ifndef FUNC_USE_SMALL_REGISTRY
// macros used in each class to quickly register 4 different template instantiations
#define FUNC_REGISTER_EACH_ULUT_IMPL(classname)\
  FUNC_REGISTER_ULUT_IMPL(classname,double,double); \
  FUNC_REGISTER_ULUT_IMPL(classname,float,double);  \
  FUNC_REGISTER_ULUT_IMPL(classname,double,float);  \
  FUNC_REGISTER_ULUT_IMPL(classname,float,float);

#define FUNC_REGISTER_EACH_TEMPLATED_ULUT_IMPL(classname,other...) \
  FUNC_REGISTER_TEMPLATED_ULUT_IMPL(classname,double,double,other); \
  FUNC_REGISTER_TEMPLATED_ULUT_IMPL(classname,float,double,other);  \
  FUNC_REGISTER_TEMPLATED_ULUT_IMPL(classname,double,float,other);  \
  FUNC_REGISTER_TEMPLATED_ULUT_IMPL(classname,float,float,other)

#else // just build double -> double tables

#define FUNC_REGISTER_EACH_ULUT_IMPL(classname)\
  FUNC_REGISTER_ULUT_IMPL(classname,double,double); \

#define FUNC_REGISTER_EACH_TEMPLATED_ULUT_IMPL(classname,other...) \
  FUNC_REGISTER_TEMPLATED_ULUT_IMPL(classname,double,double,other); \

#endif // FUNC_USE_SMALL_REGISTRY

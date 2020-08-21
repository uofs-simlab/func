/*
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

#include <map>
#include <memory>
#include <vector>
#include <functional>
#include <fstream>

template <typename IN_TYPE>
struct UniformLookupTableParameters
{
  IN_TYPE minArg;
  IN_TYPE maxArg;
  IN_TYPE stepSize;

  // support initializer lists
  UniformLookupTableParameters(IN_TYPE min, IN_TYPE max, IN_TYPE step) :
    minArg(min), maxArg(max), stepSize(step) {}
  UniformLookupTableParameters(){}
};

static constexpr unsigned int alignments[] = {0,1,2,4,4,8,8,8,8};

template <typename OUT_TYPE, unsigned int N>
struct alignas(sizeof(OUT_TYPE)*alignments[N]) polynomial{
  static const unsigned int num_coefs = N;
  OUT_TYPE coefs[N];
};

/* Use to inherit UniformLookupTable's member variables 
   without having to put "this->" everywhere. */
#define INHERIT_UNIFORM_LUT(IN_TYPE,OUT_TYPE) \
  using UniformLookupTable<IN_TYPE,OUT_TYPE>::m_grid; \
  using UniformLookupTable<IN_TYPE,OUT_TYPE>::m_numIntervals; \
  using UniformLookupTable<IN_TYPE,OUT_TYPE>::m_numTableEntries; \
  using UniformLookupTable<IN_TYPE,OUT_TYPE>::m_stepSize; \
  using UniformLookupTable<IN_TYPE,OUT_TYPE>::m_stepSize_inv; \
  using UniformLookupTable<IN_TYPE,OUT_TYPE>::m_tableMaxArg

template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE>
class UniformLookupTable : public EvaluationImplementation<IN_TYPE,OUT_TYPE>
{
protected:
  INHERIT_EVALUATION_IMPL(IN_TYPE,OUT_TYPE);
  std::unique_ptr<IN_TYPE[]> m_grid;  // pointers to grid and evaluation data
  // a polynomial (above) array needs to be provided by each implementation

  // get the ith polynomial's jth coefficient
  virtual OUT_TYPE get_table_entry(unsigned int i, unsigned int j)=0;
  // get the number of coefficients used
  virtual unsigned int get_num_coefs()=0;

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
    m_grid.reset(new IN_TYPE[m_numIntervals]);
  }

  /* Set every generic member variable from a json file */
  UniformLookupTable(FunctionContainer<IN_TYPE,OUT_TYPE> *func_container,
      std::string filename) :
    EvaluationImplementation<IN_TYPE,OUT_TYPE>(func_container->standard_func, "uniform_lookup_table")
  {
    // Assuming we're still using the same function container as before
    // Also kinda inefficient since the derived class will also make its own
    // jsonStats but oh well
    std::ifstream file_reader(filename);
    using nlohmann::json;
    json jsonStats;
    file_reader >> jsonStats;
    m_name = jsonStats["name"].get<std::string>();
    m_minArg = jsonStats["minArg"].get<IN_TYPE>();
    m_maxArg = jsonStats["maxArg"].get<IN_TYPE>();
    m_stepSize = jsonStats["stepSize"].get<IN_TYPE>();
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
    // the array of polynomials must now be built by each table implementation
  }

  virtual ~UniformLookupTable(){};

  /* public access of protected data */
  IN_TYPE step_size(){ return m_stepSize; };
  unsigned num_table_entries(){ return m_numTableEntries; };
  unsigned num_intervals(){ return m_numIntervals; };

  void print_details(std::ostream &out)
  {
    out << m_name << " " << m_minArg << " " << m_maxArg << " "
        << m_stepSize << " " << m_numIntervals << " ";
  }

  virtual std::pair<IN_TYPE,IN_TYPE> arg_bounds_of_interval(unsigned intervalNumber)
  {
    return std::make_pair(m_minArg + intervalNumber*m_stepSize,m_minArg + (intervalNumber+1)*m_stepSize);
  }


  /* Write table data to the provided ostream in the form of json. Virtual in case derived
     classes need to save their own member vars. eg. NonuniformLUTs will store their
     grid and rebuild their transfer function */
  virtual void print_details_json(std::ostream& out)
  {
    // TODO if FunC gets a namespace we should consider renaming our json
    // functions to to_json() and from_json functions as seen here
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
  TODO desperately need a dev flag for this

  Usage example:
  // given the fc is an initialized function container
  // and par is an initialized UniformLookupTableParameters
  std::unique_ptr<EvaluationImplementation<double>> p_table = 
    UniformLookupTableFactory<double>::Create("UniformCubicPrecomputedInterpolationTable",fc, par);

//////////////////////////////////////////////////////////////////////////// */
template <typename IN_TYPE, typename OUT_TYPE = IN_TYPE, class OTHER = UniformLookupTableParameters<IN_TYPE>>
class UniformLookupTableFactory
{
public:
  // Only ever hand out unique pointers
  static std::unique_ptr<UniformLookupTable<IN_TYPE,OUT_TYPE>> Create(
      std::string name, FunctionContainer<IN_TYPE,OUT_TYPE> *fc,
      OTHER args)
  {
    // Create a UniformLookupTable
    UniformLookupTable<IN_TYPE,OUT_TYPE> *instance = nullptr;

    // find the name in the registry and call factory method.
    auto it = get_registry().find(name);
    if(it != get_registry().end())
      instance = it->second(fc, args);

    // wrap instance in a unique ptr and return (if created)
    if(instance == nullptr)
      throw std::invalid_argument(name + " not found in registry.");
    return std::unique_ptr<UniformLookupTable<IN_TYPE,OUT_TYPE>>(instance);
  }

  // Actual registration function
  static void RegisterFactoryFunction(std::string name,
      std::function<UniformLookupTable<IN_TYPE,OUT_TYPE>*(
        FunctionContainer<IN_TYPE,OUT_TYPE>*, 
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
  static std::map<std::string, std::function<UniformLookupTable<IN_TYPE,OUT_TYPE>*(
			FunctionContainer<IN_TYPE,OUT_TYPE>*, OTHER)>>& get_registry()
  {
    // Get the singleton instance of the registry map
    static std::map<std::string, std::function<UniformLookupTable<IN_TYPE,OUT_TYPE>*(
               FunctionContainer<IN_TYPE,OUT_TYPE>*,OTHER)>> registry;
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
template <class T, typename IN_TYPE, typename OUT_TYPE, class OTHER = UniformLookupTableParameters<IN_TYPE>>
class UniformLookupTableRegistrar {
public:
  UniformLookupTableRegistrar(std::string className)
  {
    UniformLookupTableFactory<IN_TYPE,OUT_TYPE>::RegisterFactoryFunction(className,
      [](FunctionContainer<IN_TYPE,OUT_TYPE> *fc,
         OTHER args
      ) -> UniformLookupTable<IN_TYPE,OUT_TYPE>* { return new T(fc, args); }
    );
  }
};

/*
   Macros for class registration:
   - FUNC_REGISTER_LUT goes inside class definition
   - FUNC_REGISTER_DOUBLE_AND_FLOAT_LUT_IMPLS goes underneath the class definition.
   Several different versions of this macro exist for registering templated classes
   - other... is for template parameters unrelated to the tables IN_TYPE and OUT_TYPE
*/
#define FUNC_REGISTER_LUT(classname) \
private: \
  static const UniformLookupTableRegistrar<classname,IN_TYPE,OUT_TYPE> registrar; \
  static const UniformLookupTableRegistrar<classname,IN_TYPE,OUT_TYPE,std::string> str_registrar

#define FUNC_STR_EXPAND(x...) #x
#define FUNC_STR(x...) FUNC_STR_EXPAND(x)

#define FUNC_REGISTER_LUT_IMPL(classname,IN_TYPE,OUT_TYPE) \
  template<> const \
  UniformLookupTableRegistrar<classname<IN_TYPE,OUT_TYPE>,IN_TYPE,OUT_TYPE> \
    classname<IN_TYPE,OUT_TYPE>::registrar(FUNC_STR(classname))

#define FUNC_REGISTER_TEMPLATED_LUT_IMPL(classname,IN_TYPE,OUT_TYPE,other...) \
  template<> const \
    UniformLookupTableRegistrar<classname<IN_TYPE,OUT_TYPE,other>,IN_TYPE,OUT_TYPE> \
    classname<IN_TYPE,OUT_TYPE,other>::registrar(FUNC_STR(classname<other>))

// macros used in each class to quickly register 4 different template instantiations
#define FUNC_REGISTER_DOUBLE_AND_FLOAT_LUT_IMPLS(classname)\
  FUNC_REGISTER_LUT_IMPL(classname,double,double); \
  FUNC_REGISTER_LUT_IMPL(classname,float,double);  \
  FUNC_REGISTER_LUT_IMPL(classname,double,float);  \
  FUNC_REGISTER_LUT_IMPL(classname,float,float);

#define FUNC_REGISTER_TEMPLATED_DOUBLE_AND_FLOAT_LUT_IMPLS(classname,other...) \
  FUNC_REGISTER_TEMPLATED_LUT_IMPL(classname,double,double,other); \
  FUNC_REGISTER_TEMPLATED_LUT_IMPL(classname,float,double,other);  \
  FUNC_REGISTER_TEMPLATED_LUT_IMPL(classname,double,float,other);  \
  FUNC_REGISTER_TEMPLATED_LUT_IMPL(classname,float,float,other)

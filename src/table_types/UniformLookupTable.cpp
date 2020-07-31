/* Implementation of a Uniform Lookup table */
#include "UniformLookupTable.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

template <typename IN_TYPE, typename OUT_TYPE>
UniformLookupTable::UniformLookupTable(FunctionContainer *func_container, UniformLookupTableParameters par) :
  EvaluationImplementation(func_container->standard_func, "uniform_lookup_table")
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

template <typename IN_TYPE, typename OUT_TYPE>
std::pair<IN_TYPE,IN_TYPE> UniformLookupTable::arg_bounds_of_interval(unsigned intervalNumber)
{
  return std::make_pair(m_grid[intervalNumber],m_grid[intervalNumber+1]);
}

template <typename IN_TYPE, typename OUT_TYPE>
void UniformLookupTable::print_details(std::ostream &out)
{
  out << m_name << " " << m_minArg << " " << m_maxArg << " "
      << m_stepSize << " " << m_numIntervals << " ";
}

/* ////////////////////////////////////////////////////////////////////////////
   Factory and registration management
//////////////////////////////////////////////////////////////////////////// */
template <typename IN_TYPE, typename OUT_TYPE>
std::unique_ptr<UniformLookupTable> UniformLookupTableFactory::Create(std::string name,
                              FunctionContainer *fc,
                              UniformLookupTableParameters par)
{
  // Create a UniformLookupTable
  UniformLookupTable * instance = nullptr;

  // find the name in the registry and call factory method.
  auto it = get_registry().find(name);
  if(it != get_registry().end())
    instance = it->second(fc,par);

  // wrap instance in a unique ptr and return (if created)
  if(instance == nullptr)
    throw std::invalid_argument(name + " not found in registry.");
  return std::unique_ptr<UniformLookupTable>(instance);
}

template <typename IN_TYPE, typename OUT_TYPE>
std::vector<std::string> UniformLookupTableFactory::get_registry_keys()
{
  // copy all keys from the registry map into a vector
  std::vector<std::string> keys;
  for (auto const& elem : get_registry() ) {
    keys.push_back(elem.first);
  }
  return keys;
}

template <typename IN_TYPE, typename OUT_TYPE>
void UniformLookupTableFactory::RegisterFactoryFunction(std::string name,
std::function<UniformLookupTable*(FunctionContainer*,UniformLookupTableParameters)> classFactoryFunction)
{
  // register a derived class factory function
  get_registry()[name] = classFactoryFunction;
}

template <typename IN_TYPE, typename OUT_TYPE>
std::map<std::string, std::function<UniformLookupTable*(
			       FunctionContainer*,UniformLookupTableParameters)>>& UniformLookupTableFactory::get_registry()
{
  // Get the singleton instance of the registry map
  static std::map<std::string, std::function<UniformLookupTable*(
			       FunctionContainer*,UniformLookupTableParameters)>> registry;
  return registry;
}

template class UniformLookupTable;

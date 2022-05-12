/* Implementation of a Uniform Lookup table */
#include "UniformLookupTable.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

UniformLookupTable::UniformLookupTable(FunctionContainer *func_container, UniformLookupTableParameters par) :
  EvaluationImplementation(func_container->double_func, "uniform_lookup_table")
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
  m_grid.reset(new double[m_numIntervals]);

}

std::pair<double,double> UniformLookupTable::arg_bounds_of_interval(unsigned intervalNumber)
{
  return std::make_pair(m_grid[intervalNumber],m_grid[intervalNumber+1]);
}

void UniformLookupTable::print_details(std::ostream &out)
{
  out << m_name << " " << m_minArg << " " << m_maxArg << " "
      << m_stepSize << " " << m_numIntervals << " ";
}

/* ////////////////////////////////////////////////////////////////////////////
   Factory and registration management
//////////////////////////////////////////////////////////////////////////// */
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

std::vector<std::string> UniformLookupTableFactory::get_registry_keys()
{
  // copy all keys from the registry map into a vector
  std::vector<std::string> keys;
  for (auto const& elem : get_registry() ) {
    keys.push_back(elem.first);
  }
  return keys;
}

void UniformLookupTableFactory::RegisterFactoryFunction(std::string name,
std::function<UniformLookupTable*(FunctionContainer*,UniformLookupTableParameters)> classFactoryFunction)
{
  // register a derived class factory function
  get_registry()[name] = classFactoryFunction;
}

std::map<std::string, std::function<UniformLookupTable*(
			       FunctionContainer*,UniformLookupTableParameters)>>& UniformLookupTableFactory::get_registry()
{
  // Get the singleton instance of the registry map
  static std::map<std::string, std::function<UniformLookupTable*(
			       FunctionContainer*,UniformLookupTableParameters)>> registry;
  return registry;
}

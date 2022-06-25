/*
  Factory for Lookup Tables
  - LookupTableFactory<TIN,TOUT>::create(str_name, fc, par) generates table types derived from LookupTable<TIN,TOUT>
  - Note: New implementations must be added to the registry by adding to the ::initialize() member function
*/
#pragma once
#include "TableIncludes.hpp"
#include "config.hpp" // FUNC_USE_BOOST, FUNC_USE_ARMADILLO
#include <memory>     // unique_ptr
#include <map>
#include <vector>
#include <string>

/*
 * Two-stage string expansion macros
 */
#define FUNC_STR_EXPAND(x...) #x
#define FUNC_STR(x...) FUNC_STR_EXPAND(x)

/*
 * Macros to easily add table types into the registry
 * - Call this inside ::initialize_registry() for all desired table types
 */
#define FUNC_REGISTER_ONE(classname)                                                                                   \
  registry.insert({FUNC_STR(classname), [](FunctionContainer<TIN, TOUT> *fc, OTHER args) -> LookupTable<TIN, TOUT> * { \
                     return new classname<TIN, TOUT>(fc, args);                                                        \
                   }})

#define FUNC_REGISTER_TEMPLATE(classname, templates...)                                                                \
  registry.insert(                                                                                                     \
       {FUNC_STR(classname<templates>), [](FunctionContainer<TIN, TOUT> *fc, OTHER args) -> LookupTable<TIN, TOUT> * { \
         return new classname<TIN, TOUT, templates>(fc, args);                                                         \
       }})

// This will have to change slightly depending on the maximum number of template arguments,
// but it's easy enough to understand, prevents name polution, and is needed in FunctionContainer
// anyways: https://stackoverflow.com/a/11763277/14090389
#define FUNC_REGISTRY_GET_MACRO(_1,_2,_3,NAME,...) NAME

// Call with FUNC_ADD_TABLE_TO_REGISTRY(classname,template_args_if_they_exist)
#define FUNC_ADD_TABLE_TO_REGISTRY(...) \
  FUNC_REGISTRY_GET_MACRO(__VA_ARGS__, FUNC_REGISTER_TEMPLATE, FUNC_REGISTER_TEMPLATE, FUNC_REGISTER_ONE, )(__VA_ARGS__)

// nonuniform LUTs require Boost
//#ifdef FUNC_USE_BOOST
//#define FUNC_ADD_TABLE_TO_REGISTRY(classname,...)    \
//  FUNC_ADD_TABLE_TO_REGISTRY(Uniform##classname,...);         \
//  FUNC_ADD_TABLE_TO_REGISTRY(NonUniform##classname,...);      \
//  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformPseudo##classname,...)
//#else
//#define FUNC_ADD_TABLE_TO_REGISTRY(classname,...)    \
//  FUNC_ADD_TABLE_TO_REGISTRY(Uniform##classname,...)
//#endif

template <typename TIN, typename TOUT = TIN, class OTHER = LookupTableParameters<TIN>> class LookupTableFactory {
public:

  /*
   * The map type that holds the registry
   */
  using registry_t =
      std::map<std::string, std::function<LookupTable<TIN, TOUT> *(FunctionContainer<TIN, TOUT> *, OTHER)>>;

  /*
   * Constructor initializes registry, default destructor.
   */
  LookupTableFactory() { initialize_registry(); };
  ~LookupTableFactory() = default;

  /*
   * Create a lookup table from
   * - string_name - Stringified table type
   * - fc          - FunctionContainer holding the function that the table evaluates
   * - args        - Additional arguments needed for construcing the table
   */
  std::unique_ptr<LookupTable<TIN, TOUT>> create(std::string string_name, FunctionContainer<TIN, TOUT> *fc, OTHER args);

  /*
   * Return a container of the keys for table types that have been registered
   */
  std::vector<std::string> get_registered_keys();

private:

  /*
   * Hold mapping from strings to constructors for derived table types.
   */
  registry_t registry;

  /*
   *  Register the desired table types (construct a valid map for the `registry`)
   */
  void initialize_registry();
};

/* --------------------------------------------------------------------------
 * --------------------------------------------------------------------------
 *      Implementation
 * --------------------------------------------------------------------------
 * -------------------------------------------------------------------------- */

/*
 *  Initialize the registry
 *  - New implementations of table types must be added to the registry here
 *  - TODO specialize OTHER=LUTparameters for a case where boost is not available (they're trying to generate a LUT in that case)
 */
template <typename TIN, typename TOUT, class OTHER> void LookupTableFactory<TIN, TOUT, OTHER>::initialize_registry() {

  // TODO our Taylor/Pade tables don't have nonuniform variants yet
  FUNC_ADD_TABLE_TO_REGISTRY(UniformConstantTaylorTable);
#ifdef FUNC_USE_BOOST
  FUNC_ADD_TABLE_TO_REGISTRY(UniformLinearTaylorTable);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformQuadraticTaylorTable);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformCubicTaylorTable);

#ifdef FUNC_USE_ARMADILLO
  // Pade tables need both Boost and Armadillo to build
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,1,1);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,2,1);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,3,1);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,4,1);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,5,1);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,6,1);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,2,2);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,3,2);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,4,2);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,5,2);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,3,3);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,4,3);
#endif
#endif

  FUNC_ADD_TABLE_TO_REGISTRY(UniformCubicHermiteTable);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformCubicPrecomputedInterpolationTable);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformQuadraticPrecomputedInterpolationTable);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformLinearPrecomputedInterpolationTable);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformLinearInterpolationTable);

#ifdef FUNC_USE_ARMADILLO
  FUNC_ADD_TABLE_TO_REGISTRY(UniformArmadilloPrecomputedInterpolationTable,4);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformArmadilloPrecomputedInterpolationTable,5);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformArmadilloPrecomputedInterpolationTable,6);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformArmadilloPrecomputedInterpolationTable,7);
#endif
}


template <typename TIN, typename TOUT> class LookupTableFactory<TIN,TOUT,std::string>{
  std::map<std::string, std::function<LookupTable<TIN, TOUT> *(FunctionContainer<TIN, TOUT> *, std::string)>> registry;
  void initialize_registry(); // only this function changes
public:  
  std::unique_ptr<LookupTable<TIN, TOUT>> create(std::string string_name, FunctionContainer<TIN, TOUT> *fc, std::string args);
  std::vector<std::string> get_registered_keys(); // does this actually need to change?
};

/* We always want to be able to read these types from a json file */
template <typename TIN, typename TOUT> void LookupTableFactory<TIN, TOUT, std::string>::initialize_registry() {
#define OTHER std::string // this is a big kludge
  // TODO our Taylor/Pade tables don't have nonuniform variants yet
  FUNC_ADD_TABLE_TO_REGISTRY(UniformCubicTaylorTable);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformQuadraticTaylorTable);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformLinearTaylorTable);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformConstantTaylorTable);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,1,1);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,2,1);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,3,1);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,4,1);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,5,1);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,6,1);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,2,2);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,3,2);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,4,2);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,5,2);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,3,3);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformPadeTable,4,3);

  FUNC_ADD_TABLE_TO_REGISTRY(UniformCubicHermiteTable);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformCubicPrecomputedInterpolationTable);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformQuadraticPrecomputedInterpolationTable);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformLinearPrecomputedInterpolationTable);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformLinearInterpolationTable);

  FUNC_ADD_TABLE_TO_REGISTRY(UniformArmadilloPrecomputedInterpolationTable,4);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformArmadilloPrecomputedInterpolationTable,5);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformArmadilloPrecomputedInterpolationTable,6);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformArmadilloPrecomputedInterpolationTable,7);
#undef OTHER
}

/*
 *  Return a vector of the keys that have been registered
 */
template <typename TIN, typename TOUT, class OTHER>
std::vector<std::string> LookupTableFactory<TIN, TOUT, OTHER>::get_registered_keys() {
  // copy all keys from the registry map into a vector
  std::vector<std::string> keys;
  for (auto const &elem : registry)
    keys.push_back(elem.first);
  return keys;
}

/*
 *  Create a new lookup table. Throw exception asking for an unregistered table.
 */
template <typename TIN, typename TOUT, class OTHER>
std::unique_ptr<LookupTable<TIN, TOUT>>
LookupTableFactory<TIN, TOUT, OTHER>::create(std::string name, FunctionContainer<TIN, TOUT> *fc, OTHER args) {
  // Create a LookupTable
  LookupTable<TIN, TOUT> *instance = nullptr;

  // find the name in the registry and call factory method.
  auto it = registry.find(name);
  if (it != registry.end())
    instance = it->second(fc, args);

  // wrap instance in a unique ptr and return (if created)
  // TODO typos are super common so suggest a "close match" to name in the thrown error
  if (instance == nullptr)
    throw std::invalid_argument(name + " not found in registry.");
  return std::unique_ptr<LookupTable<TIN, TOUT>>(instance);
}

// Legacy func typedef (TODO is this ever used anymore???)
//template <typename TIN, typename TOUT = TIN, class OTHER = LookupTableParameters<TIN>>
//using UniformLookupTableFactory = LookupTableFactory<TIN,TOUT,OTHER>;


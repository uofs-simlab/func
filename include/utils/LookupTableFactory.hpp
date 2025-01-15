#pragma once
#include "config.hpp" // FUNC_USE_BOOST, FUNC_USE_ARMADILLO, FUNC_DECLARE_TEMPLATE_AS_EXTERN
#include "tables.hpp"
#include <memory>     // unique_ptr
#include <map>
#include <vector>
#include <string>
#include <stdexcept>

/*
 * Two-stage string expansion macros
 */
#define FUNC_STR_EXPAND(...) #__VA_ARGS__
#define FUNC_STR(...) FUNC_STR_EXPAND(__VA_ARGS__)

/*
 * Macros to easily add table types into the registry
 * - Call this inside ::initialize_registry() for all desired table types
 */
#define FUNC_REGISTER_ONE(classname)                                                                                                                    \
  registry.insert({FUNC_STR(classname), [](const FunctionContainer<TIN, TOUT>& fc, const LookupTableParameters<TIN>& args, const nlohmann::json& jsonStats) -> LookupTable<TIN, TOUT> * { \
                     return new classname<TIN, TOUT>(fc, args, jsonStats);                                                                              \
                   }})

#define FUNC_REGISTER_TEMPLATE(classname, ...)                                                                                                 \
  registry.insert(                                                                                                                                      \
       {FUNC_STR(classname<__VA_ARGS__>), [](const FunctionContainer<TIN, TOUT>& fc, const LookupTableParameters<TIN>& args, const nlohmann::json& jsonStats) -> LookupTable<TIN, TOUT> * { \
         return new classname<__VA_ARGS__, TIN, TOUT>(fc, args, jsonStats);                                                                               \
       }})

// This will have to change slightly depending on the maximum number of template arguments,
// but it's easy enough to understand, prevents name polution, and is needed in FunctionContainer
// anyways: https://stackoverflow.com/a/11763277/14090389
#define FUNC_REGISTRY_GET_MACRO(_1,_2,_3,NAME,...) NAME

// Call with FUNC_ADD_TABLE_TO_REGISTRY(classname,template_args_if_they_exist)
#define FUNC_ADD_TABLE_TO_REGISTRY(...) \
  FUNC_REGISTRY_GET_MACRO(__VA_ARGS__, FUNC_REGISTER_TEMPLATE, FUNC_REGISTER_TEMPLATE, FUNC_REGISTER_ONE, )(__VA_ARGS__)

namespace func {

/**
  Factory for Lookup Tables
  - LookupTableFactory<TIN,TOUT>::create(str_name, fc, par) generates table types derived from LookupTable<TIN,TOUT>
  - Note: New implementations must be added to the registry by adding to the ::initialize() member function
*/
template <typename TIN, typename TOUT = TIN> class LookupTableFactory {
public:

  /* The map type that holds the registry */
  using registry_t =
      std::map<std::string, std::function<LookupTable<TIN, TOUT> *(const FunctionContainer<TIN, TOUT>&, const LookupTableParameters<TIN>&, const nlohmann::json& jsonStats)>>;

  /* Constructor initializes registry, default destructor. */
  LookupTableFactory() { initialize_registry(); };
  ~LookupTableFactory() = default;

  /*
   * Create a lookup table from
   * - string_name - Stringified table type
   * - fc          - FunctionContainer holding the function that the table evaluates
   * - args        - Additional arguments needed for construcing the table
   */
  std::unique_ptr<LookupTable<TIN, TOUT>> create(std::string string_name, const FunctionContainer<TIN, TOUT>& fc, const LookupTableParameters<TIN>& args,
      const nlohmann::json& jsonStats=nlohmann::json());

  /* Return a container of the keys for table types that have been registered */
  std::vector<std::string> get_registered_keys();

private:

  /* Hold mapping from strings to constructors for derived table types. */
  registry_t registry;

  /* Register every LookupTable implementation we officially support (construct the registry) */
  void initialize_registry();
};


/* --------------------------------------------------------------------------
 * --------------------------------------------------------------------------
 *      Implementation
 * --------------------------------------------------------------------------
 * -------------------------------------------------------------------------- */

/**  Initialize the registry
 *  - New implementations of table types must be added to the registry here
 */
template <typename TIN, typename TOUT>
void LookupTableFactory<TIN, TOUT>::initialize_registry() {
  FUNC_ADD_TABLE_TO_REGISTRY(UniformTaylorTable,1);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformTaylorTable,2);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformTaylorTable,3);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformTaylorTable,4);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformTaylorTable,5);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformTaylorTable,6);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformTaylorTable,7);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformCubicHermiteTable);

  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformTaylorTable,1);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformTaylorTable,2);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformTaylorTable,3);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformTaylorTable,4);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformTaylorTable,5);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformTaylorTable,6);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformTaylorTable,7);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformCubicHermiteTable);

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

  FUNC_ADD_TABLE_TO_REGISTRY(UniformChebyInterpTable,1);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformChebyInterpTable,2);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformChebyInterpTable,3);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformChebyInterpTable,4);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformChebyInterpTable,5);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformChebyInterpTable,6);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformChebyInterpTable,7);

  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformChebyInterpTable,1);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformChebyInterpTable,2);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformChebyInterpTable,3);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformChebyInterpTable,4);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformChebyInterpTable,5);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformChebyInterpTable,6);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformChebyInterpTable,7);

  FUNC_ADD_TABLE_TO_REGISTRY(UniformExactInterpTable,0);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformExactInterpTable,1);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformExactInterpTable,2);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformExactInterpTable,3);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformExactInterpTable,4);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformExactInterpTable,5);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformExactInterpTable,6);

  FUNC_ADD_TABLE_TO_REGISTRY(UniformLinearRawInterpTable);

  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformExactInterpTable,1);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformExactInterpTable,2);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformExactInterpTable,3);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformExactInterpTable,4);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformExactInterpTable,5);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformExactInterpTable,6);

}

/**
 *  Return a vector of the keys that have been registered
 */
template <typename TIN, typename TOUT>
std::vector<std::string> LookupTableFactory<TIN, TOUT>::get_registered_keys() {
  // copy all keys from the registry map into a vector
  std::vector<std::string> keys;
  for (auto const &elem : registry)
    keys.push_back(elem.first);
  return keys;
}

/**
 *  Create a new lookup table. Throw exception asking for an unregistered table.
 */
template <typename TIN, typename TOUT>
std::unique_ptr<LookupTable<TIN, TOUT>>
LookupTableFactory<TIN, TOUT>::create(std::string name, const FunctionContainer<TIN, TOUT>& fc, const LookupTableParameters<TIN>& args, const nlohmann::json& jsonStats) {
  // Create a LookupTable
  LookupTable<TIN, TOUT> *instance = nullptr;

  // find the name in the registry and call factory method.
  auto it = registry.find(name);
  if (it != registry.end())
    instance = it->second(fc, args, jsonStats); // we found the constructor corresponding to name

  // wrap instance in a unique ptr and return (if created)
  // TODO typos are super common so suggest a "close match" to name in the thrown error
  if (instance == nullptr)
    throw std::invalid_argument(name + " not found in registry.");
  return std::unique_ptr<LookupTable<TIN, TOUT>>(instance);
}


/** from_json for unique_ptr<LookupTable> unlocks this fancy syntax:
 *
 * ```cpp
 *   nlohmann::json jsonStats;
 *   std::ifstream(filename) >> jsonStats;
 *   auto lut = jsonStats.get<std::unique_ptr<func::LookupTable<TIN,TOUT>>>(); // call the constructor (or the from_json) referred to by "name"
 * ```
 * This is only possible because std::unique_ptr<T> is default constructable
*/
template <typename TIN, typename TOUT>
void from_json(const nlohmann::json& jsonStats, std::unique_ptr<LookupTable<TIN,TOUT>>& lut) {
  LookupTableFactory<TIN,TOUT> factory;
  lut = factory.create(jsonStats);
}

FUNC_DECLARE_TEMPLATE_AS_EXTERN(LookupTableFactory)

} // namespace func

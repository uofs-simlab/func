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
#include <stdexcept>

/*
 * Two-stage string expansion macros
 */
#define FUNC_STR_EXPAND(x...) #x
#define FUNC_STR(x...) FUNC_STR_EXPAND(x)

/*
 * Macros to easily add table types into the registry
 * - Call this inside ::initialize_registry() for all desired table types
 */
#define FUNC_REGISTER_ONE(classname)                                                                                                                    \
  registry.insert({FUNC_STR(classname), [](FunctionContainer<TIN, TOUT> *fc, LookupTableParameters<TIN> args, const nlohmann::json& jsonStats) -> LookupTable<TIN, TOUT> * { \
                     return new classname<TIN, TOUT>(fc, args, jsonStats);                                                                              \
                   }})

#define FUNC_REGISTER_TEMPLATE(classname, templates...)                                                                                                 \
  registry.insert(                                                                                                                                      \
       {FUNC_STR(classname<templates>), [](FunctionContainer<TIN, TOUT> *fc, LookupTableParameters<TIN> args, const nlohmann::json& jsonStats) -> LookupTable<TIN, TOUT> * { \
         return new classname<TIN, TOUT, templates>(fc, args, jsonStats);                                                                               \
       }})

// This will have to change slightly depending on the maximum number of template arguments,
// but it's easy enough to understand, prevents name polution, and is needed in FunctionContainer
// anyways: https://stackoverflow.com/a/11763277/14090389
#define FUNC_REGISTRY_GET_MACRO(_1,_2,_3,NAME,...) NAME

// Call with FUNC_ADD_TABLE_TO_REGISTRY(classname,template_args_if_they_exist)
#define FUNC_ADD_TABLE_TO_REGISTRY(...) \
  FUNC_REGISTRY_GET_MACRO(__VA_ARGS__, FUNC_REGISTER_TEMPLATE, FUNC_REGISTER_TEMPLATE, FUNC_REGISTER_ONE, )(__VA_ARGS__)

namespace func {

template <typename TIN, typename TOUT = TIN> class LookupTableFactory {
public:

  /*
   * The map type that holds the registry
   */
  using registry_t =
      std::map<std::string, std::function<LookupTable<TIN, TOUT> *(FunctionContainer<TIN, TOUT> *, LookupTableParameters<TIN>, const nlohmann::json& jsonStats)>>;

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
  std::unique_ptr<LookupTable<TIN, TOUT>> create(std::string string_name, FunctionContainer<TIN, TOUT> *fc, LookupTableParameters<TIN> args,
      const nlohmann::json& jsonStats=nlohmann::json());

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
   *  - templated because we have to specialize for std::string
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
 */
template <typename TIN, typename TOUT>
void LookupTableFactory<TIN, TOUT>::initialize_registry() {
  // TODO our Taylor/Pade tables don't have nonuniform variants yet
  FUNC_ADD_TABLE_TO_REGISTRY(UniformConstantTaylorTable);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformLinearTaylorTable);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformQuadraticTaylorTable);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformCubicTaylorTable);

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
  FUNC_ADD_TABLE_TO_REGISTRY(UniformCubicInterpolationTable);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformQuadraticInterpolationTable);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformLinearInterpolationTable);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformLinearRawInterpolationTable);

  FUNC_ADD_TABLE_TO_REGISTRY(UniformArmadilloInterpolationTable,4);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformArmadilloInterpolationTable,5);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformArmadilloInterpolationTable,6);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformArmadilloInterpolationTable,7);

  FUNC_ADD_TABLE_TO_REGISTRY(UniformChebyInterpolationTable,1);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformChebyInterpolationTable,2);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformChebyInterpolationTable,3);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformChebyInterpolationTable,4);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformChebyInterpolationTable,5);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformChebyInterpolationTable,6);
  FUNC_ADD_TABLE_TO_REGISTRY(UniformChebyInterpolationTable,7);


  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformCubicHermiteTable);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformCubicInterpolationTable);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformQuadraticInterpolationTable);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformLinearInterpolationTable);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformLinearRawInterpolationTable);

  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformArmadilloInterpolationTable,4);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformArmadilloInterpolationTable,5);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformArmadilloInterpolationTable,6);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformArmadilloInterpolationTable,7);

  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformPseudoCubicHermiteTable);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformPseudoCubicInterpolationTable);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformPseudoQuadraticInterpolationTable);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformPseudoLinearInterpolationTable);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformPseudoLinearRawInterpolationTable);

  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformPseudoArmadilloInterpolationTable,4);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformPseudoArmadilloInterpolationTable,5);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformPseudoArmadilloInterpolationTable,6);
  FUNC_ADD_TABLE_TO_REGISTRY(NonUniformPseudoArmadilloInterpolationTable,7);
}

/*
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

/*
 *  Create a new lookup table. Throw exception asking for an unregistered table.
 */
template <typename TIN, typename TOUT>
std::unique_ptr<LookupTable<TIN, TOUT>>
LookupTableFactory<TIN, TOUT>::create(std::string name, FunctionContainer<TIN, TOUT> *fc, LookupTableParameters<TIN> args, const nlohmann::json& jsonStats) {
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


/* from_json for unique_ptr<LookupTable> unlocks this fancy syntax:
```c++
  nlohmann::json jsonStats;
  std::ifstream(filename) >> jsonStats;
  auto lut = jsonStats.get<std::unique_ptr<func::LookupTable<TIN,TOUT>>>(); // call the constructor (or the from_json) referred to by "name"
```
Only drawback of this function is we lose access to the function that generated that LUT
this works because std::unique_ptr<T> is default constructable
*/
template <typename TIN, typename TOUT>
void from_json(const nlohmann::json& jsonStats, std::unique_ptr<LookupTable<TIN,TOUT>>& lut) {
  std::string name = jsonStats.at("name");
  LookupTableFactory<TIN,TOUT> factory;

  lut = factory.create(name, nullptr, LookupTableParameters<TIN>{0,0,0}, jsonStats);
}

FUNC_DECLARE_TEMPLATE_AS_EXTERN(LookupTableFactory)

} // namespace func

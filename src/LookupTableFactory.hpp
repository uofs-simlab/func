/*
  Factory for Lookup Tables
  - LookupTableFactory<TIN,TOUT>::create(str_name, fc, par) generates table types derived from LookupTable<TIN,TOUT>
  - Note: New implementations must be added to the registry by adding to the ::initialize() member function
*/
#pragma once
#include "TableIncludes.hpp"
#include "config.hpp" // FUNC_USE_QUADMATH
#include <memory>     // unique_ptr

#define FUNC_STR_EXPAND(x...) #x
#define FUNC_STR(x...) FUNC_STR_EXPAND(x)

/*
 * Macro to easily add table types into the registry
 * - Call this inside ::initialize_registry() for all desired table types
 */
#define FUNC_ADD_TO_REGISTRY(X)                                                                                        \
  registry.insert({FUNC_STR(X), [](FunctionContainer<TIN, TOUT> *fc, OTHER args) -> LookupTable<TIN, TOUT> * {         \
                     return new X<TIN, TOUT>(fc, args);                                                                \
                   }});

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
 *      Implementations
 * --------------------------------------------------------------------------
 * -------------------------------------------------------------------------- */

/*
 *  Initialize the registry
 *  - New implementations of table types must be added to the registry here
 */
template <typename TIN, typename TOUT, class OTHER> void LookupTableFactory<TIN, TOUT, OTHER>::initialize_registry() {

  FUNC_ADD_TO_REGISTRY(UniformCubicHermiteTable);
  FUNC_ADD_TO_REGISTRY(UniformCubicPrecomputedInterpolationTable);
  FUNC_ADD_TO_REGISTRY(UniformCubicTaylorTable);
  FUNC_ADD_TO_REGISTRY(UniformLinearInterpolationTable);
  FUNC_ADD_TO_REGISTRY(UniformLinearPrecomputedInterpolationTable);
  FUNC_ADD_TO_REGISTRY(UniformLinearTaylorTable);
  FUNC_ADD_TO_REGISTRY(UniformQuadraticPrecomputedInterpolationTable);
  FUNC_ADD_TO_REGISTRY(UniformQuadraticTaylorTable);

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
  if (instance == nullptr)
    throw std::invalid_argument(name + " not found in registry.");
  return std::unique_ptr<LookupTable<TIN, TOUT>>(instance);
}

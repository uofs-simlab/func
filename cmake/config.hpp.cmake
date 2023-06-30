#pragma once
// Let the hpp files know which dependencies CMake found
#cmakedefine FUNC_USE_BOOST
#cmakedefine FUNC_USE_ARMADILLO

#if (__cplusplus > 201700L || _MSVC_LANG > 201700L)
#define FUNC_IF_CONSTEXPR if constexpr
#else
#define FUNC_IF_CONSTEXPR if
#define FUNC_NO_CXX17_IF_CONSTEXPR
#endif

// This is a cmake variable defined in CMakeLists.txt
// and will define the macro FUNC_DECLARE_TEMPLATE_AS_EXTERN(classname)
@FUNC_DECLARE_TEMPLATE_AS_EXTERN@


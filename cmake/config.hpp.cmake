#pragma once
// Let our hpp files know which dependencies we could find
#cmakedefine FUNC_USE_BOOST
#cmakedefine FUNC_USE_ARMADILLO

// This is a cmake variable defined in CMakeLists.txt
// and will define the macro FUNC_DECLARE_TEMPLATE_AS_EXTERN(classname)
@FUNC_DECLARE_TEMPLATE_AS_EXTERN@


#pragma once
// Let our hpp files know which dependencies we could find
#define FUNC_USE_BOOST
#define FUNC_USE_ARMADILLO

// This is a cmake variable defined in CMakeLists.txt
// and will define the macro FUNC_DECLARE_TEMPLATE_AS_EXTERN(classname)
#define FUNC_DECLARE_TEMPLATE_AS_EXTERN(classname)\
  extern template class classname<double,double>;\


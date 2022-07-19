/* CMake will automatically fill in CLASSNAME, and PAIR (typename pair)
   for each class that we want cpp files for */
#include "@CLASSNAME@.hpp"

/* Remember to declare extern template ... in the corresponding hpp file
   otherwise this file is useless */
template class func::@CLASSNAME@<@PAIR@@TEMPLATES@>;

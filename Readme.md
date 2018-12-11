FunC
====

`FunC` (Function Comparator) is a C++ tool for quickly profiling the performance of various different abstracted implementations of mathematical function evaluations for given:
- computing environment
- mathematical function
- evaluation tolerance

Currently, only uniform LUTs (interpolation, Taylor, and Hermite) and direct evaluations are supported.


Requirements
------------

- C++11 compliant compiler
- Boost (tested with 1.58-1.65)
- quadmath (for tolerance-based table generation capabilities)

Build:

- CMake version >= 3.1

    mkdir build && cd build/
	cmake -DCMAKE_INSTALL_PREFIX=<install-dir>  ..
	make install


Notes
-----

This tool is split up into two separate libraries:
- `func_impls` provides the implementation types used in the library (does not require quadmath to build/use)
- `func` combines the implementation types with
  - `UniformLookupTableGenerator` for generating LUTs according to various criteria (requires quadmath)
  - `ImplementationComparator` for comparing the performance of various implementations


References
----------

TODO reference to SISC paper


Copyright
---------

TODO decide on copyright

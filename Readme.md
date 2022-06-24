FunC
====

`FunC` (Function Comparator) is a C++ tool for approximating any single variable, real valued function $f:\mathbb{R}\to\mathbb{R}$ with a lookup table (LUT). This tool aims to steamline the process of finding an optimal LUT of $f$ for a given program. This includes factors such as
- error tolerance
- domain usage (ie the inputs to $f$ during runtime)
- evaluation frequency (ie how much work is being done in between calls to $f$)

`FunC` can build LUTs using interpolation (up to degree 7), Taylor polynomials, Pade approximants, or Hermite interpolation over a uniform or automatically generated nonuniform mesh.


Requirements
------------

- C++11 compliant compiler (tested with gcc and clang)
- Boost 1.71.0 or newer*
- Armadillo (tested with 9.1-10.1)*

\*These libraries are only required for _table generation_. They are _not_ required if every table is being read from a json file (which is an optimization that most production level code will make).

### Build:

- CMake version >= 3.1
```
mkdir build && cd build/
cmake -DCMAKE_INSTALL_PREFIX=<install-dir>  ..
make install
```

After make install, linking to the library (outside of cmake build) requires:
- `<install-dir>/lib` is in your `LD_LIBRARY_PATH` environment variable,
- `<install-dir>/include/func` is in your include flags, and
- `-lfunc` are in your link flags

Notes
-----
- `FunC` is "pseudo header only." All of its functionality is available through headers only. Link with `libfunc.so` to greatly speed up the compile time of any program using our `LookupTableGenerator`.


References
----------

TODO reference to SISC paper


Copyright
---------

TODO decide on copyright

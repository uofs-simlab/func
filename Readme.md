FunC
====

`FunC` (Function Comparator) is a C++ tool for approximating any single-variable, pure function $f:\texttt{TIN}\to \texttt{TOUT}$ with a lookup table (LUT) over a closed interval $[a,b]$. $\texttt{TIN}$ and $\texttt{TOUT}$ are types with overloads for $\texttt{operator+,-}$, and there is a commutative $\texttt{operator*:TIN}\times\texttt{TOUT}\to\texttt{TOUT}$ (so $\texttt{TOUT}$ forms a vector space over the field $\texttt{TIN}$). We take a LUT as any piecewise approximation of $f$, so a LUT of $f$ takes the following form.

```math
\begin{equation}
L(x) = \begin{cases} p_0(x), & \text{ if } x_0 \leq x < x_1, \\
p_1(x), & \text{ if } x_1 \leq x < x_2,
\\ \quad \vdots \\
p_{N-1}(x), & \text{ if } x_{N-1} \leq x \leq x_N,\end{cases}
\end{equation}
```

`FunC` can build LUTs where each $p_k$ in the equation above are interpolating polynomials (up to degree 7 with Chebyshev nodes or degree 3 with equally spaced nodes), Taylor polynomials (up to degree 7), Pade approximants, or degree 3 Hermite interpolating polynomials. The $x_k$ in the equation above partition $[a,b]$. They can form a uniform partition (so $x_k=a+k(b-a)/N$) or be an automatically generated nonuniform partition.

`FunC` aims to streamline finding a good LUT of $f$ for a user's application. To do so, we must measure factors such as
- absolute and relative tolerances for error
- domain usage (i.e. the inputs to $f$ during the user's program's runtime)
- evaluation frequency (i.e. how much work is being done in between calls to $f$)
`FunC`'s DirectEvaluation class measures the first two, and a LookupTableGenerator optimizes a LUT's step size according to maximum tolerances for error.


Requirements
------------

- C++14 compliant compiler (tested with g++ and clang++)
- Boost 1.71.0 or newer*
- Armadillo (tested with 9.1-12.6)*

\*Boost and Armadillo are only required for _generating LUTs_. They are _not_ required if every LUT used is read from a JSON file (which is an optimization that most production-level code will make).

### Build:

- CMake version >= 3.13
```
mkdir build && cd build/
cmake -DCMAKE_INSTALL_PREFIX=<install-dir> ..
make install
```

After make install, linking to the library (outside of cmake build) requires:
- `<install-dir>/lib` is in your `LD_LIBRARY_PATH` environment variable,
- `<install-dir>/include/func` is in your include flags, and
- `-lfunc` is used when linking

Example usage
-------------
The following code shows an MWE of building a LUT of $f$.
```

#include <boost/math/special_functions/jacobi_elliptic.hpp>
template <typename T> T f(T x){ return boost::math::jacobi_cn(0.5, x); } // user's function

/* If FUNC_DEBUG is defined before including func.hpp, then any DirectEvaluation or FailureProofTable will have a histogram that stores each argument passed to their operator() during program runtime */
// #define FUNC_DEBUG
#include <func/func.hpp>
#include <iostream>
int main(){
    /* FUNC_SET_F is a macro */
    func::FunctionContainer<double> fc {FUNC_SET_F(f,double)};
    func::DirectEvaluation<double> de {fc};
    /* build a LUT of f over [0,2] with a step size of h=0.1. Each subinterval will use degree 3 Chebyshev interpolating polynomials */
    func::UniformChebyInterpTable<3,double> lut {fc, {0.0,2.0,0.1}};

    std::cout << de(1.011) << "\n"; // Call the function on line 7 with T=double and x=1.011
    std::cout << lut(1.011) << "\n"; // return an approximation of f(1.011) with a piecewise cubic polynomial
    std::cout << lut << "\n"; // print info about lut
    std::cout << de << std::endl; // print histogram to stdout if FUNC_DEBUG is defined
}
```


Notes
-----
- `FunC` can be used as if it is a header-only library (so templates TIN/TOUT can be set to ); however, files that take a long time to compile (eg. LookupTableFactory.hpp) are compiled for various user-provided types (float, double, long double). As such, linking with `libfunc.so` is not necessary but will greatly speed up the compile time of user code.


References
----------

TODO reference to SISC paper


Copyright
---------

TODO decide on copyright

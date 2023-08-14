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

\*Boost and Armadillo are only required for _generating LUTs_. They are _not_ required if every LUT used is read from a JSON file (which is a good optimization to make if $f$ is particularly slow to evaluate).

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
The following code shows an MWE of building a DirectEvaluation and a LUT of an elliptic function.
A DirectEvaluation can record every argument passed to its operator() in a histogram which is useful for determining useful bounds $a,b$ for a LUT. Furthermore, a DirectEvaluation can perturb its arguments 

```c++
/* user's function here. To use automatic differentiation, then the definition of f must be templated, and any function f
 * calls must have overloads for Boost's autodiff_fvar */
#include <boost/math/special_functions/jacobi_elliptic.hpp>
template <typename T> T f(T x){ return boost::math::jacobi_cn(0.5, x); }

/* If FUNC_DEBUG is defined before including func.hpp, then any DirectEvaluation or FailureProofTable will have a histogram
 * that stores each argument passed to their operator() during program runtime */
// #define FUNC_DEBUG
#include <func/func.hpp>
#include <iostream>
int main(){
    /* FUNC_SET_F is a macro required to take advantage of Boost's automatic differentiation.
     * - If f is templated on two types, call as FUNC_SET_F(f,TIN,TOUT)
     * - If f cannot be templated as shown  FunctionContainer could be constructed with f<double>, but */
    func::FunctionContainer<double> fc {FUNC_SET_F(f,double)};

    /* Arguments to a DirectEvaluation are (FunctionContainer fc, TIN min=0, TIN max=1, uint nbins=10, TOUT aerr=0, TIN rerr=0)
     * where min,max are used as bounds for the histogram */
    func::DirectEvaluation<double> de {fc,0,2,10,1,1};
    /* Call the function on line 7 with T=double and x=1.011. If FUNC_DEBUG is defined then f(x)(1+R*rerr)+A*aerr is returned
     * instead where A,R are random numbers sampled from a uniformly distributed random variable over [0,1] */
    std::cout << de(1.011) << "\n"; 
    std::cout << de << std::endl; // print histogram to stdout if FUNC_DEBUG is defined

    /* build a LUT of f over [0,2] with a step size of h=0.1. Each subinterval will use degree 3 Chebyshev interpolating polynomials */
    func::UniformChebyInterpTable<3,double> lut {fc, {0.0,2.0,0.1}};
    std::cout << lut(1.011) << "\n"; // return an approximation of f(1.011) with a piecewise cubic polynomial

    std::cout << lut << "\n"; // print info about lut
}
```


Notes
-----
- `FunC` can be used as a header-only library. So, TIN/TOUT can be any types such that TIN approximates a field and TOUT forms a vector space over TIN. Files that take a long time to compile (eg. LookupTableFactory.hpp, LookupTableComparator, etc) are compiled for various common types (float, double, long double). As such, linking with `libfunc.so` is unnecessary but will greatly speed up the compile time of user-level code.


References
----------

TODO reference to SISC paper


Copyright
---------

TODO decide on copyright

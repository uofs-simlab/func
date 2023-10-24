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

After `make install`, linking to the library (outside of cmake build) requires:
- `<install-dir>/lib` is in your `LD_LIBRARY_PATH` environment variable,
- `<install-dir>/include/func` is in your include flags, and
- `-lfunc` is used when linking


FunC classes and examples
-------------
The following code shows an MWE building various FunC objects out of an elliptic function. Each of FunC's LUTs are threadsafe.
A DirectEvaluation's operator() directly calls the function it was constructed with. It can
1. Record every argument passed to its operator() in a histogram. This is useful for determining useful bounds $a,b$ for a LUT of f.
2. Randomly perturb its outputs to simulate error.

Each of FunC's LUTs uses some piecewise function to approximate the user's function. LUTs are all constructed with a
FunctionContainer<TIN,TOUT> and a LookupTableParameters<TIN>. See the example below.

By default, LUTs do not perform bounds checking. To perform bounds checking, call a 
```c++
/* user's function here. To use automatic differentiation, the definition of f must be templated, and any function that f
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
     * - If f cannot be templated as shown  FunctionContainer could be constructed with f<double>, but then any LUT that requires derivative information cannot be built */
    func::FunctionContainer<double> fc {FUNC_SET_F(f,double)};

    /* DirectEvaluation's constructor takes arguments (FunctionContainer fc, TIN min=0, TIN max=1, uint nbins=10, TOUT aerr=0, TIN rerr=0)
     * - min, max, nbins, aerr, rerr are only used if the preprocessor macro FUNC_DEBUG is defined. 
     * - if FUNC_DEBUG is defined then min and max are used as bounds for the histogram. nbins is the number of buckets used
     * - if FUNC_DEBUG is defined then the DirectEvaluation returns f(x)(1+R*rerr)+A*aerr instead of f(x) where R,A are
     *   sampled from a uniformly distributed random variable over [0,1] */
    func::DirectEvaluation<double> de {fc,0,2,10,1,1};

    /* Call the function on line 7 with T=double and x=1.011.  */
    std::cout << de(1.011) << "\n"; 
    std::cout << de << std::endl; // print histogram of arguments to stdout if FUNC_DEBUG is defined

    /* build a LUT of f over [0,2] with a step size of h=0.1. Each subinterval will use degree 3 Chebyshev interpolating polynomials */
    func::UniformChebyInterpTable<3,double> cheb3 {fc, {0.0, 2.0, 0.1}};
    std::cout << cheb3(1.011) << "\n"; // return an approximation of f(1.011) with a piecewise cubic polynomial
    std::cout << cheb3 << "\n"; // print info about cheb3

    /* FunC can also build fast LUTs over nonuniform partitions of [a,b] with the following. These LUTs perform interval search
     * in 6 FLOPs and 4 IOPs which is far better than an O(logn) binary search. Due to a quirk of this optimization, FunC's
     * nonuniform LUTs have lower absolute error only if |f'(a)| or |f'(b)| are large. Nonuniform LUTs use the same interface as
     * uniform LUTs because the user has no control over the constructed partition. */
    func::NonUniformChebyInterpTable<3,double> nonunif {fc, {0.0, 2.0, 0.1}};
    std::cout << nonunif(1.011) << "\n"; // return an approximation of f(1.011) with a piecewise cubic polynomial
    std::cout << nonunif << "\n"; // print info about nonunif

    /* By default, a LUT does not perform bounds checking, so calling cheb3(2.5) is undefined (likely a segfault).
     * If the user wishes to perform bounds checking every time they call operator() then they can build a
     * FailureProofTable. The constructor takes arguments (FunctionContainer fc, LookupTableParameters<TIN>, min=1, max=-1, nbins=10). */
    func::FailureProofTable<func::UniformChebyInterpTable<3,double>,double> F(fc, {0.0, 2.0, 0.1});
    std::cout << F(2.5) << "\n";

    /* the following line prints a histogram of out-of-bounds arguments to stdout if FUNC_DEBUG is defined. If max<min then just the largest
     * and smallest arguments are displayed */
    std::cout << F << "\n";
}
```
Running `examples/experiment.cpp` will show each of FunC's currently supported LUT types. Users can easily benchmark each LUT with
their own required bounds & step size by modifying examples/experiment.cpp to use their own function.


Notes
-----
- `FunC` can be used as a header-only library. So, TIN/TOUT can be any types such that TIN approximates a field and TOUT forms a vector space over TIN. Files that take a long time to compile (eg. LookupTableFactory.hpp, LookupTableComparator, etc) are compiled for various common types (float, double, long double). As such, linking with `libfunc.so` is unnecessary but will greatly speed up the compile time of user-level code.


References
----------

TODO reference to SISC paper


Copyright
---------

TODO decide on copyright

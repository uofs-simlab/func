# Example usage

Before delving into the details of each feature in FunC, we provide several
examples of replacing a mathematical function with a LUT. We then generalize
this process and summarize it as a workflow diagram. After this, we give a
brief overview of the class structure of FunC.

[comment]: # (To learn more about optimizing LUTs, read the paper.)


## How to use FunC to replace a mathematical function with a LUT

The following example illustrates how to use FunC to build a cubic
LUT for an exponential integral over \f$[0.01,2]\f$ with step size \f$h=0.1\f$.

```cpp
#include <boost/math/special_functions/expint.hpp>
#include <func/func.hpp>
#include <iostream>

template <typename T> T f(T x){ return boost::math::expint(1,x); } // user's function
int main(){ // build an approximation of f over [0.01,2] with uniform stepsize h=0.1
  func::FunctionContainer<double> fc {FUNC_SET_F(f,double)};
  func::UniformExactInterpTable<3,double> LUT {fc, {0.01,2.0,0.1}};
  double x; std::cin >> x;
  std::cout << "f(" << x << ")~" << LUT(x) << std::endl; // print piecewise cubic approx.
}
```

We observe the following:
- The LUT is only used once.
- The approximation is built according to a step size. We do not currently know how much error the approximation has.
- It is impossible in principle to know what the user will input, so the subdomain \f$[0.01,2]\f$ is likely insufficient.

In this case, using a LUT is a poor choice because the overhead from
building the LUT is not balanced out by repeatedly calling the LUT, the LUT
introduces an unknown amount of error, and is only valid over a small subset of
the original domain of \f$f\f$. By default, `FunC`'s LUTs do not perform bounds
checking because doing so introduces nontrivial slowdown when calling
`operator()`. such, any user input outside the range \f$[0.01,2]\f$ is
undefined behavior. It is up to the user to guarantee this undefined behavior
is not possible. This is generally done by using a `LookupTable` container
(defined shortly), or making the interval \f$[a,b]\f$ larger.

We now consider the following improved example using the same \f$f\f$ as before.
```cpp
int main(){
  // build an approximation of f over [0.01,2] with uniform stepsize h=0.1
  func::FunctionContainer<double> fc {FUNC_SET_F(f,double)};
  func::FailureProofTable<func::UniformEqSpaceInterpTable<3,double>> lut {fc, {0.01,2.0,0.1}};
  // compute max error of LUT
  func::LookupTableGenerator<double> gen(fc,{});
  std::cout << "error=" << gen.error_of_table(lut,1.0) << "\n";

  // take two numbers from the user: max and nevals
  double max; std::cin >> max; int nevals; std::cin >> nevals;
  // print nevals random numbers in the range [0.01, max].
  std::mt19937 mt(0); std::uniform_real_distribution<double> unif(0.01,0.01+std::abs(max));
  for(int i=0; i<nevals; i++){
    double x = unif(mt);
    if(x < 2.0) std::cout << " f(" << x << ") ~ " << lut(x) << "\n";
    else std::cout << " f(" << x << ") = " << lut(x) << "\n";
  } std::cout << std::endl;
}
```

Sample input and output:

```
error=0.0216333
3 3 # user input
 f(1.78261) ~ 0.0663324
 f(2.53435) = 0.0238136
 f(2.57526) = 0.0225693
```


This example is better suited for a LUT because it can involve numerous
repeated applications. The undefined behavior is gone because the
`FailureProofTable` resorts back to \f$f\f$'s defining mathematical formula if
\f$x\f$ is out of the bounds of the `UniformExactInterpTable`. The
`LookupTableGenerator` provides an estimate of
\f[
	E(L) = \max_{x\in[a,b]}\frac{|f(x) - L(x)|}{a_{\mathrm{tol}} + r_{\mathrm{tol}}|f(x)|},
\f]
(by default with \f$a_{\mathrm{tol}}=r_{\mathrm{tol}}=1\f$ but these values are
adjustable). Whether this is an acceptable amount of error depends on the use
case. If the user provides \f$(a_{\mathrm{tol}},r_{\mathrm{tol}}) = (1,0)\f$,
then \f$E(L)\f$ is the absolute error of \f$f\f$. Similarly, if
\f$(a_{\mathrm{tol}},r_{\mathrm{tol}}) = (0,1)\f$, then \f$E(L)\f$ is the
relative error of \f$f\f$. The user can provide any positive values for
\f$a_{\mathrm{tol}}\f$ and \f$r_{\mathrm{tol}}\f$, and if they are both nonzero
then \f$E(L)<1\f$ does not necessarily imply \f$L\f$ satisfies _both_ the relative
and absolute error tolerances individually.


The following code shows an MWE of building a `DirectEvaluation` and a LUT of a special function.
A `DirectEvaluation` can record every argument passed to its `operator()` in
a histogram which is critical for determining useful bounds \f$a,b\f$ for a
LUT. A DirectEvaluation can easily simulate error in a LUT by perturbing its
arguments by \f$r_{\mathrm{tol}},a_{\mathrm{tol}}\f$. So, the
`DirectEvaluation` can return `rtol*R*f(x) + A*atol` where \f$R\f$ and \f$A\f$
are uniformly distributed random numbers in \f$[-1,1]\f$.

```cpp
/* User's function here. Some LUT types require derivatives of the user's
 * function, and this is provided though Boost's automatic differentiation
 * library. To use automatic differentiation, the definition of f must be
 * templated, and any function f calls must have overloads for Boost's autodiff_fvar */
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
     * instead where A,R are random numbers sampled from a uniformly distributed random variable over [-1,1] */
    std::cout << de(1.011) << "\n"; 
    std::cout << de << std::endl; // print histogram to stdout if FUNC_DEBUG is defined

    /* build a LUT of f over [0,2] with a step size of h=0.1. Each subinterval will use degree 3 Chebyshev interpolating polynomials */
    func::UniformChebyInterpTable<3,double> lut {fc, {0.0,2.0,0.1}};
    std::cout << lut(1.011) << "\n"; // return an approximation of f(1.011) with a piecewise cubic polynomial

    std::cout << lut << "\n"; // print info about lut
}
```

Instead of building a LUT according to a step size, it is better to
build a LUT according to tolerances for error that are as large as
possible for the user's purpose.
```cpp
int main(){
		auto tableTol = 1e-2;
    func::FunctionContainer<double> fc {FUNC_SET_F(f,double)};
    LookupTableGenerator<TYPE> gen(func_container,0.1,2.0);
    /* Use a LUT factory to build according to a string. Using a NonUniform LUT */
    auto lut = gen.generate_by_tol("NonUniformChebyInterpTable<3>",tableTol);
}
```

Any member function of the
`LookupTableGenerator` class that returns a 
`std::unique_ptr<LookupTable>` (`generate_by_step`, `generate_by_impl_size`,
`generate_by_tol`) can take an optional `std::string filename`. When given a
`filename`, that member function returns the result of
`generate_by_file(filename)` if `filename` exists.  Otherwise, that member
function saves its result to `filename` before returning. As such, the code
used to generate a LUT is automatically optimized for future runs of the user's
program.


### A general workflow

The following is the general workflow we suggest for replacing a
mathematical function with a LUT using `FunC`.
- The user identifies a mathematical function whose evaluation
  consumes a substantive proportion of total runtime. Further, the
  mathematical function itself must be complicated enough to warrant
  the use of a LUT. For example, it is unlikely that a LUT will speed
  up an elementary function such as \f$\sin(x)\f$, \f$e^x\f$, etc. because
  those functions have been continually optimized for decades. Good
  candiates for a LUT include deeply nested function compositions,
  long summations/products, and special functions. Also, the process
  of building a LUT is not without expense. `FunC` must either
  evaluate the user's function a set number of times or read the LUT
  from a file (both are potentially slow). So, the user's code must
  evaluate \f$f\f$ a sufficiently large number of times for construction
  of a LUT to be appropriate.
- The user determines an interval of approximation \f$[a,b]\f$ and tolerances
	for error in their LUT. LUTs over smaller domains and with coarser
	error tolerances perform faster and use less memory.
- The user should experiment with a wide variety of different LUTs
  in isolation, each constructed according to the user's tolerances for
	error. `FunC` provides the `LookupTableGenerator` to test each LUT
	in isolation.
- After determining 1--2 ideal LUTs, the user should still benchmark
	their code after using those LUTs. Use of a LUT necessarily increases
	memory usage and that can result in overall slowdown compared to the
	original code.

We note that the above process does not put any emphasis on the
specific approximation used on each subinterval. Again, most LUTs
perform similarly because they share much of the same source code. The
main determining factor for a LUT's performance is its _order_ of
convergence. For a user to determine a suitable order for their application,
they must experiment with several LUTs in isolation.




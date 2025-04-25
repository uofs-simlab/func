Lookup Tables over nonuniform partitions of [a,b]
---------------------------------------------------
`FunC` provides two methods for constructing a LUT over a nonuniform partition
of \f$[a,b]\f$. First, many of `FunC`'s LookupTable implementations have built-in
support for a nonuniform partition, but such a construction does not allow for
custom partitions of \f$[a,b]\f$ from the user. This way, we ensure `FunC` can
hash its nonuniform LUTs in \f$6\f$ FLOPs and zero comparisons. If a custom
partition is required then the other option is to use a `CompositeLookuptable`,
although the resulting LUT will be much slower.

Let \f$L_1,L_2,...,L_M\f$ be LUTs of \f$f\f$ over the pairwise disjoint
intervals \f$[a_1,b_1], [a_2,b_2], ..., [a_M,b_M]\f$, respectively (not
necessarily partitioning the domain of \f$f\f$). A `CompositeTable` of
\f$L_1,L_2,...,L_M\f$ is a piecewise function of the form
\f[
	C(x) =
	\cases{
		L_1(x), & if $a_1\leq x\leq b_1$, \cr
		L_2(x), & if $a_2\leq x\leq b_2$, \cr
		\hspace*{3mm} . & \hspace*{11mm} . \cr
		\hspace*{3mm} . & \hspace*{11mm} . \cr
		\hspace*{3mm} . & \hspace*{11mm} . \cr
		L_M(x), & if $a_M\leq x\leq b_M$, \cr
		f(x), & otherwise.
	}
\f]
A `CompositeLookupTable` allows users to build LUTs over
arbitrary nonuniform partitions. This is useful if \f$f\f$ has
discontinuities, \f$f\f$ is difficult to accurately approximate on some
subset of its domain, or the user's program requires that \f$f\f$ is exact
at certain other special points (roots, extrema, inflection points, etc). The
downside is that this requires \f$O(\log n)\f$ comparisons each time the class's
`operator()` is called.

We can reduce the relative error in a LUT by building a `CompositeLookupTable`
over \f$f\f$'s roots. Doing so with \f$f(x) = \ln|\Gamma(x)|\f$ over
\f$[0.1,3]\f$, \f$L=\texttt{UniformExactInterpTable<3>}\f$, and \f$30\f$
subintervals reduces \f$E(L)\f$ with \f$a_{\mathrm{tol}}=r_{\mathrm{tol}}=1\f$
from \f$1.19805\times10^{-4}\f$ to \f$5.7253\times10^{-6}\f$ (\f$21\f$ times
less error). 

As for the nonuniform LUTs, they tend to perform best when \f$f'\f$ is largest
at its endpoints \f$a,b\f$. For example, the nonuniform LUT will have almost
the exact same partition as a uniform LUT for the function \f$f(x)=e^{x^2}\f$
(because) \f$f'(a)=-f'(b)\approx10^{-10}\f$ is very small.
To remedy this issue, we can build a `CompositeLookupTable` over
\f$f\f$'s inflection points (extremum of \f$f'\f$) as shown in the following
figure. The constituent nonuniform LUTs use a nontrivial partition of
\f$[-5,5]\f$, and the overall composite LUT is \f$28\f$ times more accurate than a
single nonuniform LUT and has the same memory usage the other LUTs. 


Building a `CompositeLookupTable` of \f$e^{-10x^2}\f$ and including inflection points in the partition of \f$[-5,5]\f$:
```cpp
FunctionContainer<double> func_container{FUNC_SET_F(MyFunction,double)}; auto step = 0.05;
UniformExactInterpTable<3,double>       uniformlut(func_container, {min,max,step});
NonUniformExactInterpTable<3,double> nonuniformlut(func_container, {min,max,step});

/* Build a Composite LUT over the inflection points of exp(-10*x*x) with the same error
 * as uniformlut. Using rtol = atol = 1.0 with E(L) */
LookupTableGenerator<double> gen(func_container, min, max);
auto err = gen.error_of_table(uniformlut); auto a = 0.035;
CompositeLookupTable<double> nonunifcom(func_container, {
	/* {tableKey, left, right, atol, rtol}, */
	{"NonUniformExactInterpTable<3>",min,-1.0/sqrt(2.0*10.0),                       a*err,a*err},
	{"NonUniformExactInterpTable<3>",    -1.0/sqrt(2.0*10.0),1.0/sqrt(2.0*10.0),    a*err,a*err},
	{"NonUniformExactInterpTable<3>",                        1.0/sqrt(2.0*10.0),max,a*err,a*err},
});
/* Verify nonuniformlut and uniformlut are approx. equal and compare with the composite LUT */
std::cout << "Error in uniform LUT:    " << err << std::endl;
std::cout << "Error in nonuniform LUT: " << err << " + " 
	<< gen.error_of_table(nonuniformlut) - err << std::endl;
std::cout << "Error in nonuniform composite LUT: " << gen.error_of_table(nonunifcom) << std::endl;
std::cout << "Memory usage of uniform LUT: " << uniformlut.size() << std::endl;
std::cout << "Memory usage of nonuniform composite LUT: " << nonunifcom.size() << std::endl;
```

Output:

```
Error in uniform LUT:    1.05229e-06
Error in nonuniform LUT: 1.05229e-06 + -3.83676e-14
Error in nonuniform composite LUT: 3.52733e-08
Memory usage of uniform LUT: 6432
Memory usage of nonuniform composite LUT: 6496
```

The following figure shows how long it takes to apply each LUT from the
previous figure to a random vector of length \f$1\,000\,000\f$. We see that the
`CompositeLookupTable`'s `operator()` is about \f$12\f$ times slower than
individual LUTs. The `CompositeLookupTable`'s improvement in accuracy comes at
a cost.

Average time to apply the `operator()` of each LUT from the previous
figure ten times to a random vector of size \f$1\,000\,000\f$
```
----------------------------------------------------------------------------
Table input and output types: d -> d
Number of trials performed: 10
Number of evaluations used: 1 000 000
----------------------------------------------------------------------------
| LookupTable:      NonUniformExactInterpTable<3> -5 5 0.05 200
| Memory usage (B): 6432
| Timings:          Min 0.00310616s Max 0.00320847s Mean 0.00313021s
----------------------------------------------------------------------------
| LookupTable:      UniformExactInterpTable<3> -5 5 0.05 200
| Memory usage (B): 6432
| Timings:          Min 0.00259926s Max 0.00672951s Mean 0.0030255s
----------------------------------------------------------------------------
| LookupTable:      CompositeLookupTable -5 5 0.0154212 200
| Memory usage (B): 6496
| Timings:          Min 0.0371598s Max 0.0372521s Mean 0.0371988s
----------------------------------------------------------------------------
```

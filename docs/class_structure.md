### FunC's class structure

With this workflow in mind, we now present a brief summary of `FunC`'s
debugging tools and important classes. To visualize how each of these
classes relate to one another, `FunC`'s UML class diagram is provided. Each of
`FunC`'s classes are related to `LookupTable` because they either
implement `LookupTable`, encapsulate a `LookupTable`
implementation, or construct `LookupTable` implementations.

- Classes implementing the `LookupTable` interface implement a
	useful set of functions for approximating a mathematical function with a
	piecewise function. The most important member function of a
	`LookupTable` implementation is its `operator()` (because it
	returns approximations of \f$f(x)\f$).
- The `MetaTable` class provides all the mechanisms required to
	approximate a mathematical function with an array of
	`Polynomial`. `MetaTable` exists to reduce code
	redundancy and as such is templated on several parameters: the
	number of polynomial coefficients for each subinterval,
	`TIN`,`TOUT`, and whether the partition of \f$[a,b]\f$ is
	uniform. Currently, every class that _constructs_ a piecewise
	approximation of \f$f\f$ inherits from `MetaTable`.
- The `LookupTableGenerator` class uses the factory design
	pattern and provides several member functions for building any
	supported LUT according to a step size, data size, or tolerances for
	relative and absolute error.
- The `DirectEvaluation` class is used for profiling or
	debugging as per the preprocessor macro `FUNC_DEBUG`.
	When this macro is defined, it helps determine useful LUT bounds
	(by recording each argument it is given before returning \f$f(x)\f$)
	and tolerances for relative and absolute error (by perturbing its
	outputs).
- The `FailureProofTable` class passes each of its arguments \f$x\f$ to
	the LUT it encapsulates after checking whether \f$x\f$ is within its LUT's bounds.
	If \f$x\f$ is not within the LUT's bounds, then the `FailureProofTable`
	computes \f$f(x)\f$ using the defining mathematical formula for \f$f\f$. This
	makes LUTs safe and straightforward to incorporate into existing code,
	especially if it is impossible/impractical to ensure each argument
	lies within a LUT's bounds.
- The `ArgumentRecord` class only exists in `FunC` if the preprocessor macro
  `FUNC_DEBUG` is defined. If `FUNC_DEBUG` is defined, then every argument
	passed to a `DirectEvaluation` and any out of bounds arguments passed to a
	`FailureProofTable` are also passed to `ArgumentRecord` to save in a
	histogram before computing \f$f(x)\f$. When the destructor of an
	`ArgumentRecord` is called, it prints its histogram to a provided
	`std::ostream*` (but does nothing if the pointer is null).
- The `CompositeLookupTable` class builds a LUT of \f$f\f$ over several
  pairwise disjoint intervals. This enables users to build a LUT over custom
	partitions. When performing interval search, a `CompositeLookupTable` must
	perform binary search over a sorted tree of endpoints in \f$O(\log(N))\f$ time.
	If binary search fails, \f$f(x)\f$ is returned. This class is ideal for
	piecewise continuous functions and can interact nicely with nonuniform LUTs.
- The `LookupTableComparator` class can compare the time
	taken to apply a set of LUTs to a uniformly distributed random
	vector.

Every class in `FunC` with an `operator()` implements `LookupTable`. Only three
classes that implement `LookupTable` do not also inherit from `MetaTable`. Of
these, `FailureProofTable` and `CompositeLookupTable` are what we call
_LUT containers_. They resort back to the defining mathematical formula
for \f$f\f$ if their interval search fails. LUT containers notably slow
the LUT they encapsulate only because their interval search is slower. The only
other class implementing `LookupTable` that does not also inherit
`MetaTable` is the `DirectEvaluation` class. A `DirectEvaluation` does
not provide an approximation of \f$f\f$ for any
inputs. Rather, it evaluates \f$f\f$ using the defining mathematical formula. The
`DirectEvaluation` class is provided for the convenience of debugging
and profiling mathematical functions (depending on whether the preprocessor
macro `FUNC_DEBUG` is defined).



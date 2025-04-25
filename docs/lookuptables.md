### FunC's LookupTable implementations

The practical differences between any two LUTs is the error in the
approximation they use for each subinterval (the \f$p_k\f$ from
the equation on the main page) and the number of coefficients it uses to represent
its \f$p_k\f$. Categorizing every `LookupTable` implementation based on how
they compute their \f$p_k\f$ results in the following four distinct families and
two special cases. The following table shows each family of `LookupTable`
supported in `FunC` along with some properties.

|Name                 | Necessarily continuous? | Requires derivative? | NonUniform Support |
|---------------------|-------------------------|----------------------|--------------------|
|ChebyInterpTable<n>  | No        | No  | Yes |
|ExactInterpTable<n>  | \f$C^0\f$ | No  | Yes |
|TaylorTable<n>       | No        | Yes | Yes |
|PadeTable<m,n>       | No        | Yes | No  |
|CubicHermiteTable    | \f$C^1\f$ | Yes | Yes |
|LinearRawInterpTable | \f$C^0\f$ | No  | No  |

Notes:
- Every nonuniform LUT requires derivatives of \f$f\f$ to construct their nonuniform partition.
- An `ExactInterpTable<0>` is a piecewise constant LUT,
  so it is likely not continuous.


Overview of each LUTs
- By default, a `ChebyInterpTable<n>` takes each \f$p_k\f$ to
  be the polynomial of degree \f$n\f$ interpolating \f$f\f$ on \f$n+1\f$ Chebyshev
  nodes of the first kind in \f$[x_k,x_{k+1}]\f$. This class depends on
	Armadillo to solve Vandermonde systems. Armadillo only supports matrices
  with entries in `float` or `double`. `FunC` takes the
  working precision to be `double` because that is the most
  accurate type available to Armadillo. When using a `ChebyInterpTable<n>`
  with other types, `FunC` statically casts `TIN` or `TOUT` to
  `double`, solves the linear system with Armadillo, and then
  casts back to the original types. 
- An `ExactInterpTable<n>` takes each \f$p_k\f$ to be the
	polynomial of degree-\f$n\f$ interpolating \f$f\f$ on \f$n+1\f$ Chebyshev nodes
	of the _second kind_ in \f$[x_k,x_{k+1}]\f$ (these are
	equally spaced nodes for \f$n\in\{0,1,2\}\f$). These LUTs are
	deemed "exact" because their source code contains a hard-coded symbolic
	expression for \f$V^{-1}\f$. The coefficients are accurate to the
	precision of `TIN`. `ExactInterpTable<n>` is the only
	`LookupTable` implementation that can be built over any types where `TOUT`
	forms an approximate vector space over an approximate field `TIN`.
- A `TaylorTable<n>` takes each \f$p_k\f$ to be a degree \f$n\f$
	truncated Taylor series of \f$f\f$ centered at the midpoint of
	\f$[x_k,x_{k+1}]\f$. This class uses Boost's automatic differentiation
	library to compute derivatives from the source code defining \f$f\f$.
- A `PadeTable<m,n>` takes each \f$p_k\f$ to be the \f$[m/n]\f$
	Pad\'e approximant of \f$f\f$ centered at the midpoint of
	\f$[x_k,x_{k+1}]\f$. We require \f$m \geq n > 0\f$. This class
	depends on both Armadillo and Boost. As such, the working precision
	is `double`.

There are two classes that are not templated on an unsigned integer.
- A `CubicHermiteTable` takes each \f$p_k\f$ to be the cubic
  Hermite spline over the data
	\f[
		\{(x_k,f(x_k)), (x_k,f'(x_k)), (x_{k+1}, f(x_{k+1})), (x_{k+1},f'(x_{k+1}))\}.
	\f]
	These LUTs are \f$C^1[a,b]\f$ and depend on Boost to compute derivatives.
- A `UniformLinearRawInterpTable` builds the same \f$p_k\f$ as a
	`UniformExactInterpTable<1>` with the same parameters, but a
	`UniformLinearRawInterpTable` saves memory by only storing
	\f$f(x_0), f(x_1),...,f(x_n)\f$. Its `operator()` must then
	compute the two coefficients of \f$p_k\f$ from \f$f(x_k)\f$ and \f$f(x_{k+1})\f$
	before returning \f$p_k(x)\f$. For comparison, a `UniformExactInterpTable<1>`
	stores both coefficients of each \f$p_k\f$, so it uses approximately twice as
	much memory as a `UniformLinearRawInterpTable` with the same parameters.  The
	overhead in a `UniformLinearRawInterpTable`'s `operator()` is not large. There
	is no nonuniform variant of `UniformLinearRawInterpTable` because it does not
	allow for quick interval search.

We recommend users try to choose a particular family with properties that are
conducive to best approximate $f$. For example, if the LUT must be continuous,
then an `ExactInterpTable` or `CubicHermiteTable` are best.

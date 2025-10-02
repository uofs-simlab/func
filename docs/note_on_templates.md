### Note on templates

One can theoretically use LUTs with any types such that `TOUT` forms an
approximate vector space over `TIN` (addition and scalar multiplication
accurate to machine epsilon exists). This generality is possible for the
`ExactInterpTable` because the solution to its Vandermonde system is
hard-coded. We cannot currently offer this level of generality for `PadeTable`s
and `ChebyInterpTable`s because they depend on Armadillo to solve linear
systems of equations (possibly also LUTs requiring derivatives because it is
difficult to work out the theory in that case). Armadillo only supports
matrices with entries in `float` or `double` (which is typical
for high-performance linear algebra libraries). The generality allowed by
`FunC`'s templates is typical of header-only libraries, but `FunC` is not a
header-only library. Template values `TIN`\f$=\f$`TOUT`\f$=\f$`float` and
`TIN`\f$=\f$`TOUT`\f$=\f$`double` are explicitly instantiated for the
`LookupTableFactory`. These are then compiled into a dynamic library.
This way, the user can link their project with `libfunc.so` (avoiding the
increase in compile time from templates) if they use LUTs with the two most
common numeric types. Doing this has resulted in a \f$1.75\f$ times speedup compared
to headers only when compiling all the example code on our GitHub. This test
does not include the time taken to compile `libfunc.so`, but the user
need only compile `libfunc.so` once anyways. Linking with a dynamic
library makes it much easier to quickly experiment with different LUTs and
debug the user code. If user code instantiates other values of `TIN` and
`TOUT`, then those LUTs are compiled from scratch with those template values
at the same time as the user's code. Similarly, the rest of `FunC`'s code is
header-only.

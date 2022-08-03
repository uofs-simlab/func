/* Interpolate f over Chebyshev points then (optionally based on a template) refine this approximation using Remez' algorithm
 * Refining the points will generate a best possible polynomial approximation if
 * - f' and f'' are available and not too slow (note: they will be accurate enough b/c we use autodiff)
 * - will need to use a higher precision type (float128?) when doing any refinement (need to find global max error w/ Brent or Newton's)
 * - max absolute pointwise error is known from algorithm (but our LUTGenerator only cares about max relative pointwise error)
 * - For error approximation, the error computation from LookupTableGenerator could be made into its own class
 * */

/* Interpolate f over Chebyshev points then (optionally based on a template) refine this approximation using Remez' algorithm
 * Refining the points will generate a best possible polynomial approximation if
 * - f' and f'' are available and not too slow (note: they will be accurate enough b/c we use autodiff)
 * - will need to use a higher precision type (float128?) when doing any refinement (need to find global max error w/ Brent or Newton's)
 * - max absolute pointwise error is known from algorithm (but our LUTGenerator only cares about max relative pointwise error)
 * - For error approximation, the error computation from LookupTableGenerator could be made into its own class
 * */

/* function to return every chebyshev node in [a,b]
 * corresponding to the positive integer N */
//template <typename TIN, unsigned int N>
//std::array<TIN,N> cheb_nodes(TIN a, TIN b){
//  std::array<TIN,N> chebary;
//  for(unsigned int k=1; k<=N; i++)
//    chebary[k-1] = (a+b)/2 + (b-a)*cos(pi*(2*k-1)/(2*N))/2
//  return chebary;
//}

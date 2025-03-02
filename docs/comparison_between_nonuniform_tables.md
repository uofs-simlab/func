Lookup Tables over nonuniform partitions of $[a,b]$
---------------------------------------------------
`FunC` provides two methods for constructing a LUT over a nonuniform partition
of $[a,b]$. First, many of `FunC`'s LookupTable implementations have built-in
support for a nonuniform partition, but such a construction does not allow any
input from the user. This construction ensures `FunC` can hash its nonuniform
LookupTables in $6$ FLOPs and zero comparisons. Another option is to use a
`CompositeLookuptable`.

Let $L_1,L_2,\dotsc,L_M$ be LUTs of $f$ over the pairwise disjoint intervals
$[a_1,b_1], [a_2,b_2]$, $\dotsc$, $[a_M,b_M]$, respectively (not necessarily
partitioning the domain of $f$). A \texttt{CompositeTable} of
$L_1,L_2,\dotsc,L_M$ is a piecewise function of the form
\[
	C(x) =
	\begin{cases}
		L_1(x), & \text{if } a_1\leq x\leq b_1, \\
		L_2(x), & \text{if } a_2\leq x\leq b_2, \\
		\hspace*{3mm} \vdots & \hspace*{11mm} \vdots \\
		L_M(x), & \text{if } a_M\leq x\leq b_M, \\
		f(x), & \text{otherwise.}
	\end{cases}
\]
A \texttt{CompositeLookupTable} allows users to build LUTs over
arbitrary nonuniform partitions. This is useful if $f$ has
discontinuities, $f$ is difficult to accurately approximate on some
subset of its domain, or the user's program requires that $f$ is exact
at certain other special points (roots, extrema, inflection points, etc). The
downside is that this requires $O(\log n)$ comparisons each time the class's
`operator()` is called.

We can reduce the relative error in a LUT by building a
\texttt{CompositeLookupTable} over $f$'s roots. Doing so with $f(x) =
\ln|\Gamma(x)|$ over $[0.1,3]$, $L=\texttt{UniformExactInterpTable<3>}$,
and $30$ subintervals reduces $E(L)$ with $\abstol=\reltol=1$ from
$1.19805\times10^{-4}$ to $5.7253\times10^{-6}$ ($21$ times less error).

To remedy this issue, we can build a
\texttt{CompositeLookupTable} over $f$'s inflection points as shown in
Figure \ref{fig:InflectionCompositeLookupTable}. The constituent
nonuniform LUTs use a nontrivial partition of $[-5,5]$, and the
overall composite LUT is $28$ times more accurate than a single
nonuniform LUT. 
Figure \ref{fig:InflectionCompositeLookupTableTiming}
shows how long it takes to apply each LUT from Figure
\ref{fig:InflectionCompositeLookupTable} to a random vector of length
$1\,000\,000$. We see that the \texttt{CompositeLookupTable}'s
\texttt{operator()} is about $12$ times slower than individual~LUTs.

\begin{figure}[htbp]
	\caption{Building a \texttt{CompositeLookupTable} of $e^{-10x^2}$ and including inflection points in the partition of $[-5,5]$}
	\label{fig:InflectionCompositeLookupTable}
\footnotesize
\centering
\begin{lstlisting}
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
\end{lstlisting}
%/* Verify nonuniformlut and uniformlut are approx. equal and compare with the composite LUT */
%std::cout << "Error in uniform LUT:    " << err << std::endl;
%std::cout << "Error in nonuniform LUT: " << err << " + " 
%	<< gen.error_of_table(nonuniformlut) - err << std::endl;
%std::cout << "Error in nonuniform composite LUT: " << gen.error_of_table(nonunifcom) << std::endl;
%std::cout << "Memory usage of uniform LUT: " << uniformlut.size() << std::endl;
%std::cout << "Memory usage of nonuniform composite LUT: " << nonunifcom.size() << std::endl;
Output:
\begin{lstlisting}
Error in uniform LUT:    1.05229e-06
Error in nonuniform LUT: 1.05229e-06 + -3.83676e-14
Error in nonuniform composite LUT: 3.52733e-08
Memory usage of uniform LUT: 6432
Memory usage of nonuniform composite LUT: 6496
\end{lstlisting}
The LUTs \texttt{uniformlut} and \texttt{nonuniformlut} use almost
identical partitions of $[-5,5]$ because $e^{-10x^2}$ is approximately $0$ at
the endpoints $-5,5$. The composite LUT's uses the same amount of memory as
\texttt{uniformlut}, but its consituent nonuniform LUTs use nontrivial transfer
functions. Overall, it has $28$ times less error than \texttt{uniformlut}.
\end{figure}

\begin{figure}[htbp]
  \caption{Average time to apply the \texttt{operator()} of each LUT
    from Figure \ref{fig:InflectionCompositeLookupTable} ten times to
    a random vector of size $1\,000\,000$}
	\label{fig:InflectionCompositeLookupTableTiming}
\footnotesize
\centering
\begin{lstlisting}
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
\end{lstlisting}
\end{figure}



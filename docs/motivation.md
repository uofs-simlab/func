# Overview of FunC and its structure

Before delving into the details of each feature in FunC, we provide two minimal
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
`operator()` (see Section). Depending on the order
of the LUT, the minimum slowdown can be between \f$3\%\f$ to \f$33\%\f$ for the machine
in Figure. As such, any
user input outside the range \f$[0.01,2]\f$ is undefined behavior. It is up to the
user to guarantee this undefined behavior is not possible. This is generally
done by using a `LookupTable` container (defined shortly), or making the
interval \f$[a,b]\f$ larger.

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
\f$x\f$ is out of the bounds of the `UniformExactInterpTable`. 

The `LookupTableGenerator` provides an estimate of
\f[
	E(L) = \max_{x\in[a,b]}\frac{|f(x) - L(x)|}{a_{\mathrm{tol}} + r_{\mathrm{tol}}|f(x)|},
\f]
(with $\abstol=\reltol=1$). Whether this is an acceptable amount of error
depends on the use case. If the user provides $(\abstol,\reltol) = (1,0)$, then
$E(L)$ is the absolute error of $f$. Similarly, if $(\abstol,\reltol) = (0,1)$,
then $E(L)$ is the relative error of $f$. The user can provide any positive
values for $\abstol$ and $\reltol$, and if they are both nonzero then $E(L)<1$
does not necessarily imply $L$ satisfies \emph{both} the relative and absolute
error tolerances individually.


The following code shows an MWE of building a DirectEvaluation and a LUT of a special function.
A DirectEvaluation can record every argument passed to its operator() in
a histogram which is critical for determining useful bounds $a,b$ for a
LUT. A DirectEvaluation can easily simulate error in a LUT by perturbing its
arguments by $r_{\mathrm{tol}},a_{\mathrm{tol}}$. So, the DirectEvaluation can
return `rtol*R*f(x) + A*atol` where $R$ and $A$ are uniformly distributed
random numbers in $[-1,1]$.

```c++
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
build a LUT according to tolerance for error that are as large as
possible.
```c++
int main(){
    func::FunctionContainer<double> fc {FUNC_SET_F(f,double)};
    LookupTableGenerator<TYPE> gen(func_container,0.0,2.0);
    /* Use a LUT factory to build according to a string. Using a NonUniform LUT */
    auto lut = gen.generate_by_tol("NonUniformChebyInterpTable<3>",tableTol);
}
```


The following is the general workflow we suggest for replacing a
mathematical function with a LUT using \FunC. We use this workflow in
Section~\ref{sec:examples} to create a LUT for use with the Canadian
Hydrological Model.
\begin{enumerate}
\item The user identifies a mathematical function whose evaluation
  consumes a substantive proportion of total runtime. Further, the
  mathematical function itself must be complicated enough to warrant
  the use of a LUT. For example, it is unlikely that a LUT will speed
  up an elementary function such as $\sin(x)$, $e^x$, etc. because
  those functions have been continually optimized for decades. Good
  candiates for a LUT include deeply nested function compositions,
  long summations/products, and special functions. Also, the process
  of building a LUT is not without expense. \FunC{} must either
  evaluate the user's function a set number of times or read the LUT
  from a file (both are potentially slow). So, the user's code must
  evaluate $f$ a sufficiently large number of times for construction
  of a LUT to be appropriate.
\item The user determines an interval of approximation $[a,b]$ and tolerances
	for error in their LUT. As discussed in Subsection~\ref{sec:performance},
	LUTs over smaller \revFour{intervals} and coarser error tolerances perform
	faster and use less memory.
\item The user should experiment with a wide variety of different LUTs
  in isolation, each constructed according to the user's tolerances for
	error. \FunC{} provides the \texttt{LookupTableGenerator} to test each LUT
	in isolation.
\item After determining 1--2 ideal LUTs, the user should still benchmark
	their code after using those LUTs. Use of a LUT necessarily increases
	memory usage and that can result in overall slowdown compared to the
	original code.
\end{enumerate}
In Section~\ref{sec:performance}, we explore several other important
factors that determine the performance of a LUT. This includes more
reasons why a LUT could slow down the overall user code despite being
several times faster in isolation.

\begin{figure}[htbp]
\caption{Suggested \FunC{} workflow}
\label{fig:workflow}
\footnotesize
\centering
\begin{tikzpicture}[node distance=0.8cm]
	\node (start) [process] {Use a profiler to identify a bottleneck};
	\node (bound) [decision, below right=of start, text width=30mm]
		{Determine LUT \\ bounds $a,b$ and \\ tolerances $\abstol,\reltol$ with a \texttt{DirectEvaluation}};
	\node (failu) [process, above right=of bound, text width=30mm] {Use a \texttt{LookupTable} Container};
	\node (prove) [process, below right=of bound, text width=30mm] {Guarantee the suggested bounds are large enough};
	\node (gener) [process, below right=of failu, text width=30mm] {Generate several \texttt{LookupTable} implementations according to $\abstol$ and $\reltol$ with a \texttt{LookupTableGenerator}};
	\node (time) [process, below=of gener, text width=30mm] {Time each LUT with a \texttt{LookupTableComparator}};

	\draw [arrow] (start) -- (bound);
	\draw [arrow] (bound) -- (failu);
	\draw [arrow] (bound) -- (prove);
	\draw [arrow] (prove) -- (gener);
	\draw [arrow] (failu) -- (gener);
	\draw [arrow] (gener) -- (time);
\end{tikzpicture}

%Notes:
%\begin{enumerate}
%\item It is only safe to use a LUT outside of a LUT container (see
%	Section \ref{sec:containers}) if the user can guarantee their LUT's
%	bounds are large enough; otherwise, an out-of-bounds array access
%	can occur, resulting in a segmentation fault or bad output.
%\item The choice between which \texttt{LookupTable} container to use
%	(either \texttt{FailureProofTable} or \texttt{CompositeLookupTable})
%  depends on how the user's function is best represented.
%\item Users can avoid several function evaluations by reading their
%  LUT from a file.
%\end{enumerate}
\end{figure}

We note that the above process does not put any emphasis on the
specific approximation used on each subinterval. Again, most LUTs
perform similarly because they share much of the same source code. The
main determining factor for a LUT's performance is its \emph{order} of
convergence (explained more in Section~\ref{sec:math}). For a user to
determine a suitable order for their application, they must experiment
with several LUTs in isolation.

With this workflow in mind, we now present a brief summary of \FunC's
debugging tools and important classes. To visualize how each of these
classes relate to one another, \FunC's UML class diagram is provided
in Figure \ref{fig:class_diagram}. Each of \FunC's classes are related
to \texttt{LookupTable} because they either implement
\texttt{LookupTable}, encapsulate a \texttt{LookupTable}
implementation, or construct \texttt{LookupTable} implementations.
\begin{itemize}
\item Classes implementing the \texttt{LookupTable} interface implement a
	useful set of functions for approximating a mathematical function with a
	piecewise function. The most important member function of a
	\texttt{LookupTable} implementation is its \texttt{operator()} (because it
	returns approximations of $f(x)$).
  \item The \texttt{MetaTable} class provides all the mechanisms required to
    approximate a mathematical function with an array of
    \texttt{Polynomial}. \texttt{MetaTable} exists to reduce code
    redundancy and as such is templated on several parameters: the
    number of polynomial coefficients for each subinterval,
    \tin{},\tout{}, and whether the partition of $[a,b]$ is
    uniform. Currently, every class that \emph{constructs} a piecewise
    approximation of $f$ inherits from \texttt{MetaTable}.
\item The \texttt{LookupTableGenerator} class uses the factory design
	pattern and provides several member functions for building any
	supported LUT according to a step size, data size, or tolerances for
	relative and absolute error.
\item The \texttt{DirectEvaluation} class is used for profiling or
	debugging as per the preprocessor macro \texttt{FUNC\_DEBUG}.
	When this macro is defined, it helps determine useful LUT bounds
	(by recording each argument it is given before returning $f(x)$)
	and tolerances for relative and absolute error (by perturbing its
	outputs).
\item The \texttt{FailureProofTable} class passes each of its arguments $x$ to
	the LUT it encapsulates after checking whether $x$ is within its LUT's bounds.
	If $x$ is not within the LUT's bounds, then the \texttt{FailureProofTable}
	computes $f(x)$ using the defining mathematical formula for $f$. This
	makes LUTs safe and straightforward to incorporate into existing code,
	especially if it is impossible/impractical to ensure each argument
	lies within a LUT's bounds.
\item The \texttt{ArgumentRecord} class only exists in \FunC{}
	if the preprocessor macro \texttt{FUNC\_DEBUG} is defined. If 
	\texttt{FUNC\_DEBUG} is defined, then every argument passed to a
	\texttt{DirectEvaluation} and any out of bounds arguments passed to a
	\texttt{FailureProofTable} are also passed to \texttt{ArgumentRecord}
	to save in a histogram before computing $f(x)$. When the destructor of
	an \texttt{ArgumentRecord} is called, it prints its histogram to
	a provided \texttt{std::ostream*} (but does nothing if the pointer is
	null).
\item The \texttt{CompositeLookupTable} class builds a LUT of $f$ over several
	pairwise disjoint intervals. This enables users to build a LUT over custom
	partitions. When performing interval search, a \texttt{CompositeLookupTable}
	must perform binary search over a sorted tree of endpoints in $O(\log(N))$
	time. If binary search fails, $f(x)$ is returned. This class is ideal for
	piecewise continuous functions and can interact nicely with nonuniform LUTs.
\item The \texttt{LookupTableComparator} class can compare the time
	taken to apply a set of LUTs to a uniformly distributed random
	vector.
\end{itemize}
Every class in \FunC{} with an \texttt{operator()} implements
\texttt{LookupTable}. Only three classes that
implement\texttt{LookupTable} do not also inherit from
\texttt{MetaTable}. Of these, \texttt{FailureProofTable} and
\texttt{CompositeLookupTable} are what we call \emph{LUT
  containers}. They resort back to the defining mathematical formula
for $f$ if their interval search fails. LUT containers notably slow
the LUT they encapsulate only because their interval search is slower
(which we discuss further in Section~\ref{sec:performance}). The only
other class implementing \texttt{LookupTable} that does not also
inherit \texttt{MetaTable} is the \texttt{DirectEvaluation} class. A
\texttt{DirectEvaluation} does not provide an approximation of $f$ for
any inputs. Rather, it evaluates $f$ using the defining mathematical
formula. The \texttt{DirectEvaluation} class is provided for the
convenience of debugging and profiling mathematical functions
(depending on whether the preprocessor macro \texttt{FUNC\_DEBUG} is
defined). To learn more about how to use \FunC, see the user manual
attached to \FunC's Github repo \cite{Func2023web}.



\begin{figure}[htbp]
\caption{\FunC's Class Diagram}
\label{fig:class_diagram}
\footnotesize
\begin{tikzpicture}[node distance=0.5cm]
	\node (lutint) [class, text width=2cm,fill=blue!30] {\texttt{<{}<Interface>{}> LookupTable}};

	\node (metata) [class, above=2.3cm of lutint,fill=green!50] {\texttt{MetaTable}};

	\node (eqspac) [class, left=of metata, xshift=-0.2cm,fill=green!30] {\texttt{ExactInterpTable}};
	\node (chebyt) [class, above=0.2cm of eqspac,xshift=0.5cm,fill=green!30] {\texttt{ChebyInterpTable}};
	\node (linraw) [class, above=0.2cm of chebyt,xshift=0.7cm,fill=green!30] {\texttt{LinearRawInterpTable}};
	\node (hermit) [class, right=of metata,xshift=0.2cm,fill=green!30] {\texttt{CubicHermiteTable}};
	\node (taylor) [class, above=0.2cm of hermit,xshift=-0.7cm,fill=green!30] {\texttt{TaylorTable}};
	\node (padeta) [class, above=0.2cm of taylor,xshift=-0.9cm,fill=green!30] {\texttt{PadeTable}};
	\node (polyno) [class, below left=of metata,fill=green!15] {\texttt{Polynomial}};
	\node (transf) [class, below right=of metata,fill=green!15] {\texttt{TransferFunction}};

	\node (failur) [class, left=1.2cm of lutint,fill=orange!30] {\texttt{FailureProofTable}};
	\node (direct) [class, above=of failur, xshift=0.3cm,fill=orange!30] {\texttt{DirectEvaluation}};
	\node (compos) [class, below=of failur, xshift=0.3cm,fill=orange!30] {\texttt{CompositeLookupTable}};

	\node (argrec) [class, left=of direct, yshift=-0.2cm,fill=orange!15] {\texttt{ArgumentRecord}};

	\node (lutgen) [class, right=1.6cm of lutint,fill=red!30] {\texttt{LookupTableGenerator}};
	\node (lutcom) [class, below=of lutgen, xshift=-0.3cm,fill=red!30] {\texttt{LookupTableComparator}};
	\node (lutfac) [class, above=of lutgen=-0.6cm,fill=red!15] {\texttt{LookupTableFactory}};

	

	\draw [thick,-{Triangle[open]}] (eqspac) -- (metata);
	\draw [thick,-{Triangle[open]}] (chebyt.south east) -- (metata.north west);
	\draw [thick,-{Triangle[open]}] (linraw.south east) -- (metata);
	\draw [thick,-{Triangle[open]}] (hermit) -- (metata);
	\draw [thick,-{Triangle[open]}] (taylor.south west) -- (metata.north east);
	\draw [thick,-{Triangle[open]}] (padeta) -- (metata);
	\draw [thick,-{Diamond[open]}] (polyno.north east) -- (metata.south west);
	\draw [thick,-{Diamond[open]}] (transf.north west) -- (metata.south east);

	\draw [thick,-{Diamond[open]}] (argrec) -- (direct.west);
	\draw [thick,-{Diamond[open]}] (argrec) -- (failur.north west);

	\draw [thick,-{Diamond[open]}] (lutint) -- (failur.south east);
	\draw [thick,-{Diamond[open]}] (lutint) -- (compos.south east);

	\draw [thick,-{Diamond[open]}] (lutfac) -- (lutgen);
	\draw [thick,-{Diamond[open]}] (lutint.south east) -- (lutcom.north west);

	\draw [dashed, thick, ->] (metata) -- (lutint);
	\draw [dashed, thick, ->] (direct.south east) -- (lutint.north west);
	\draw [dashed, thick, ->] (failur.north east) -- (lutint);
	\draw [dashed, thick, ->] (compos.north east) -- (lutint.south west);

	\draw [thick, -] (lutgen) -- node[text width=1cm] {factory pattern}(lutint);
\end{tikzpicture}

\raggedright
Arrow legend:
\begin{itemize}
	\item Arrows with a triangular tip mean "inherits." For example, \texttt{ChebyInterpTable} inherits \texttt{MetaTable}.
	\item Arrows with a diamond tip mean "is a member variable of." For example, \texttt{Polynomial} is a member variable of \texttt{MetaTable}.
	\item Dotted arrows mean "implements the interface." For example, \texttt{MetaTable} implements the interface \texttt{LookupTable}.
\end{itemize}
\end{figure}


Currently, there are $45$ supported classes inheriting from
\texttt{MetaTable} for every combination of \tin{} and \tout{}. Of
these, $22$ \texttt{LookupTable} implementations can be built over a
nonuniform partition as described in Subsection \ref{sec:meshes}
(specified by prepending the names with \texttt{Uniform} or
\texttt{NonUniform}). Template meta-programming enables \FunC{} to use
the same source code to define the \texttt{operator()} of most
\texttt{LookupTable} implementations.  Because of this, the practical
differences between any two LUTs is the error in their~$p_k$ from
equation~\eqref{eqn:lut} and the number of coefficients it uses to
represent its $p_k$. Categorizing every \texttt{LookupTable}
implementation based on how they compute their $p_k$ results in the
following four distinct families and two special cases. Figure
\ref{fig:families} shows each family of \texttt{LookupTable} supported
in \FunC~\version{} along with some properties.

\begin{figure}[htbp]
\caption{Families of LUT Families}
\label{fig:families}
\footnotesize
\centering
\begin{tabular}{|c|c|c|c|}
	\hline
	Name & \unaligned{Necessarily \\ continuous?} & \unaligned{Requires \\ derivative?} & \unaligned{NonUniform \\ Support?} \\
	\hline
	ChebyInterpTable<n> & No & No & Yes \\
	\hline
	ExactInterpTable<n> & $C^0$ & No & Yes \\
	\hline
	TaylorTable<n> & No & Yes & Yes \\
	\hline
	PadeTable<m,n> & No & Yes & No \\
	\hline
	CubicHermiteTable & $C^1$ & Yes & Yes \\
	\hline
	LinearRawInterpTable & $C^0$ & No & No \\
	\hline
\end{tabular}

Notes:
\begin{itemize}
%\item To generate a LUT that requires derivatives of $f$, the
%  definition of $f$ must have overloads for Boost's
%  \texttt{autodiff\_fvar}.
\item Every nonuniform LUT requires derivatives of $f$.
\item An \texttt{ExactInterpTable<0>} is a piecewise constant LUT,
  so it is likely not continuous.
%\item \FunC{} only uses Armadillo \cite{Arma1} to build
%  \texttt{ChebyInterpTable}s and \texttt{PadeTable}s. The highest
%  precision type supported by Armadillo is \texttt{double}, so
%  \texttt{ChebyInterpTable}s require that \tin{} can be cast to
%  \texttt{double}, and \texttt{PadeTable}s require that both \tin{}
%  and \tout{} can be cast to \texttt{double}.
\end{itemize}
\end{figure}

\begin{itemize}
\item By default, a \texttt{ChebyInterpTable<n>} takes each $p_k$ to
  be the polynomial of degree $n$ interpolating $f$ on $n+1$ Chebyshev
  nodes of the first kind in $[x_k,x_{k+1}]$. This class depends on
	Armadillo to solve Vandermonde systems. Armadillo only supports matrices
  with entries in \texttt{float} or \texttt{double}. \FunC{} takes the
  working precision to be \texttt{double} because that is the most
  accurate available type. When using a \texttt{ChebyInterpTable<n>}
  with other types, \FunC{} statically casts \tin{} or \tout{} to
  \texttt{double}, solves the linear system with Armadillo, and then
  casts back to the original types. 
	%\texttt{ChebyInterpTable<n>} is
  %the only \texttt{LookupTable} implementation that uses the
  %\texttt{special\_points} field in the \texttt{LookupTableParameters}
  %struct. Given any special points of $f$ in $[x_k,x_{k+1}]$
  %(e.g., roots, critical points, inflection points, etc.), the nodes on
  %some subintervals are perturbed to preserve these points to high
  %accuracy.
	\item An \texttt{ExactInterpTable<n>} takes each $p_k$ to be the
		polynomial of degree-$n$ interpolating $f$ on $n+1$ Chebyshev nodes
		of the \emph{second kind} in $[x_k,x_{k+1}]$ (recall that these are
		equally spaced nodes for $n\in\{0,1,2\}$). These LUTs are
		\texttt{Exact} because their source contains a hard-coded symbolic
		expression for $V^{-1}$. The coefficients are accurate to the
		precision of \tin{}. \texttt{ExactInterpTable<n>} is the only
		\texttt{LookupTable} implementation where \tin{} and \tout{} are
		guaranteed to be usable whenever $\tout{}$ forms an approximate
		vector space over an approximate field $\tin{}$.
	\item A \texttt{TaylorTable<n>} takes each $p_k$ to be a degree $n$
		truncated Taylor series of $f$ centered at the midpoint of
		$[x_k,x_{k+1}]$. This class uses Boost's automatic differentiation
		library to compute derivatives from the source code defining $f$.
	\item A \texttt{PadeTable<m,n>} takes each $p_k$ to be the $[m/n]$
		Pad\'e approximant of $f$ centered at the midpoint of
		$[x_k,x_{k+1}]$. We require $m \geq n > 0$. This class
		depends on both Armadillo and Boost. As such, the working precision
		is \texttt{double}.
\end{itemize}
The two special cases are as follows.
\begin{enumerate}
\item A \texttt{CubicHermiteTable} takes each $p_k$ to be the cubic
  Hermite spline over the data
	\[
		\{(x_k,f(x_k)), (x_k,f'(x_k)), (x_{k+1}, f(x_{k+1})),
		(x_{k+1},f'(x_{k+1}))\}.
	\]
	These LUTs are $C^1[a,b]$ and depend on Boost to compute derivatives.
\item A \texttt{UniformLinearRawInterpTable} builds the same $p_k$ as a
	\texttt{UniformExactInterpTable<1>} with the same parameters, but a
	\texttt{UniformLinearRawInterpTable} saves memory by only storing
	$f(x_0), f(x_1),\dotsc,f(x_n)$. Its \texttt{operator()} must then
	compute the two coefficients of $p_k$ from $f(x_k)$ and $f(x_{k+1})$
	before returning $p_k(x)$. For comparison, a
	\texttt{UniformExactInterpTable<1>} stores both coefficients of each
	$p_k$, so it uses approximately twice as much memory as a
	\texttt{UniformLinearRawInterpTable} with the same parameters.
	In Subsection \ref{sec:meshes}, we compare each of the linear
	interpolation LUTs provided by \FunC{} (uniform, nonuniform, and raw).
	Figure \ref{fig:linear-error} shows that the overhead in a
	\texttt{UniformLinearRawInterpTable}'s \texttt{operator()} is not
	large. There is no nonuniform variant of
	\texttt{UniformLinearRawInterpTable} because it does not allow for
	quick interval search.
\end{enumerate}
We recommend that users try to choose a particular family with
properties that are conducive to best approximate $f$. For example, if
the LUT must be continuous, then an \texttt{ExactInterpTable} or
\texttt{CubicHermiteTable} are best.




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

We now consider the main function of the following improved example using the same \f$f\f$ as before.
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


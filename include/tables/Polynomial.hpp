/* - polynomial is essentially a typedef for std::array<TOUT,N> but polynomials are sometimes aligned
 * - Our convention for writing polynomials is:
 *     p(x) = m_table[x0].coefs[0] + m_table[x0].coefs[1]*x + ... + m_table[x0].coefs[N-1]*x^(N-1)
 * Users should perform Horner's style evaluations
 *
 * Polynomials could store all sorts of things:
 * - 2D LUTs may have coefs for x & y dimensions of each subsquare,
 * - LUTs could store derivative coefs */

#include <string>
#include <ostream> // operator<<

#pragma once

namespace func {

static constexpr unsigned int alignments[] = {0,1,2,4,4,8,8,8,8,16,16,16,16,16,16,16,16};

template <typename TOUT, unsigned int N, bool B> struct polynomial_helper;

template <typename TOUT, unsigned int N>
struct alignas(sizeof(TOUT)*alignments[N]) polynomial_helper<TOUT,N,true> {
  constexpr unsigned int size() const noexcept { return N; }
  TOUT coefs[N];
};

template <typename TOUT, unsigned int N>
struct polynomial_helper<TOUT,N,false> {
  constexpr unsigned int size() const noexcept { return N; }
  TOUT coefs[N];
};

template <typename TOUT, unsigned int N>
using polynomial = polynomial_helper<TOUT,N,std::is_floating_point<TOUT>::value>;


/* methods for computing with polynomials */
constexpr unsigned int factorial(unsigned int n){
  if(n == 0u) return 1u;
  return n*factorial(n-1u);
}

constexpr unsigned int permutation(unsigned int n, unsigned int k){
  if(k == 0u) return 1u;
  return n*permutation(n-1u,k-1u);
}

template <unsigned int N, typename TOUT, typename TIN = TOUT>
TOUT polynomial_diff(polynomial<TOUT,N> p, TIN x, unsigned s){
  TOUT sum = static_cast<TOUT>(0);
  for(unsigned int k=N; k>s; k--)
    sum = p.coefs[k-1]*permutation(k-1,s) + sum*x;
  return sum;
}

/* convenient debugging method for printing a polynomial */
template <unsigned int N, typename TOUT>
std::string polynomial_print(const polynomial<TOUT,N>& p){
  std::string sum = "";
  for(unsigned int k=N; k>0; k--)
    sum = sum + std::to_string(p.coefs[k-1]) + "x^" + std::to_string(k);
  return sum;
}

/* print basic info about a polynomials */
template <unsigned int N, typename TOUT>
std::ostream& operator<<(std::ostream& out, const polynomial<TOUT,N>& p){
  for(unsigned int k=N; k>1; k--)
    out << p.coefs[k-1] << "x^" << k << " + ";
  out << p.coefs[0] << "\n";
  return out;
}

/* wraps operator<< */
template <unsigned int N, typename TOUT>
inline std::string to_string(const polynomial<TOUT,N>& p) {
  std::ostringstream ss;
  ss << p;
  return ss.str();
}

}

/** \file Polynomial.hpp
 * \brief Define the polynomial class and some convenience functions (eg. Horner's eval, derivatives, printing)
 * \defgroup Polynomial
 * */

#include <string>
#include <ostream> // operator<<

#pragma once

namespace func {

static constexpr unsigned int alignments[] = {0,1,2,4,4,8,8,8,8,16,16,16,16,16,16,16,16};

/** \brief A typedef for std::array<TOUT,N> but polynomial arrays are sometimes aligned (always aligned for float or double)
 *  \ingroup Polynomial
 * 
 * \note
 * - Our convention for writing polynomials is:
 *   \f[p(x) = m_table[x0].coefs[0] + m_table[x0].coefs[1]*x + ... + m_table[x0].coefs[N-1]*x^{N-1}\f]
 *
 *  \tparam B determines whether an array of polynomials over TOUT are aligned with alignas(sizeof(TOUT)*alignments[N])
 *
 * Polynomials could store all sorts of things:
 * - 3D LUTs may have coefs for x & y dimensions of each subrectangle,
 * - LUTs could store derivative coefs */
template <typename TOUT, unsigned int N, bool B> struct polynomial_helper;

/** \brief Arrays of this type of polynomial are aligned
 * \ingroup Polynomial */
template <typename TOUT, unsigned int N>
struct alignas(sizeof(TOUT)*alignments[N]) polynomial_helper<TOUT,N,true> {
  constexpr unsigned int size() const noexcept { return N; }
  TOUT coefs[N];
};

/** \brief Arrays of this type of polynomial are not aligned
 *  \ingroup Polynomial */
template <typename TOUT, unsigned int N>
struct polynomial_helper<TOUT,N,false> {
  constexpr unsigned int size() const noexcept { return N; }
  TOUT coefs[N];
};

/** \ingroup Polynomial */
template <typename TOUT, unsigned int N>
using polynomial = polynomial_helper<TOUT,N,std::is_floating_point<TOUT>::value>;


constexpr unsigned int factorial(unsigned int n){
  if(n == 0u) return 1u;
  return n*factorial(n-1u);
}

constexpr unsigned int permutation(unsigned int n, unsigned int k){
  if(k == 0u) return 1u;
  return n*permutation(n-1u,k-1u);
}

/** \brief Compute p^(s)(x), the sth derivative of p at x.
 *
 * \note p cannot be empty
 * */
template <unsigned int N, typename TOUT, typename TIN = TOUT>
inline TOUT polynomial_diff(polynomial<TOUT,N> p, TIN x, unsigned s){
  //TOUT sum = static_cast<TOUT>(0.0);
  //TOUT com = static_cast<TOUT>(0.0);
  //for(unsigned int k=N; k>s; k--){
  //  TOUT y = p.coefs[k-1]*permutation(k-1,s) - com;
  //  TOUT t = y + sum*x;
  //  com = (t - sum*x) - y;
  //  sum = t;
  //}
  //return sum;
  //for(unsigned int k=N; k>s; k--){
  //  sum = p.coefs[k-1]*permutation(k-1,s) + sum*x;
  //}

  TOUT sum = static_cast<TIN>(permutation(N-1,s))*p.coefs[N-1];
  for(unsigned int k=N-1; k>s; k--){
    sum *= x;
    sum += static_cast<TIN>(permutation(k-1,s))*p.coefs[k-1];
  }
  return sum;
}

/** \brief Given a polynomial p:[a,b]->\R, compute the coefficients of
 *   q:[c,d]->\R such that q(x) = p( [(b-a)x + (ad-bc)]/(d-c) ) by
 *   expanding p in a Taylor series.
 *   
 *   This is used all over FunC (for example, special case for rightmost
 *   interval, Nonuniform LUTs, and Taylor/Pade tables). Optimizations
 *   are very welcome!
 * */
template <unsigned int N, typename TOUT, typename TIN = TOUT>
inline polynomial<TOUT,N> taylor_shift(polynomial<TOUT,N> p, TIN a, TIN b, TIN c, TIN d){
  polynomial<TOUT,N> q = p;
  for(unsigned int k=0; k<N; k++)
    q.coefs[k] = polynomial_diff(p, (a*d-b*c)/(d-c), k)*static_cast<TIN>(pow((b-a)/(d-c),k))/static_cast<TIN>(factorial(k));
  return q;
}

/** \brief Compute p(x)
 *
 * \note p cannot be empty
 * */
template <unsigned int N, typename TOUT, typename TIN = TOUT>
inline TOUT eval(polynomial<TOUT,N> p, TIN x){
  //TOUT sum = static_cast<TOUT>(0);
  //for(unsigned int k=N; k>0; k--)
  //  sum = p.coefs[k-1] + sum*x;
  //return sum;
  TOUT sum = p.coefs[N-1];
  for(unsigned int k=N-1; k>0; k--){
    sum *= x;
    sum += p.coefs[k-1];
  }
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

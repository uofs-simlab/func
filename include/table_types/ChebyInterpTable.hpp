/* 1st to 7th degree polynomial interpolation LUT over Chebyshev nodes (precomputed coefficients)
 *  t_s = (a+b)/2 + (b-a)\cos(\frac{2s-1}{2n}\pi)/2, \quad s=1,\dotsc,n.
 *
 * Usage example for a 4th degree interpolant:
 *   ChebyInterpTable<double,double,4> look(&function,LookupTableParameters<double>{0,10,0.0001});
 *   double val = look(0.87354);
 *
 *  Notes:
  - This class only works if TOUT and TIN can both be cast to double. 
    Armadillo Mat<T>'s `is_supported_elem_type<T>` will only let us do arithmetic
    with float or double (not even long double!) and arma::field is useless.
    TODO use a symbolic linear algebra library instead...

  - the template implementation is only registered for N=1,2,3,4,5,6,7
    but users can manually construct this class with larger N

    TODO
    - take a continuity parameter (use 2 nodes for C0, clamped spline for C1)
    - replace nearest points with user-provided roots of the function
 */


#pragma once
#include "MetaTable.hpp"
#include "config.hpp"
#include <stdexcept>

#ifdef FUNC_USE_ARMADILLO
#include <armadillo>
#endif

namespace func {

template <unsigned int N, typename TIN, typename TOUT=TIN, GridTypes GT=GridTypes::UNIFORM>
class ChebyInterpTable final : public MetaTable<N+1,TIN,TOUT,GT>
{
  INHERIT_META(N+1,TIN,TOUT,GT);
public:
  ChebyInterpTable() = default;
  ChebyInterpTable(const MetaTable<N+1,TIN,TOUT,GT>& L): MetaTable<N+1,TIN,TOUT,GT>(L) {}

  // build the LUT from scratch or look in filename for an existing LUT
  ChebyInterpTable(const FunctionContainer<TIN,TOUT>& func_container, const LookupTableParameters<TIN>& par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<N+1,TIN,TOUT,GT>(func_container, par, jsonStats) {
#ifndef FUNC_USE_ARMADILLO
    if(jsonStats.empty())
      throw std::invalid_argument("Error in func::ChebyTable: Chebyshev LUTs need Armadillo to be generated but Armadillo is not available");
#else
    if(!jsonStats.empty())
      return; // all our work is already done

    /* Base class default variables */
    m_name = grid_type_to_string<GT>() + "ChebyInterpTable<" + std::to_string(N) + ">";
    m_numTableEntries = m_numIntervals+1;
    m_order = N+1; // N is the degree of the polynomial interpolant so the order is N+1
    m_dataSize = sizeof(m_table[0]) * m_numTableEntries;

    auto fun = func_container.standard_fun;

    /* Thoughts:
     * - Vandermonde system for cheby nodes over [0,1] has condition number about 1000 times larger than [-1,1]
     *   but we don't do that because the formula we currently use for shifting polynomials is too poorly conditioned.
     * - Using inv here produces objectively worse results! however, Armadillo (and every other high performance numeric
     *   linear algebra library Shawn has looked into) does not have a solver for general vector spaces.
     *   That is, most libraries cannot solve Ax=b where the entries of A are not float or double and the entries of
     *   b aren't the same type as the entries of A.
     * - Directly solve the system when using a type supported by armadillo
     * */
    arma::mat Van = arma::ones(N+1, N+1);
    Van.col(1) = (1 + arma::cos(arma::datum::pi*(2*arma::linspace(1,N+1,N+1)-1) / (2*(N+1))))/2;
    for(unsigned int i=2; i<N+1; i++)
      Van.col(i) = Van.col(i-1) % Van.col(1); // the % does elementwise multiplication

    if(!std::is_floating_point<TOUT>::value)
      Van = arma::inv(Van);

    /* Allocate and set table */
    m_table.reset(new polynomial<TOUT,N+1>[m_numTableEntries]);
    FUNC_BUILDPAR
    for(unsigned int ii=0;ii<m_numTableEntries-1;++ii) {
      TIN x;
      TIN h = m_stepSize;
      // (possibly) transform the uniform grid into a nonuniform grid
      FUNC_IF_CONSTEXPR(GT == GridTypes::UNIFORM)
        x = m_minArg + ii*m_stepSize;
      else{
        x = m_transferFunction(m_minArg + ii*m_stepSize);
        h = m_transferFunction(m_minArg + (ii+1)*m_stepSize) - x;
      }

      // build the vector of coefficients from function values
      auto a = static_cast<double>(x);
      auto b = static_cast<double>(x+h);
      arma::vec xvec = (a+b)/2 + (b-a)*arma::cos(arma::datum::pi*(2*arma::linspace(1,N+1,N+1)-1)/(2*(N+1)))/2;

      FUNC_IF_CONSTEXPR(std::is_floating_point<TOUT>::value){
        arma::vec y(N+1);
        for (unsigned int k=0; k<N+1; k++)
          y[k] = fun(static_cast<TIN>(xvec[k]));

        y = arma::solve(Van, y);
        for(unsigned int k=0; k<N+1; k++)
          m_table[ii].coefs[k] = y[k];
      }else{
        std::vector<TOUT> y(N+1);
        for (unsigned int k=0; k<N+1; k++)
          y[k] = fun(static_cast<TIN>(xvec[k]));

        for(unsigned int k=0; k<N+1; k++){
          m_table[ii].coefs[k] = static_cast<TIN>(Van(k,0))*y[0];
          for(unsigned int s=1; s<N+1; s++){
            m_table[ii].coefs[k] += static_cast<TIN>(Van(k,s))*y[s];
          }
        }
      }

      /* TODO This formula is too unstable for this table type as given in this form when N>2 and h is small. */
      FUNC_IF_CONSTEXPR(GT == GridTypes::NONUNIFORM){
        auto p = m_table[ii];
        for(unsigned int k=0; k<N+1; k++)
          m_table[ii].coefs[k] = polynomial_diff(p,-x/h,k)/static_cast<TIN>(pow(h,k))/static_cast<TIN>(factorial(k));
      }
    }
    // special case to make lut(tableMaxArg) work
    m_table[m_numTableEntries-1].coefs[0] = fun(m_tableMaxArg);
    for (unsigned int k=1; k<N+1; k++)
      m_table[m_numTableEntries-1].coefs[k] = static_cast<TIN>(0)*m_table[m_numTableEntries-1].coefs[0];
#endif
  }
};

// define friendlier names
template <unsigned int N, typename TIN, typename TOUT=TIN>
using UniformChebyInterpTable = ChebyInterpTable<N,TIN,TOUT,GridTypes::UNIFORM>;
template <unsigned int N, typename TIN, typename TOUT=TIN>
using NonUniformChebyInterpTable = ChebyInterpTable<N,TIN,TOUT,GridTypes::NONUNIFORM>;

} // namespace func

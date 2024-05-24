#pragma once
#include "MetaTable.hpp"
#include "config.hpp"
#include <stdexcept>
#include <boost/math/special_functions/binomial.hpp>

#ifdef FUNC_USE_ARMADILLO
#include <armadillo>
#endif

namespace func {

/** \brief LUT using degree 1 to 7 polynomial interpolation over Chebyshev nodes on each subinterval (all coefficients are precomputed)
 *  \ingroup MetaTable
 *
 *  \f[ t_s = (a+b)/2 + (b-a)\cos(\frac{2s-1}{2n}\pi)/2, \quad s=1,\dotsc,n. \f]
 *
 * \code{.cpp}
 * // ChebyInterpTable does not require Boost's autodiff? TODO it might in the future?
 * template <typename T>
 * T foo(T x){ return x; }
 *
 * int main(){
 *   auto min = 0.0; auto max = 10.0; auto step = 0.0001;
 *   UniformChebyInterpTable<1,double> L1({FUNC_SET_F(foo,double)}, {min, max, step}); // degree 1
 *   UniformChebyInterpTable<2,double> L2({FUNC_SET_F(foo,double)}, {min, max, step}); // degree 2
 *   UniformChebyInterpTable<3,double> L3({FUNC_SET_F(foo,double)}, {min, max, step});
 *   UniformChebyInterpTable<4,double> L4({FUNC_SET_F(foo,double)}, {min, max, step});
 *   UniformChebyInterpTable<5,double> L5({FUNC_SET_F(foo,double)}, {min, max, step});
 *   UniformChebyInterpTable<6,double> L6({FUNC_SET_F(foo,double)}, {min, max, step});
 *   UniformChebyInterpTable<7,double> L7({FUNC_SET_F(foo,double)}, {min, max, step});
 *   double val = look(0.87354);
 * }
 * \endcode
 *
 *
 * \note
 * - This is currently the only LookupTableImplementation using the
 *   special_points field in LookupTableParameters. special_points is a vector of 3-tuples:
 *   {(x_1,s_1,f^{(s_1)}(x_1)),...,(x_n,s_n,f^{(s_n)}(x_n))}.
 *   n Cheby nodes are replaced with the nearest nodes {x_k} in this list, and
 *   f (or its derivate) is exact at those nodes. Doing this can drastically
 *   reduce the relative error in f and/or its derivatives. Surely the error in
 *   the resulting LUT is somehow related to the Chebyshev polynomial of the 1st kind.
 * - ChebyTable only works if we can cast both TOUT and TIN to double. This requirement
 *   exists because Armadillo Mat<T>'s `is_supported_elem_type<T>` will only let us do arithmetic
 *   with float or double (not even long double!). You might think "generic types
 *   is what arma::field" is made for but that class does nothing.
 * - the template implementation is only registered for N=1,2,3,4,5,6,7
 *   but users can manually construct this class with larger N if they wish (but we make no promises on convergence/error).
 */
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
    if(fun == nullptr)
      throw std::invalid_argument("Error in func::ChebyInterpTable: Given an invalid FunctionContainer");

    // build the default Vandermonde matrix (make this a function?)
    arma::mat Van = arma::ones(N+1,N+1);
    Van.col(1) = (1 + arma::cos( arma::datum::pi*(2*arma::linspace(1,N+1,N+1)-1) / (2*(N+1)) ))/2.0;
    for(unsigned int i=2; i<N+1; i++)
      Van.col(i) = Van.col(i-1) % Van.col(1); // % does elementwise multiplication

    /* special_points is a vector of 3-tuples:
     * {(x_1,s_1,f^{(s_1)}(x_1)),...,(x_n,s_n,f^{(s_n)}(x_n))}. Organize points in a dictionary
     * where each key (unsigned int) is the subinterval that point belongs to */
    std::map<unsigned,std::vector<std::tuple<TIN,unsigned,TOUT>>> d_intervals;
    for (std::tuple<TIN,unsigned,TOUT> tup : par.special_points){
      if(std::get<1>(tup) > N){
        std::cerr << "Warning: Given f^{(" << std::get<1>(tup) << ")}(" << std::get<0>(tup) <<
          ") = " << std::get<2>(tup) << " but a ChebyTable<" << N << ">" <<
          " can only accommodate derivatives of order at most " << N+1 << ".";
        continue;
      }
      unsigned int x0 = MetaTable<N+1,TIN,TOUT,GT>::template hash<GT>(std::get<0>(tup)).first;
      d_intervals[x0].emplace_back(tup);
    }

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

      // Chebyshev nodes over [x,x+h]:
      arma::vec xvec = x + h*(1 + arma::cos( arma::datum::pi*(2*arma::linspace(1,N+1,N+1)-1)/(2*(N+1)) ))/2.0;
      // y will contain our desired coefficients after the following if statement
      arma::vec y(N+1);

      /* check if this subinterval contain any points the user wants _exact_ */
      auto iter = d_intervals.find(ii);
      if(iter == d_intervals.end()){
        for (unsigned int k=0; k<N+1; k++)
          y[k] = fun(static_cast<TIN>(xvec[k]));
        y = arma::solve(Van,y); /* solve_opts strictly waste time with cheby nodes */
      }else{
        /* replace each nearest chebyshev node in this interval with the user's nodes.*/
        arma::vec  xvec2 = xvec;
        arma::Col<unsigned> svec(N+1);
        for(unsigned int k=0; k < std::min(N+1,static_cast<unsigned>(iter->second.size())); k++){
          auto tup = iter->second[k];
          TIN t_k = std::get<0>(tup);
          svec[k] = std::get<1>(tup);
          arma::uword i = arma::index_min(abs(xvec2 - t_k));
          xvec(i) = t_k;
          y[i] = std::get<2>(tup);
          xvec2[i] = std::numeric_limits<TIN>::quiet_NaN();
        }

        /* compute f(x) for every point that wasn't provided */
        for(unsigned int k=0; k<N+1; k++){
          if(!std::isnan(xvec2[k]))
            y[k] = fun(static_cast<TIN>(xvec2[k]));
        }

        arma::mat Van2 = arma::ones(N+1,N+1);
        Van2.col(1) = (xvec-x)/h;
        for(unsigned int i=2; i<N+1; i++)
          Van2.col(i) = Van2.col(i-1) % Van2.col(1);

        /* differentiate rows of Van2. This row will not be entirely 0 because N+1>s */
        for(unsigned int k=0; k<N+1; k++){
          auto s = svec[k];
          if(s == 0) continue;

          arma::vec D(N+1, arma::fill::zeros);
          for(unsigned int r=0; r<N+1-s; r++)
            D[r] = permutation(s+r,s);
          Van2.row(k) = D % Van2.row(k);
        }

        /* reorder so the data is in increasing order to improve the conditioning of V */
        arma::uvec idxvec = arma::sort_index(xvec, "ascend");
        //xvec = xvec(idxvec);
        y = y(idxvec);
        Van2 = Van2.rows(idxvec);

        /* TODO solve_opts? Any good ones for a very poorly conditioned matrix? */
        y = arma::solve(Van2,y);
      }

      for(unsigned int k=0; k<N+1; k++)
        m_table[ii].coefs[k] = y[k];

      FUNC_IF_CONSTEXPR(GT == GridTypes::NONUNIFORM){
        auto p = m_table[ii];
        for(unsigned int s=0; s<N+1; s++)
          m_table[ii].coefs[s] = polynomial_diff(p,-x/h,s)/static_cast<TIN>(pow(h,s))/static_cast<TIN>(factorial(s));
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

#pragma once
#include "MetaTable.hpp"
#include "config.hpp"
#include <stdexcept>
#include <boost/math/special_functions/binomial.hpp>

#ifdef FUNC_USE_ARMADILLO
#include <armadillo>
#endif

namespace func {

/** \brief LUT using degree 1 to 7 polynomial interpolation over Chebyshev nodes on each subinterval
 *  \ingroup MetaTable
 *
 * Chebyshev nodes are a partition of the interval \f$[a,b]\f$ such that
 *  \f[ t_s = (a+b)/2 + (b-a)\cos\left(\frac{2s-1}{2n}\pi\right)/2, \quad s=1,...,n. \f]
 *
 * Example usage
 * \code{.cpp}
 * // return x^9
 * template <typename T>
 * T foo(T x){ return (x*x*x)*(x*x*x)*(x*x*x); }
 *
 * int main(){
 *   auto min = -1.0; auto max = 1.0; auto step = 0.001;
 *   UniformChebyInterpTable<1,double> L1({FUNC_SET_F(foo,double)}, {min, max, step}); // degree 1
 *   UniformChebyInterpTable<2,double> L2({FUNC_SET_F(foo,double)}, {min, max, step}); // degree 2
 *   UniformChebyInterpTable<3,double> L3({FUNC_SET_F(foo,double)}, {min, max, step});
 *   UniformChebyInterpTable<4,double> L4({FUNC_SET_F(foo,double)}, {min, max, step});
 *   UniformChebyInterpTable<5,double> L5({FUNC_SET_F(foo,double)}, {min, max, step});
 *   UniformChebyInterpTable<6,double> L6({FUNC_SET_F(foo,double)}, {min, max, step});
 *   UniformChebyInterpTable<7,double> L7({FUNC_SET_F(foo,double)}, {min, max, step});
 *   // similar for nonuniform case:
 *   NonUniformChebyInterpTable<1,double> L7({FUNC_SET_F(foo,double)}, {min, max, step}); // deg 1 nonuniform
 *   NonUniformChebyInterpTable<4,double> L7({FUNC_SET_F(foo,double)}, {min, max, step});
 *   double val = look(0.87354);
 *   // The following preserves the root of $f$ at $x=0$.
 *   UniformChebyInterpTable<2,double> L2({FUNC_SET_F(foo,double)}, {min, max, step, {{0, 0, 0.0}}};
 * }
 * \endcode
 *
 * \note The template implementation is only registered in the factory for
 *   \f$N=1,2,3,4,5,6,7\f$ but users could construct this class with larger
 *   \f$N\f$ if they wish. We make no promises on convergence/error in this
 *   case.
 * \note ChebyTable only works if we can cast both TOUT and TIN to double. This requirement
 *   exists because Armadillo `Mat<T>`'s `is_supported_elem_type<T>` will only let us do arithmetic
 *   with float or double (not even long double!). You might think "generic types
 *   is what arma::field is made for" but that class appears to do nothing.
 * \note This is currently the only LookupTableImplementation using the
 *   special_points field in LookupTableParameters -- a vector of 3-tuples:
 *   \f[{(x_1,s_1,f^{(s_1)}(x_1)),...,(x_n,s_n,f^{(s_n)}(x_n))}.\f]
 *   Then, Chebyshev nodes are replaced with the nearest nodes \f${x_k}\f$ in this list, and
 *   \f$f\f$ (or its derivate) is exact at those nodes. Doing this can
 *   reduce the relative error in \f$f\f$ and/or its derivatives.
 * \warning Attempting to make \f$f\f$'s derivatives exact at some given
 *   special points can result in a singular Vandermonde matrix. We don't have
 *   a way to handle this issue at the moment. 
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
      double x;
      double h = static_cast<double>(m_stepSize);
      // (possibly) transform the uniform grid into a nonuniform grid
      FUNC_IF_CONSTEXPR(GT == GridTypes::UNIFORM)
        x = static_cast<double>(m_minArg + ii*m_stepSize);
      else{
        x = static_cast<double>(m_transferFunction(m_minArg + ii*m_stepSize));
        h = static_cast<double>(m_transferFunction(m_minArg + (ii+1)*m_stepSize) - x);
      }

      // Chebyshev nodes over [x,x+h]:
      arma::vec xvec = x + h*(1.0 + arma::cos( arma::datum::pi*(2*arma::linspace(1,N+1,N+1)-1)/(2*(N+1)) ))/2.0;
      // y will contain our desired coefficients after the following if statement
      arma::vec y(N+1);

      /* check if this subinterval contain any points the user wants _exact_ */
      auto iter = d_intervals.find(ii);
      if(iter == d_intervals.end()){
        for (unsigned int k=0; k<N+1; k++)
          y[k] = static_cast<double>(fun(static_cast<TIN>(xvec[k])));
        y = arma::solve(Van,y); /* solve_opts strictly waste time with cheby nodes */
      }else{
        /* replace each nearest chebyshev node in this interval with the user's nodes.*/
        arma::vec  xvec2 = xvec;
        arma::Col<unsigned> svec(N+1);
        for(unsigned int k=0; k < std::min(N+1,static_cast<unsigned>(iter->second.size())); k++){
          auto tup = iter->second[k];
          double t_k = static_cast<double>(std::get<0>(tup));
          svec[k] = std::get<1>(tup);
          arma::uword i = arma::index_min(abs(xvec2 - t_k));
          xvec(i) = t_k;
          y[i] = static_cast<double>(std::get<2>(tup));
          xvec2[i] = std::numeric_limits<double>::quiet_NaN();
        }

        /* compute f(x) for every point that wasn't provided */
        for(unsigned int k=0; k<N+1; k++){
          if(!std::isnan(xvec2[k]))
            y[k] = static_cast<double>(fun(static_cast<TIN>(xvec2[k])));
        }

        arma::mat Van2 = arma::ones(N+1,N+1);
        Van2.col(1) = (xvec-x)/h;
        for(unsigned int i=2; i<N+1; i++)
          Van2.col(i) = Van2.col(i-1) % Van2.col(1);

        /* differentiate rows of Van2. This row will not be entirely 0 because N+1 > s */
        for(unsigned int k=0; k<N+1; k++){
          auto s = svec[k];
          if(s == 0) continue;

          arma::rowvec D(N+1, arma::fill::zeros);
          for(unsigned int r=0; r<N+1-s; r++)
            D[s+r] = permutation(r+s,s)*Van2(k,r);
            //D[s+r] = permutation(1+s+r,1+s)*Van2(k,r);
          Van2.row(k) = D;
        }

        /* reorder so the data is in increasing order to improve the conditioning of V */
        arma::uvec idxvec = arma::sort_index(xvec, "ascend");
        y = y(idxvec);
        Van2 = Van2.rows(idxvec);

        //std::cout << Van2 << " " << y << std::endl;

        /* TODO solve_opts? Any good ones for a probably poorly conditioned matrix? */
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
    /* special case to make lut(tableMaxArg) work. Move the second last polynomial into the last interval (shifting the domain for uniform LUTs) */
    FUNC_IF_CONSTEXPR(GT == GridTypes::UNIFORM)
      m_table[m_numTableEntries-1] = taylor_shift(m_table[m_numTableEntries-2], static_cast<TIN>(1), static_cast<TIN>(2), static_cast<TIN>(0), static_cast<TIN>(1));
    else
      m_table[m_numTableEntries-1] = m_table[m_numTableEntries-2];
#endif
  }
};

// define friendlier names
template <unsigned int N, typename TIN, typename TOUT=TIN>
using UniformChebyInterpTable = ChebyInterpTable<N,TIN,TOUT,GridTypes::UNIFORM>;
template <unsigned int N, typename TIN, typename TOUT=TIN>
using NonUniformChebyInterpTable = ChebyInterpTable<N,TIN,TOUT,GridTypes::NONUNIFORM>;

} // namespace func

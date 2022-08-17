/*
  4th to 7th degree polynomial interpolation LUT with uniform sampling (precomputed
  coefficients solved using Armadillo matrices)

  Usage example for a 4th degree interpolant:
    ArmadilloPrecomputedInterpolationTable<double,double,4> look(&function,LookupTableParameters<double>{0,10,0.0001});
    double val = look(0.87354);

  Notes:
  - This class only works if TOUT and TIN can both be cast to double. 
    Armadillo Mat<T>'s `is_supported_elem_type<T>` will only let us do arithmetic
    with float or double (not even long double!) and arma::field is useless.
    TODO disable constructor if TOUT cannot be cast to double

  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
  - the template implementation is only registered for N=4,5,6,7
  (ie, available polynomial interpolation is of degrees 4 up to degree 7)
  but users can manually construct this class with larger N
*/
#pragma once
#include "MetaTable.hpp"
#include "config.hpp"
#include <stdexcept>

#ifdef FUNC_USE_ARMADILLO
#include <armadillo>
#endif

// TODO what default lets us higher precision types to their full extent?????
#define FUNC_ARMA_LU_SOLVE
#ifdef FUNC_ARMA_LU_SOLVE
#define FUNC_ARMA_SOLVE_OPTS arma::solve_opts::none
#else
//#define FUNC_ARMA_SOLVE_OPTS arma::solve_opts::none
#define FUNC_ARMA_SOLVE_OPTS arma::solve_opts::refine
#endif

namespace func {

template <typename TIN, typename TOUT, unsigned int N, GridTypes GT=GridTypes::UNIFORM>
class ArmadilloPrecomputedInterpolationTable final : public MetaTable<TIN,TOUT,N+1,GT>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,N+1,GT);

  static const std::string classname;
public:
  // build the LUT from scratch or look in filename for an existing LUT
  ArmadilloPrecomputedInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, LookupTableParameters<TIN> par,
      const nlohmann::json& jsonStats=nlohmann::json()) :
    MetaTable<TIN,TOUT,N+1,GT>(jsonStats.empty() ? // use the default move constructor for MetaTable (probably not elided...)
      std::move(MetaTable<TIN,TOUT,N+1,GT>(func_container, par)) :
      std::move(MetaTable<TIN,TOUT,N+1,GT>(jsonStats, classname, func_container)))
  {
#ifndef FUNC_USE_ARMADILLO
    // throw a descriptive exception. This could theoretically be a compile time error; however, that will only stop us from
    // registering this table (in which case, users would just get a vague "table not registered" error)
    if(jsonStats.empty())
      throw std::invalid_argument("Armadillo tables need Armadillo to be generated");
#else
    if(!jsonStats.empty())
      return; // all our work is already done

    /* Base class default variables */
    m_name = classname;
    m_numTableEntries = m_numIntervals+1;
    m_order = N+1; // N is the degree of the polynomial interpolant so the order is N+1
    m_dataSize = static_cast<unsigned>(sizeof(m_table[0]) * (m_numTableEntries));

    /* build the vandermonde system for finding the interpolating polynomial's coefficients */
    arma::mat Van = arma::ones(N+1, N+1);
    Van.col(1) = arma::linspace(0,1,N+1);
    for(unsigned int i=2; i<N+1; i++)
      Van.col(i) = Van.col(i-1) % Van.col(1); // the % does elementwise multiplication

#ifdef FUNC_ARMA_LU_SOLVE
    // LU factor the matrix we just built
    arma::mat L, U, P;
    arma::lu(L,U,P,Van);
#endif

    /* Allocate and set table */
    m_grid.reset(new TIN[m_numTableEntries]);
    m_table.reset(new polynomial<TOUT,N+1>[m_numTableEntries]);
    for (unsigned int ii=0;ii<m_numTableEntries-1;++ii) {
      TIN x;
      TIN h = m_stepSize;
      // (possibly) transform the uniform grid into a nonuniform grid
      if (GT == GridTypes::UNIFORM)
        x = m_minArg + ii*m_stepSize;
      else{
        x = m_transferFunction.g(m_minArg + ii*m_stepSize);
        h = m_transferFunction.g(m_minArg + (ii+1)*m_stepSize) - x;
      }
      // grid points
      m_grid[ii] = x;
      // build the vector of coefficients from function values
      arma::vec xvec = arma::linspace(static_cast<double>(x),static_cast<double>(x+h),N+1);
      arma::vec y(N+1);
      // TODO is this performed in simd: y.for_each([this](Mat<double>::elem_type& xk) { xk = static_cast<double>(m_func(static_cast<TIN>(xk))); });
      for (unsigned int k=0; k<N+1; k++)
        y[k] = static_cast<double>(m_func(static_cast<TIN>(xvec[k])));

      // make y the coefficients of the polynomial interpolant
#ifdef FUNC_ARMA_LU_SOLVE
      y = arma::solve(arma::trimatu(U), arma::solve(arma::trimatl(L), P*y));
#else
      y = arma::solve(Van, y, FUNC_ARMA_SOLVE_OPTS);
#endif

      // move this back into the m_table array
      for (unsigned int k=0; k<N+1; k++)
        m_table[ii].coefs[k] = static_cast<TOUT>(y[k]);
    }
    // special case to make lut(tableMaxArg) work
    m_grid[m_numTableEntries-1] = m_tableMaxArg;
    m_table[m_numTableEntries-1].coefs[0] = m_func(m_tableMaxArg);
    for (unsigned int k=1; k<N+1; k++)
      m_table[m_numTableEntries-1].coefs[k] = 0;
#endif
  }

  // operator() comes straight from the MetaTable
};

template <typename TIN, typename TOUT, unsigned int N, GridTypes GT>
const std::string ArmadilloPrecomputedInterpolationTable<TIN,TOUT,N,GT>::classname = grid_type_to_string<GT>() + "ArmadilloPrecomputedInterpolationTable<" + std::to_string(N) + ">";

// define friendlier names
template <typename TIN, typename TOUT, unsigned int N>
using UniformArmadilloPrecomputedInterpolationTable = ArmadilloPrecomputedInterpolationTable<TIN,TOUT,N,GridTypes::UNIFORM>;
template <typename TIN, typename TOUT, unsigned int N>
using NonUniformArmadilloPrecomputedInterpolationTable = ArmadilloPrecomputedInterpolationTable<TIN,TOUT,N,GridTypes::NONUNIFORM>;
template <typename TIN, typename TOUT, unsigned int N>
using NonUniformPseudoArmadilloPrecomputedInterpolationTable = ArmadilloPrecomputedInterpolationTable<TIN,TOUT,N,GridTypes::NONUNIFORM_PSEUDO>;

} // namespace func

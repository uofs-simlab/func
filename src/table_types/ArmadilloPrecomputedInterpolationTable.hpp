/*
  4th to 7th degree polynomial interpolation LUT with uniform sampling (precomputed
  coefficients solved using Armadillo matrices)

  Usage example for a 4th degree interpolant:
    ArmadilloPrecomputedInterpolationTable<4> look(&function,0,10,0.0001);
    double val = look(0.87354);

  Notes:
  - table precomputes and stores the linear coefficient so it doesn't have to
    perform that operation every lookup (but does have to look it up)
  - static data after constructor has been called
  - evaluate by using parentheses, just like a function
  - the template implementation is only for N=4,5,6,7 (ie, available polynomial interpolation is
  of degrees 4 up to degree 7)
*/
#pragma once
#include "MetaTable.hpp"
#include "config.hpp"
#include <stdexcept>

#ifdef FUNC_USE_ARMADILLO
#define ARMA_USE_CXX11 // TODO does this affect a libraries that use FunC?
#include <armadillo>
#endif

// TODO there is some error introduced by factoring the
// vandermonde matrix before doing an LU solve.
// Is this added error worth the increase in speed?
// Should we do iterative refinement?
#define FUNC_ARMA_LU_SOLVE
#ifdef FUNC_ARMA_LU_SOLVE
#define FUNC_ARMA_SOLVE_OPTS arma::solve_opts::none
#else
//#define FUNC_ARMA_SOLVE_OPTS arma::solve_opts::none
#define FUNC_ARMA_SOLVE_OPTS arma::solve_opts::refine
#endif

template <typename TIN, typename TOUT, unsigned int N, GridTypes GT=UNIFORM>
class ArmadilloPrecomputedInterpolationTable final : public MetaTable<TIN,TOUT,N+1,GT>
{
  INHERIT_EVALUATION_IMPL(TIN,TOUT);
  INHERIT_LUT(TIN,TOUT);
  INHERIT_META(TIN,TOUT,N+1,GT);

public:
  // build the LUT from scratch or look in filename for an existing LUT
  ArmadilloPrecomputedInterpolationTable(FunctionContainer<TIN,TOUT> *func_container, LookupTableParameters<TIN> par, std::string filename = "") :
    MetaTable<TIN,TOUT,N+1,GT>((filename == "") ? // use the default move constructor for MetaTable (probably not elided...)
      std::move(MetaTable<TIN,TOUT,N+1,GT>(func_container, par)) :
      std::move(MetaTable<TIN,TOUT,N+1,GT>(func_container, filename, grid_type_to_string<GT>() + "ArmadilloPrecomputedInterpolationTable<" + std::to_string(N) + ">")))
  {
#ifndef FUNC_USE_ARMADILLO
    // throw a descriptive exception. This could theoretically be a compile time error; however, that will only stop us from
    // registering this table (in which case, users would just get a vague "table not registered" error)
    if(filename == "")
      throw std::invalid_argument("Armadillo tables need Armadillo to be generated");
#else
    /* Base class default variables */
    m_name = grid_type_to_string<GT>() + "ArmadilloPrecomputedInterpolationTable<" + std::to_string(N) + ">";
    m_numTableEntries = m_numIntervals+1;
    m_order = N+1; // N is the degree of the polynomial interpolant so the order is N+1
    m_dataSize = (unsigned) sizeof(m_table[0]) * (m_numTableEntries);

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
    m_table.reset(new polynomial<TOUT,N+1>[m_numTableEntries]);
    for (unsigned int ii=0;ii<m_numIntervals;++ii) {
      TIN x;
      TIN h = m_stepSize;
      // (possibly) transform the uniform grid into a nonuniform grid
      if (GT == UNIFORM)
        x = m_minArg + ii*m_stepSize;
      else{
        x = m_transferFunction.g(m_minArg + ii*m_stepSize);
        h = m_transferFunction.g(m_minArg + (ii+1)*m_stepSize) - x;
      }
      // grid points
      m_grid[ii] = x;
      // build the vector of coefficients from function values
      arma::vec y = arma::linspace(x,x+h,N+1);
      for (unsigned int k=0; k<N+1; k++)
        y[k] = m_func(y[k]);

      // make y the coefficients of the polynomial interpolant
#ifdef FUNC_ARMA_LU_SOLVE
      y = arma::solve(trimatu(U), arma::solve(trimatl(L), P*y));
#else
      y = arma::solve(Van, y, FUNC_ARMA_SOLVE_OPTS);
#endif

      // move this back into the m_table array
      for (unsigned int k=0; k<N+1; k++)
        m_table[ii].coefs[k] = y[k];
    }
#endif
  }

  // Look for this LUT in filename. MetaTable builds m_table and LookupTable initializes every other member variable
  // I think this can probably just be a from_json() method ...
  //ArmadilloPrecomputedInterpolationTable(std::string filename, FunctionContainer<TIN,TOUT> *func_container = nullptr) :
  //  ArmadilloPrecomputedInterpolationTable(func_container, {0,0,0}, filename)

  // operator() comes straight from the MetaTable
};

// define friendlier names
template <typename TIN, typename TOUT, unsigned int N>
using UniformArmadilloPrecomputedInterpolationTable = ArmadilloPrecomputedInterpolationTable<TIN,TOUT,N,UNIFORM>;
template <typename TIN, typename TOUT, unsigned int N>
using NonUniformArmadilloPrecomputedInterpolationTable = ArmadilloPrecomputedInterpolationTable<TIN,TOUT,N,NONUNIFORM>;
template <typename TIN, typename TOUT, unsigned int N>
using NonUniformPseudoArmadilloPrecomputedInterpolationTable = ArmadilloPrecomputedInterpolationTable<TIN,TOUT,N,NONUNIFORM_PSEUDO>;


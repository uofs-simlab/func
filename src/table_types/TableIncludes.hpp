/*
  Header file for every table type
*/
#pragma once

/* ---- Uniform Tables ---- */
#include "UniformLookupTable.hpp"

// differentiation tables
#ifdef FUNC_USE_BOOST_AUTODIFF
  #include "UniformLinearTaylorTable.hpp"
  #include "UniformQuadraticTaylorTable.hpp"
  #include "UniformCubicTaylorTable.hpp"
  #include "UniformCubicHermiteTable.hpp"

  #ifdef FUNC_USE_ARMADILLO
    #include "UniformPadeTable.hpp"
  #endif
#endif

// standard interpolation tables
#include "UniformConstantTaylorTable.hpp"
#include "UniformLinearInterpolationTable.hpp"
#include "UniformLinearPrecomputedInterpolationTable.hpp"
#include "UniformQuadraticPrecomputedInterpolationTable.hpp"
#include "UniformCubicPrecomputedInterpolationTable.hpp"
#ifdef FUNC_USE_ARMADILLO
  #include "UniformArmadilloPrecomputedInterpolationTable.hpp"
#endif


/* ---- NonUniform Tables ---- */
// Currently an armadillo and boost extension
//#ifdef FUNC_USE_ARMADILLO
//#ifdef FUNC_USE_BOOST_AUTODIFF
//  #include "NonUniformLookupTable.hpp"
//  #include "NonUniformLinearInterpolationTable.hpp"
//  #include "NonUniformPseudoLinearInterpolationTable.hpp"
//  #include "NonUniformCubicPrecomputedInterpolationTable.hpp"
//  #include "NonUniformPseudoCubicPrecomputedInterpolationTable.hpp"
//#endif
//#endif


/* ---- Table Containers ---- */
#include "FailureProofTable.hpp"
#include "CompositeLookupTable.hpp"

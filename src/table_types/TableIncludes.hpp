/*
  Header file for every table type
*/
#pragma once

/* ---- Uniform Tables ---- */
#include "UniformLookupTable.hpp"

// differentiation tables
#include "UniformConstantTaylorTable.hpp"
#include "UniformLinearTaylorTable.hpp"
#include "UniformQuadraticTaylorTable.hpp"
#include "UniformCubicTaylorTable.hpp"
#include "UniformCubicHermiteTable.hpp"
#ifdef USE_ARMADILLO
  #include "UniformPadeTable.hpp"
#endif

// standard interpolation tables
#include "UniformLinearInterpolationTable.hpp"
#include "UniformLinearPrecomputedInterpolationTable.hpp"
#include "UniformCubicPrecomputedInterpolationTable.hpp"
#include "UniformQuadraticPrecomputedInterpolationTable.hpp"
#ifdef USE_ARMADILLO
  #include "UniformArmadilloPrecomputedInterpolationTable.hpp"
#endif


/* ---- NonUniform Tables ---- */
// Currently an armadillo extension
#ifdef USE_ARMADILLO
  #include "NonUniformLookupTable.hpp"
  #include "NonUniformLinearInterpolationTable.hpp"
  #include "NonUniformPseudoLinearInterpolationTable.hpp"
  #include "NonUniformCubicPrecomputedInterpolationTable.hpp"
  #include "NonUniformPseudoCubicPrecomputedInterpolationTable.hpp"
#endif


/* ---- Meta Tables ---- */
#include "FailureProofTable.hpp"
#include "CompositeLookupTable.hpp"

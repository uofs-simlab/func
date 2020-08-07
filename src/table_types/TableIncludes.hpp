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
#include "NonUniformLookupTable.hpp"
#include "NonUniformLinearInterpolationTable.hpp"
#include "NonUniformCubicPrecomputedInterpolationTable.hpp"


/* ---- Meta Tables ---- */
#include "FailureProofTable.hpp"
#include "CompositeLookupTable.hpp"

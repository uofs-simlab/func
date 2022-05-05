/*
  Header file for every table type
*/
#pragma once

/* ---- Uniform Tables ---- */
#include "LookupTable.hpp"
#include "MetaTable.hpp"

// differentiation tables
#ifdef FUNC_USE_BOOST_AUTODIFF
  #include "LinearTaylorTable.hpp"
  #include "QuadraticTaylorTable.hpp"
  #include "CubicTaylorTable.hpp"
  #include "CubicHermiteTable.hpp"

  #ifdef FUNC_USE_ARMADILLO
    #include "PadeTable.hpp"
  #endif
#endif

// standard interpolation tables
#include "ConstantTaylorTable.hpp"
#include "LinearInterpolationTable.hpp"
#include "LinearPrecomputedInterpolationTable.hpp"
#include "QuadraticPrecomputedInterpolationTable.hpp"
#include "CubicPrecomputedInterpolationTable.hpp"
#ifdef FUNC_USE_ARMADILLO
  #include "ArmadilloPrecomputedInterpolationTable.hpp"
#endif

/* ---- Table Containers ---- */
#include "FailureProofTable.hpp"
#include "CompositeLookupTable.hpp"

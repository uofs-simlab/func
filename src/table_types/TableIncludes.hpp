/*
  Header file for every table type
*/
#pragma once

/* ---- Parent Tables ---- */
#include "LookupTable.hpp"
#include "MetaTable.hpp"

// standard interpolation tables
#include "ConstantTaylorTable.hpp"
#include "LinearInterpolationTable.hpp"
#include "LinearPrecomputedInterpolationTable.hpp"
#include "QuadraticPrecomputedInterpolationTable.hpp"
#include "CubicPrecomputedInterpolationTable.hpp"
#include "ArmadilloPrecomputedInterpolationTable.hpp"

// tables using derivatives
#include "LinearTaylorTable.hpp"
#include "QuadraticTaylorTable.hpp"
#include "CubicTaylorTable.hpp"
#include "CubicHermiteTable.hpp"
#include "PadeTable.hpp"

/* ---- Table Containers ---- */
#include "FailureProofTable.hpp"
#include "CompositeLookupTable.hpp"

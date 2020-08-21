// Sadly, c++17 is the first standard where we can use inline variables.
// So, the registrar needs all these static const definitions to be out of line.
#include "TableIncludes.hpp"
#include "config.hpp" // FUNC_USE_BOOST_AUTODIFF

#include "UniformConstantTaylorTable.hpp"
#include "UniformLinearInterpolationTable.hpp"
#include "UniformLinearPrecomputedInterpolationTable.hpp"
#include "UniformQuadraticPrecomputedInterpolationTable.hpp"
#include "UniformCubicPrecomputedInterpolationTable.hpp"

FUNC_REGISTER_DOUBLE_AND_FLOAT_LUT_IMPLS(UniformLinearInterpolationTable);
FUNC_REGISTER_DOUBLE_AND_FLOAT_LUT_IMPLS(UniformLinearPrecomputedInterpolationTable);
FUNC_REGISTER_DOUBLE_AND_FLOAT_LUT_IMPLS(UniformQuadraticPrecomputedInterpolationTable);
FUNC_REGISTER_DOUBLE_AND_FLOAT_LUT_IMPLS(UniformCubicPrecomputedInterpolationTable);
FUNC_REGISTER_DOUBLE_AND_FLOAT_LUT_IMPLS(UniformConstantTaylorTable);

#ifdef FUNC_USE_BOOST_AUTODIFF
#include "UniformLinearTaylorTable.hpp"
#include "UniformQuadraticTaylorTable.hpp"
#include "UniformCubicTaylorTable.hpp"
#include "UniformCubicHermiteTable.hpp"

FUNC_REGISTER_DOUBLE_AND_FLOAT_LUT_IMPLS(UniformLinearTaylorTable);
FUNC_REGISTER_DOUBLE_AND_FLOAT_LUT_IMPLS(UniformQuadraticTaylorTable);
FUNC_REGISTER_DOUBLE_AND_FLOAT_LUT_IMPLS(UniformCubicTaylorTable);
FUNC_REGISTER_DOUBLE_AND_FLOAT_LUT_IMPLS(UniformCubicHermiteTable);
#endif

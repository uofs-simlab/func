// Sadly, c++17 is the first standard where we can use inline variables.
// So, the registrar needs all these static const definitions to be out of line.
#include "TableIncludes.hpp"
#include "config.hpp" // FUNC_USE_BOOST_AUTODIFF

#include "UniformConstantTaylorTable.hpp"
#include "UniformLinearInterpolationTable.hpp"
#include "UniformLinearPrecomputedInterpolationTable.hpp"
#include "UniformQuadraticPrecomputedInterpolationTable.hpp"
#include "UniformCubicPrecomputedInterpolationTable.hpp"

FUNC_REGISTER_EACH_ULUT_IMPL(UniformLinearInterpolationTable);
FUNC_REGISTER_EACH_ULUT_IMPL(NonUniformLinearInterpolationTable);
FUNC_REGISTER_EACH_ULUT_IMPL(NonUniformPseudoLinearInterpolationTable);
FUNC_REGISTER_EACH_ULUT_IMPL(UniformLinearPrecomputedInterpolationTable);
FUNC_REGISTER_EACH_ULUT_IMPL(UniformQuadraticPrecomputedInterpolationTable);
FUNC_REGISTER_EACH_ULUT_IMPL(UniformCubicPrecomputedInterpolationTable);
FUNC_REGISTER_EACH_ULUT_IMPL(UniformConstantTaylorTable);

#ifdef FUNC_USE_BOOST_AUTODIFF
#include "UniformLinearTaylorTable.hpp"
#include "UniformQuadraticTaylorTable.hpp"
#include "UniformCubicTaylorTable.hpp"
#include "UniformCubicHermiteTable.hpp"

FUNC_REGISTER_EACH_ULUT_IMPL(UniformLinearTaylorTable);
FUNC_REGISTER_EACH_ULUT_IMPL(UniformQuadraticTaylorTable);
FUNC_REGISTER_EACH_ULUT_IMPL(UniformCubicTaylorTable);
FUNC_REGISTER_EACH_ULUT_IMPL(UniformCubicHermiteTable);
#endif

/* TODO Separate every LUT definition/template instantiation into its
 * own translation units so compilation will be massively parallel */

#include "ConstantTaylorTable.hpp"
#include "LinearInterpolationTable.hpp"
#include "LinearPrecomputedInterpolationTable.hpp"
#include "QuadraticPrecomputedInterpolationTable.hpp"
#include "CubicPrecomputedInterpolationTable.hpp"
#include "TaylorTable.hpp"
#include "CubicHermiteTable.hpp"

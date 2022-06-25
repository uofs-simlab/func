// Sadly, c++17 is the first standard where we can use inline variables.
// So, the registrar needs all these static const definitions to be out of line.
#include "TableIncludes.hpp"

#include "ConstantTaylorTable.hpp"
#include "LinearInterpolationTable.hpp"
#include "LinearPrecomputedInterpolationTable.hpp"
#include "QuadraticPrecomputedInterpolationTable.hpp"
#include "CubicPrecomputedInterpolationTable.hpp"

#include "TaylorTable.hpp"
#include "CubicHermiteTable.hpp"

// Global include for FunC
#pragma once
#include "config.hpp"
#include "DirectEvaluation.hpp"
#include "TableIncludes.hpp"

#include "UniformLookupTableGenerator.hpp"
#include "ImplementationComparator.hpp"
#include "RngInterface.hpp"
#include "StdRng.hpp"
#include "TransferFunctionInterface.hpp"

#ifdef FUNC_USE_ARMADILLO
#ifdef FUNC_USE_BOOST_AUTODIFF
#include "TransferFunctionSinh.hpp"
#endif
#endif

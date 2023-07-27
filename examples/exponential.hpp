#pragma once
#include <cmath>
#define FUNC(x) (exp(x))
#define FUNCNAME "exp(x)"
#define MIN_ARG -3
#define MAX_ARG 3

template<typename T>
T MyFunction(T x) { return FUNC(x); }

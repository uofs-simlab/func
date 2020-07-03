#pragma once
#include <cmath>
#define FUNC(x) (7.7-13.0287*log(x))
#define FUNCNAME "(7.7-13.0287*log(x))"

template <typename T>
T MyFunction(T x) { return FUNC(x); }

#pragma once
#include <cmath>
// remember that cmath's log is base e
#define FUNC(x) (7.7-13.0287*log(x))
#define FUNCNAME "(7.7-13.0287*log(x))"

template <typename T>
T MyFunction(T x) { return FUNC(x); }

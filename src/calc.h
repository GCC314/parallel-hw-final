#pragma once

#include "interpolate.h"
#include "input.h"
#include <vector>

void calculator(double* H, int N,
    Interpolator& itp, v_data &V, std::vector<point_data> &points,
    double dx, double dy, double dz);
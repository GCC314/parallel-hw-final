#pragma once

#include <cstring>

const int D_MODE_LAP = 1, D_MODE_SCA = 2;

void diagonize(int rank, int size, double *H, double *W, int N, 
    int mode);
#pragma once

#include "mpi.h"
#include <cstring>
#include <iostream>
using std::string;

void myabort(const string &info){
    std::cerr << "Error: " << info << "\n";
    #ifdef __MPI
        MPI_Abort(MPI_COMM_WORLD, 1);
    #else
        abort();
    #endif
}

void myassert(bool expr){
    if(!expr) myabort("my assertion failed!");
}
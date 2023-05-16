#include "mytimer.h"
#include "mpi.h"
#include <map>
#include <cstring>
#include <iostream>

std::map<std::string, double> mytimer::ts;

void mytimer::start(const std::string& name){
    double now = MPI_Wtime();
    ts[name] = now;
}

void mytimer::end(const std::string& name){
    double now = MPI_Wtime();
    std::cout << "[timer] " << name << " : ";
    std::cout << now - ts[name] << "s\n";
    ts[name] = now;
}
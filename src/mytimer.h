#ifndef __MYTIMER_H__
#define __MYTIMER_H__

#include "mpi.h"
#include <map>
#include <cstring>
#include <iostream>

class mytimer{
    public:
    static std::map<std::string, double> ts;
    static void start(const std::string& name);
    static void end(const std::string& name);
};

#endif
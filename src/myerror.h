#ifndef __MYERROR_H__
#define __MYERROR_H__

#include "mpi.h"
#include <cstring>
#include <iostream>
using std::string;

void myabort(const string &info);

void myassert(bool expr);

#endif
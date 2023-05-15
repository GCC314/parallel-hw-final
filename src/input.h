#pragma once

#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include <cctype>

#include "myerror.h"

using std::string;

typedef std::unordered_map<string, string> input_args;

struct v_data{
    int nx, ny, nz;
    double *d;
    ~v_data(){free(d);}
};

struct dist_data{
    double cutoff;
    double dr;
    int mesh;
    double *f;
    ~dist_data(){free(f);}
};

struct point_data{
    double x, y, z;
};


string ToLower(const string &src){
    string ret = "";
    for(auto c : src) ret += tolower(c);
    return ret;
}

string getarg(const input_args &M, const string &key){
    if(M.find(key) == M.end()) myabort("Missing key: " + key);
    return M[key];
}


input_args parseInput(const string &fname);
v_data parseV(const string &fname);
dist_data parseDist(const string &fname);
std::vector<point_data> parsePoints(const string &fname);


void checkInput(const input_args &args);
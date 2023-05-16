#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include <cctype>

#include "myerror.h"
#include "input.h"

using std::string;

string ToLower(const string &src){
    string ret = "";
    for(auto c : src) ret += tolower(c);
    return ret;
}

string getarg(const input_args &M, string key){
    key = ToLower(key);
    if(M.find(key) == M.end()) myabort("Missing key: " + key);
    return M.at(key);
}

bool getd(std::ifstream &ifs, double &d){
    char c;
    for(ifs >> c; \
        !ifs.eof() && !isdigit(c) && c != '.' && c != '-'; \
        ifs >> c);
    if(ifs.eof()) return false;
    ifs.putback(c);
    ifs >> d;
    return true;
}

input_args parseInput(const string &fname){
    std::ifstream ifs(fname);
    if(!ifs.is_open()) myabort("Failed to open file: " + fname);
    input_args args;
    for(string line = "";std::getline(ifs, line);){
        int ptr = line.find('#');    // find possible '#'
        if(ptr != string::npos) line = line.substr(0, ptr); // truncate
        
        std::istringstream ss(line);    // using sstream for param parsing
        string argname, argval;
        ss >> argname >> argval;
        args[ToLower(argname)] = argval;
    }
    return args;
}

template<typename F>
void parseWithFunc(const string &fname, const string &pname, F func){
    std::ifstream ifs(fname);
    if(!ifs.is_open()) myabort("Failed to open file: " + fname);
    input_args args;
    for(string line = "";std::getline(ifs, line);){
        int ptr = line.find('#');    // find possible '#'
        if(ptr != string::npos) line = line.substr(0, ptr); // truncate

        std::istringstream ss(line);    // using sstream for param parsing
        string argname, argval;
        ss >> argname;

        if(argname == pname) func(args, ifs); // call callback function
        else{
            ss >> argval;
            args[ToLower(argname)] = argval;
        }
    }
}

v_data parseV(const string &fname){
    v_data ret;
    parseWithFunc(fname, "V:", [&](const input_args& args, std::ifstream &ifs){
        ret.nx = std::stoi(getarg(args, "nx"));
        ret.ny = std::stoi(getarg(args, "ny"));
        ret.nz = std::stoi(getarg(args, "nz"));
        ret.d = (double*)malloc(sizeof(double) * ret.nx * ret.ny * ret.nz);
        for(int i = 0;i < ret.nx;i++){
            for(int j = 0;j < ret.ny;j++){
                for(int k = 0;k < ret.nz;k++){
                    getd(ifs, ret.d[i * ret.ny * ret.nz + j * ret.nz + k]);
                }
            }
        }
    });
    return ret;
}

dist_data parseDist(const string &fname){
    dist_data ret;
    parseWithFunc(fname, "f:", [&](const input_args& args, std::ifstream &ifs){
        ret.cutoff = std::stod(getarg(args, "cutoff"));
        ret.dr = std::stod(getarg(args, "dr"));
        ret.mesh = std::stod(getarg(args, "mesh"));
        myassert(getarg(args, "l") == "1");
        ret.f = (double*)malloc(sizeof(double) * ret.mesh);
        
        for(int i = 0;i < ret.mesh;i++){
            getd(ifs, ret.f[i]);
        }
    });
    return ret;
}

std::vector<point_data> parsePoints(const string &fname){
    std::ifstream ifs(fname);
    if(!ifs.is_open()) myabort("Failed to open file: " + fname);

    std::vector<point_data> ret(0);
    
    double x, y, z;
    while(1){
        bool ok = true;
        ok &= getd(ifs, x);
        ok &= getd(ifs, y);
        ok &= getd(ifs, z);
        if(!ok) break;
        ret.push_back((point_data){x, y, z});
    }
    return ret;
}

void checkInput(const input_args &args){
    myassert(getarg(args, "isHexahedral") == "0");
    myassert(getarg(args, "thetaxy") == "0");
    myassert(getarg(args, "thetayz") == "0");
    myassert(getarg(args, "thetaxz") == "0");
    myassert(getarg(args, "support_SH") == "0");
    myassert(getarg(args, "support_Periodic_Boundary") == "0");
    myassert(getarg(args, "multi_parallel_strategies") == "0");
}
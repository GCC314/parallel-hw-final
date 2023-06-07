#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include <cctype>
#include <cstdio>
#include <cstdlib>

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


#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>

const int FILE_BUFFER_SIZE = 128 * 1024 * 1024, OFFSET = 32;
// manual IO for a fast performance
v_data parseV(const string &fname){
    v_data ret;
    char *buffer = (char*) malloc(OFFSET + FILE_BUFFER_SIZE);
    int _; // anonymous temp var

    // open file
    int fd = open(fname.c_str(), O_RDONLY, 0);
    if(!~fd) myabort("Failed to open file: " + fname);

    // read nx, ny, nz;
    input_args args;
    auto parseAttr = [&](){
        string line = "";
        char c;
        for(_ = read(fd, &c, 1);c != '\n';_ = read(fd, &c, 1)) line += c;
        string argname, argval;
        std::istringstream ss(line);    // using sstream for param parsing
        ss >> argname >> argval;
        args[argname] = argval;
    };
    parseAttr(); parseAttr(); parseAttr();
    ret.nx = std::stoi(getarg(args, "nx"));
    ret.ny = std::stoi(getarg(args, "ny"));
    ret.nz = std::stoi(getarg(args, "nz"));
    ret.d = (double*)malloc(sizeof(double) * ret.nx * ret.ny * ret.nz);

    // skip header
    {
        char c;
        for(_ = read(fd, &c, 1);c != '\n';_ = read(fd, &c, 1));
    }


    // now read V !!
    int addr = 0;
    char *p_now = 0, *p_end = 0;
    bool reachend = false;
    auto loadbuf = [&](){
        p_now = buffer + OFFSET;
        ssize_t size = read(fd, buffer + OFFSET, FILE_BUFFER_SIZE);
        reachend |= (size != FILE_BUFFER_SIZE);
        p_end = buffer + OFFSET + size;
    };
    loadbuf();
    for(int i = 0;i < ret.nx;i++){
        for(int j = 0;j < ret.ny;j++){
            for(int k = 0;k < ret.nz;k++, addr++){
                if(!reachend && p_now + OFFSET >= p_end){
                    // copy and load
                    int now_offset = p_now + OFFSET - p_end;
                    memcpy(buffer, p_end - OFFSET, OFFSET);
                    loadbuf();
                    p_now = buffer + now_offset;
                }
                char *p_nxt = p_now;
                ret.d[addr] = std::strtod(p_now, &p_nxt);
                p_now = p_nxt;
            }
        }
    }

    free(buffer);
    close(fd); // close file
    return ret;
}
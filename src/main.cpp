#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstring>
#include <vector>
#include <map>
#include <unordered_map>

#include "myerror.h"
#include "input.h"

#include "mpi.h"
#include "omp.h"


int main(int argc, char **argv){
    if(argc != 2){
        std::cerr << "Usage: " << argv[0] << " [FILE]" << std::endl;
        return 1;
    }
    input_args inp = parseInput(argv[1]);
    checkInput(inp);
    
    auto points = parsePoints(getarg(inp, "points_path"));
    auto V = parseV(getarg(inp, "v_path"));
    auto dist = parseDist(getarg(inp, "distribution_path"));

    // interpolate

    // calc sum
    
    // diag

    return 0;
}
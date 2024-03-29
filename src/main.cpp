#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cstring>
#include <vector>
#include <map>
#include <unordered_map>

#include "myerror.h"
#include "input.h"
#include "interpolate.h"
#include "calc.h"
#include "mytimer.h"

#include "mpi.h"
#include "omp.h"
#include "diago.h"

int main(int argc, char **argv){
    if(argc != 3){
        std::cerr << "Usage: " << argv[0] << " [FILE] [THREAD_NUM]" << std::endl;
        return 1;
    }


    MPI_Init(&argc, &argv);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if(rank == 0){
        int thread_num = std::stoi(argv[2]);
        std::cout << "thread set: " << thread_num << '\n'; 
        omp_set_num_threads(thread_num);
    }
    
    int N = 0;
    double *H = nullptr, *W = nullptr;
    int diag_mode = 0;

    if(rank == 0){
        mytimer::start("input phase");
        input_args inp = parseInput(argv[1]);
        checkInput(inp);
        
        auto points = parsePoints(getarg(inp, "points_path"));
        auto V = parseV(getarg(inp, "v_path"));
        auto dist = parseDist(getarg(inp, "distribution_path"));

        std::cout << "Data files loaded\n";
        mytimer::end("input phase");

        N = points.size();
        std::cout << "points size: " << N << "\n";
        myassert(points.size() >= 2);

        auto itp = Interpolator(&dist);

        H = (double*)malloc(sizeof(double) * N * N);
        
        double dx = std::stod(getarg(inp, "lx")) / V.nx;
        double dy = std::stod(getarg(inp, "ly")) / V.ny;
        double dz = std::stod(getarg(inp, "lz")) / V.nz;
        calculator(H, N, itp, V, points, dx, dy, dz);

        V.unload(); itp.unload(); dist.unload();

        auto mode_str = getarg(inp, "diago_lib");
        if(mode_str == "lapack") diag_mode = D_MODE_LAP;
        else if(mode_str == "scalapack") diag_mode = D_MODE_SCA;
    }
    MPI_Bcast(&diag_mode, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank) H = (double*)malloc(sizeof(double) * N * N);
    W = (double*)malloc(sizeof(double) * N);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(H, N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(N <= 2 && !rank) for(int i = 0;i < N;i++) for(int j = 0;j < N;j++)
        std::cout << H[i * N + j] << "\n"; // for debug

    diagonize(rank, size, H, W, N, diag_mode);
    // using MPI parallelism

    free(H); free(W);
    MPI_Finalize();
    return 0;
}
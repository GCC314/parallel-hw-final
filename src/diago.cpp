#include "diago.h"
#include "myerror.h"
#include "mpi.h"
#include <lapacke.h>
#include <cstring>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "mytimer.h"

extern "C"{
/** ScaLAPACK and BLACS routines */
void sl_init_(int*, const int*, const int*);
void blacs_gridinfo_(const int*, int*, int*, int*, int*);
int  numroc_(const int*, const int*, const int*, const int*, const int*);
void descinit_(int*, const int*, const int*, const int*, const int*, const int*, const int*,
               const int*, const int*, int*);
void blacs_gridexit_(const int*);
void pdsyevd_(const char* JOBZ, const char* UPLO, const int* N, double* A, const int* IA, const int* JA,
              const int* DESCA, double* W, double* Z, const int* IZ, const int* JZ, const int* DESCZ,
              double* WORK, const int* LWORK, int* IWORK, const int* LIWORK, int* INFO);
}

void output(double *H, double *W, int N){
    {
        std::ofstream ofs("eigenvalues.log");
        if(!ofs.is_open()) myabort("Failed to write eigenvalues.log");
        ofs << "Eigenvalues:\n";
        for(int i = 0;i < N;i++) ofs << W[i] << '\n';
    }
    {
        std::ofstream ofs("eigenvectors.log");
        if(!ofs.is_open()) myabort("Failed to write eigenvectors.log");
        ofs << "Eigenvectors (row-wise):\n";
        for(int i = 0;i < N;i++) for(int j = 0;j < N;j++){
            ofs << std::setw(6) << H[i * N + j] << " \n"[j == N - 1];
        }
    }
}

void diagonize(int rank, int size, double *H, double *W, int N, 
    int mode){

    std::cout << "[diago] start\n";
    mytimer::start("diago");
    
    if(rank == 0){
        if(mode == D_MODE_LAP){ // LAPACK
            int res = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U',
                N, H, N, W);
            if(res) myabort("diago: diagonization failed.");
            output(H, W, N);
            std::cout << "[diago] end\n";
            mytimer::end("diago");
            return;
        }
        else if(mode != D_MODE_SCA) // UNKNOWN
            myabort("diago_lib: not supported.");
    }

    // SCALAPACK

    char jobz = 'V';
    char uplo = 'U';
    int IA = 1;         // Start position of submatrix
    int ierr = 0;
    int B = 4;
    int SRC = 0;
    int ICTXT = 0;
    int DESC[9];
    int info = 0;
    int proc_row = 0, proc_col = 0;
    double workN = 0.0;

    int num_proc_rows = (int)sqrt(double(size));
    int num_proc_cols = size / num_proc_rows;

    while(num_proc_cols * num_proc_rows != size)
    {   num_proc_rows--;
        num_proc_cols = size / num_proc_rows;
    }

    // Initialise process grid
    sl_init_(&ICTXT, &num_proc_rows, &num_proc_cols);
    blacs_gridinfo_(&ICTXT, &num_proc_rows, &num_proc_cols, &proc_row, &proc_col);

    // Get dimensions of local array
    int M_rows = numroc_(&N, &B, &proc_row, &SRC, &num_proc_rows);
    int M_cols = numroc_(&N, &B, &proc_col, &SRC, &num_proc_cols);

    // Initialise array descriptor (DESC) for dense matrix
    descinit_(DESC, &N, &N, &B, &B, &SRC, &SRC, &ICTXT, &M_rows, &info);

    // Allocate the output (eigenvector) matrix
    int DESC_V[9];
    descinit_(DESC_V, &N, &N, &B, &B, &SRC, &SRC, &ICTXT, &M_rows, &ierr);
    double* V = new double[M_rows * M_cols];

    double worksize;
    int lwork = -1;
    int liwork = 1;

    // Calculate workspace size
    pdsyevd_(&jobz, &uplo, &N, H, &IA, &IA, DESC, W, V, &IA, &IA,
             DESC_V, &workN, &lwork, &liwork, &liwork, &ierr);
    lwork = (int)worksize;

    double* work = new double[lwork];
    int* iwork = new int[liwork];

    pdsyevd_(&jobz, &uplo, &N, H, &IA, &IA, DESC, W, V, &IA, &IA,
            DESC_V, work, &lwork, iwork, &liwork, &ierr);

    delete[] work;
    delete[] iwork;

    if(ierr) myabort("diago: diagonization failed.");

    delete[] V;
    blacs_gridexit_(&ICTXT);
    if(!rank) output(H, W, N);
    std::cout << "[diago] end\n";
    mytimer::end("diago");
}
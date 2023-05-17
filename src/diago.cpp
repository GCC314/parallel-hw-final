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
// void sl_init_(int*, const int*, const int*);
// void blacs_gridinfo_(const int*, int*, int*, int*, int*);
int  numroc_(const int*, const int*, const int*, const int*, const int*);
void descinit_(int*, const int*, const int*, const int*, const int*, const int*, const int*,
               const int*, const int*, int*);
// void blacs_gridexit_(const int*);
void pdsyevd_(const char* JOBZ, const char* UPLO, const int* N, double* A, const int* IA, const int* JA,
              const int* DESCA, double* W, double* Z, const int* IZ, const int* JZ, const int* DESCZ,
              double* WORK, const int* LWORK, int* IWORK, const int* LIWORK, int* INFO);
void Cblacs_gridinit(int *context, char *order, int nprow, int npcol);
int Csys2blacs_handle(MPI_Comm comm);
void Cblacs_pcoord(int context, int pnum, int *prow, int *pcol);
void Cblacs_gridexit(int context);

int indxg2p_(int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);
void pdgemr2d_(int *m, int *n, double *A, int *ia, int *ja,
               int *desca, double *B, int *ib, int *jb,
               int *descb, int *gcontext);


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

    if(rank == 0) std::cout << "[diago] start\n";
    if(rank == 0) mytimer::start("diago");
    
    if(mode == D_MODE_LAP){ // LAPACK
        if(rank == 0){
            int res = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U',
                N, H, N, W);
            if(res) myabort("diago: diagonization failed.");
            output(H, W, N);
            std::cout << "[diago] end\n";
            mytimer::end("diago");
        }
        return;
    }
    else if(mode != D_MODE_SCA) // UNKNOWN
        myabort("diago_lib: not supported.");

    // SCALAPACK
    int B = 2; // block size
    int zero = 0, one = 1;
    char jobz = 'V', uplo = 'U';
    int info;
    double *A, *Z;
    int descA[9], descZ[9];

    // Initialize BLACS
    int blacs_context = Csys2blacs_handle(MPI_COMM_WORLD);
    int proc_rows = (int)sqrt(size);
    int proc_cols = size / proc_rows;
    while(proc_rows * proc_cols != size) {
        proc_rows--;
        proc_cols = size / proc_rows;
    }
    Cblacs_gridinit(&blacs_context, "Row-major", proc_rows, proc_cols);

    // Compute local array sizes
    int myrow = 0, mycol = 0;
    Cblacs_pcoord(blacs_context, rank, &myrow, &mycol);
    int local_rows = numroc_(&N, &B, &myrow, &zero, &proc_rows);
    int local_cols = numroc_(&N, &B, &mycol, &zero, &proc_cols);

    // Allocate local arrays
    A = (double*)malloc(local_rows * local_cols * sizeof(double));
    Z = (double*)malloc(local_rows * local_cols * sizeof(double));

    // Initialize matrix A
    for(int i=0; i<N; i++) {
        for(int j=0; j<N; j++) {
            int global_row = i+1;
            int global_col = j+1;
            int local_row = indxg2p_(&global_row, &B, &myrow, &zero, &proc_rows);
            int local_col = indxg2p_(&global_col, &B, &mycol, &zero, &proc_cols);
            if(local_row == myrow && local_col == mycol) {
                int local_index = (global_row-1)/proc_rows*B + (global_col-1)/proc_cols*B*local_rows + (global_row-1)%B + (global_col-1)%B*B;
                A[local_index] = H[i * N + j];
            }
        }
    }

    // Initialize array descriptors
    int lda = std::max(1, local_rows);
    descinit_(descA, &N, &N, &B, &B, &zero, &zero,
                &blacs_context, &lda, &info);
    descinit_(descZ, &N, &N, &B, &B, &zero, &zero,
                &blacs_context, &lda, &info);

    // Compute eigenvalues and eigenvectors

    int lwork = -1;
    int liwork = -1;
    double work_query;
    int iwork_query;
    pdsyevd_(&jobz,&uplo,&N,A,&one,&one,
            descA,W,Z,&one,&one,
            descZ,&work_query,&lwork,&iwork_query,&liwork,&info);

    lwork = (int)work_query;
    liwork = iwork_query;
    double* work = (double*)malloc(lwork*sizeof(double));
    int* iwork = (int*)malloc(liwork*sizeof(int));
    pdsyevd_(&jobz,&uplo,&N,A,&one,&one,
            descA,W,Z,&one,&one,
            descZ,work,&lwork,iwork,&liwork,&info);

    int rsrc = 0, csrc = 0;
    int izero = 0;
    int lld = std::max(1, local_rows);
    pdgemr2d_(&N,&N,Z,&one,&one,descZ,
            H,&one,&one,descA,
            &blacs_context);

    free(work); free(iwork);
    free(A); free(Z);

    Cblacs_gridexit(blacs_context);


    if(info) myabort("diago: diagonization failed.");

    if(!rank) output(H, W, N);
    for(int i = 0;i < N;i++) printf("%d: w[%d]=%f\n",rank,i,W[i]);
    if(rank == 0) std::cout << "[diago] end\n";
    if(rank == 0) mytimer::end("diago");
}
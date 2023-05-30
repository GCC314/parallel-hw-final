#include "diago.h"
#include "myerror.h"
#include "mpi.h"
#include <lapacke.h>
#include <cstring>
#include <iomanip>
#include <cstdio>
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

const int TAG_ZBUF = 314, TAG_ZINFO = 272;
char MAJOR_TYPE[] = "Row-major";

void output(double *H, double *W, int N){
    {
        std::ofstream ofs("output/eigenvalues.log");
        if(!ofs.is_open()) myabort("Failed to write eigenvalues.log");
        ofs << "Eigenvalues:\n";
        for(int i = 0;i < N;i++) ofs << W[i] << '\n';
    }
    {
        std::ofstream ofs("output/eigenvalues.log");
        if(!ofs.is_open()) myabort("Failed to write eigenvectors.log");
        ofs << "Eigenvectors (row-wise):\n";
        for(int i = 0;i < N;i++) for(int j = 0;j < N;j++){
            ofs << std::setw(6) << H[i * N + j] << " \n"[j == N - 1];
        }
    }
}

inline int local_to_global(int idx, int bsize, int pcoord, int plen){
    return ((idx / bsize) * plen + pcoord) * bsize + (idx % bsize);
}

// Too long and too annoying, treat in a function
void diagonize_scalapack(int rank, int size, double *H, double *W, int N, 
    int mode){
    int B = 2; // block size, matrix split into 2*2
    int zero = 0, one = 1;
    char jobz = 'V', uplo = 'U';
    int info;
    double *A, *Z;
    int descA[9], descZ[9];

    // Initialize BLACS
    int blacs_context = Csys2blacs_handle(MPI_COMM_WORLD);
    
    // split process grid, nprow * npcol == size;
    int nprow = (int)sqrt(size);
    int npcol = size / nprow;
    while(nprow * npcol != size) {
        nprow--;
        npcol = size / nprow;
    }

    Cblacs_gridinit(&blacs_context, MAJOR_TYPE, nprow, npcol);
    // Initialize the process grid

    // get coordinate in process grid
    int myrow, mycol;
    Cblacs_pcoord(blacs_context, rank, &myrow, &mycol);

    int local_nrow = numroc_(&N, &B, &myrow, &zero, &nprow);
    int local_ncol = numroc_(&N, &B, &mycol, &zero, &npcol);
    // fprintf(stderr, "rank=%d, lnr=%d, lnc=%d, mr=%d, mc=%d\n", \
    //     rank, local_nrow, local_ncol, myrow, mycol);

    // allocate the local matrix
    A = (double*)malloc(sizeof(double) * local_nrow * local_ncol);
    Z = (double*)malloc(sizeof(double) * N * N);

    // initialize the local matrix
    for (int i = 0; i < local_nrow; i++){
        for (int j = 0; j < local_ncol; j++){
            int global_i = local_to_global(i, B, myrow, nprow);
            int global_j = local_to_global(j, B, mycol, npcol);
            // fprintf(stderr, "(%d,%d)->(%d,%d)\n",i,j,global_i,global_j);
            A[i + j * local_nrow] = H[global_i + global_j * N];
        }
    }

    // for(int i = 0;i < local_nrow;i++) for(int j = 0;j < local_ncol;j++){
    //     fprintf(stderr, "rank=%d,A[%d][%d]=%.1lf\n", \
    //         rank, i, j, A[i + j * local_nrow]);
    // }

    // initialize descs
    int lda = std::max(1, local_nrow);
    descinit_(descA, &N, &N, &B, &B, &zero, &zero,
                &blacs_context, &lda, &info);
    descinit_(descZ, &N, &N, &B, &B, &zero, &zero,
                &blacs_context, &lda, &info);

    // optimize work distribution
    int lwork = -1;
    int liwork = -1;
    double work_query;
    int iwork_query;
    pdsyevd_(&jobz,&uplo,&N,A,&one,&one,
            descA,W,Z,&one,&one,
            descZ,&work_query,&lwork,&iwork_query,&liwork,&info);
    
    if(info) myabort("diago: work distribution failed");

    // Compute!
    lwork = (int)work_query;
    liwork = iwork_query;
    double* work = (double*)malloc(lwork*sizeof(double));
    int* iwork = (int*)malloc(liwork*sizeof(int));
    pdsyevd_(&jobz,&uplo,&N,A,&one,&one,
            descA,W,Z,&one,&one,
            descZ,work,&lwork,iwork,&liwork,&info);

    if(info) myabort("diago: diagonization failed");

    // reduction of Z
    double *Zbuf = (double*)malloc(sizeof(double) * N * N);
    
    if(rank){
        // distribute
        int Zinfo[4] = {local_nrow, local_ncol, myrow, mycol};
        MPI_Send(Zinfo, 4, MPI_INT, 0, TAG_ZINFO, MPI_COMM_WORLD);
        MPI_Send(Z, N * N, MPI_DOUBLE, 0, TAG_ZBUF, MPI_COMM_WORLD);
    }

    if(rank == 0){
        // reduce
        for(int id = 0;id < size;id++){
            int Zinfo[4]; // {local_nrow, local_ncol, myrow, mycol};
            if(id){
                MPI_Status sta;
                MPI_Recv(Zinfo, 4, MPI_INT, id, TAG_ZINFO, MPI_COMM_WORLD, &sta);
                MPI_Recv(Zbuf, N * N, MPI_DOUBLE, id, TAG_ZBUF, MPI_COMM_WORLD, &sta);
            }else{
                Zinfo[0] = local_nrow; Zinfo[1] = local_ncol;
                Zinfo[2] = myrow; Zinfo[3] = mycol;
                memcpy(Zbuf, Z, sizeof(double) * local_ncol * local_nrow);
            }
            
            // Write back to H
            for (int i = 0; i < Zinfo[0]; i++){
                for (int j = 0; j < Zinfo[1]; j++){
                    int global_i = local_to_global(i, B, Zinfo[2], nprow);
                    int global_j = local_to_global(j, B, Zinfo[3], npcol);
                    H[global_i + global_j * N] = Zbuf[i + j * Zinfo[0]];
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // End
    free(work); free(iwork);
    free(A); free(Z); free(Zbuf);

    Cblacs_gridexit(blacs_context);
}

void diagonize(int rank, int size, double *H, double *W, int N, 
    int mode){

    if(rank == 0){
        std::cout << "[diago] start\n";
        mytimer::start("diago");
    }
    
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

    diagonize_scalapack(rank, size, H, W, N, mode); // SCALAPACK

    if(rank == 0){
        output(H, W, N);
        std::cout << "[diago] end\n";
        mytimer::end("diago");
    }
}
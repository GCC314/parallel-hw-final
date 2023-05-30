LAPACK_ROOT_DIR = thirdparty/lapack-3.11
LAPACK_LIB_DIR = $(LAPACK_ROOT_DIR)/build
LAPACK_HEADER_DIR = $(LAPACK_ROOT_DIR)/LAPACKE/include

SCALAPACK_ROOT_DIR = thirdparty/scalapack-2.2.0
SCALAPACK_LIB = $(SCALAPACK_ROOT_DIR)/libscalapack.a

LIBS = -L$(LAPACK_LIB_DIR) -L$(SCALAPACK_ROOT_DIR) -lscalapack -llapacke -llapack -lblas
# LIBS = -L$(LAPACK_LIB_DIR) -L$(SCALAPACK_ROOT_DIR) -lscalapack -llapacke -llapack -lblacs -lpblas -lblas

CFLAGS = -std=c++11 -O3 -fopenmp
LFLAGS = -lm -lgfortran
MPIFLAGS = 

MPICXX = mpicxx
MPIRUN = mpirun
MPI_PROCESS_NUM = 2
OPENMP_THREAD_NUM = 4

EXEC_NAME = main
EXEC_DIR = build/

INPUT_FILE_DIR = data/INPUT.txt

.PHONY: main clean lapack scalapack run

run:
	$(MPIRUN) -np $(MPI_PROCESS_NUM) $(MPIRUN_FLAGS) $(EXEC_DIR)/$(EXEC_NAME) $(INPUT_FILE_DIR) $(OPENMP_THREAD_NUM)

main:
	$(MPICXX) -D __MPI -o $(EXEC_DIR)/$(EXEC_NAME) $(CFLAGS) -I$(LAPACK_HEADER_DIR) src/*.cpp $(LIBS) $(LFLAGS)

lapack:
	$(MAKE) -C $(LAPACK_ROOT_DIR) lib

scalapack:
	$(MAKE) -C $(SCALAPACK_ROOT_DIR) lib

clean:
	- $(MAKE) -C $(LAPACK_ROOT_DIR) clean
	- $(MAKE) -C $(SCALAPACK_ROOT_DIR) clean
	- rm $(EXEC_DIR)/$(EXEC_NAME)
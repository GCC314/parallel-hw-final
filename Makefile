LAPACK_ROOT_DIR = thirdparty/lapack-3.11
LAPACK_LIB_DIR = $(LAPACK_ROOT_DIR)/build
LAPACK_HEADER_DIR = $(LAPACK_ROOT_DIR)/LAPACKE/include

SCALAPACK_ROOT_DIR = thirdparty/scalapack-2.2.0
SCALAPACK_LIB = $(SCALAPACK_ROOT_DIR)/libscalapack.a

CFLAGS = -std=c++11 -O3 -fopenmp
LFLAGS = -llapacke -llapack -lm -lgfortran
MPICXX = mpicxx

EXEC_NAME = main
EXEC_DIR = build/

.PHONY: main clean lapack scalapack

main:
	$(MPICXX) -D __MPI -o $(EXEC_DIR)/$(EXEC_NAME) $(CFLAGS) $(LFLAGS) -I$(LAPACK_HEADER_DIR) -L$(LAPACK_LIB_DIR) src/*.cpp $(SCALAPACK_LIB)

lapack:
	$(MAKE) -C $(LAPACK_ROOT_DIR) lapacklib
	$(MAKE) -C $(LAPACK_ROOT_DIR) lapackelib

scalapack:
	$(MAKE) -C $(SCALAPACK_ROOT_DIR) blacslib
	$(MAKE) -C $(SCALAPACK_ROOT_DIR) scalapacklib

clean:
	- $(MAKE) -C $(LAPACK_ROOT_DIR) clean
	- $(MAKE) -C $(SCALAPACK_ROOT_DIR) clean
	- rm $(EXEC_NAME)
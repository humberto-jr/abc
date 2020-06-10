# Intel ifort + Intel MKL:
# make PES_OBJECT=filename.o
#
# Intel ifort + MAGMA:
# make use_magma PES_OBJECT=filename.o
#
# Humberto Jr
# Jun, 2020

SHELL = /bin/sh

CC = gcc
FC = ifort
PES_OBJECT =
LINEAR_ALGEBRA = MKL

#
# Intel MKL:
#

ifeq ($(LINEAR_ALGEBRA), MKL)
	MKL_LIB=$(MKLROOT)/lib/intel64
	MKL_INC=$(MKLROOT)/mkl/include

	LINEAR_ALGEBRA_LIB = -mkl=parallel -L$(MKL_LIB) -I$(MKL_INC) -Wl,--start-group $(MKL_LIB)/libmkl_intel_lp64.a $(MKL_LIB)/libmkl_intel_thread.a $(MKL_LIB)/libmkl_core.a $(MKL_LIB)/libmkl_blas95_lp64.a $(MKL_LIB)/libmkl_lapack95_lp64.a -Wl,--end-group -fopenmp -lpthread
endif

#
# MAGMA library:
#

MAGMAROOT = /usr/local/magma
CUDAROOT = /usr/local/cuda

ifeq ($(LINEAR_ALGEBRA), MAGMA)
	LINEAR_ALGEBRA_INC = -I$(CUDAROOT)/include -I$(MAGMAROOT)/include -DADD_
	LINEAR_ALGEBRA_LIB = -L$(MAGMAROOT)/lib -L$(CUDAROOT)/lib64 -lmagma -lm
endif

#
# Rules:
#

all: abc link
use_magma: abc dgemm dsyev syminv link

abc:
	$(FC) -O3 *.f90 -c
	@echo

link:
	$(FC) *.o -o abc.out $(PES_OBJECT) $(LINEAR_ALGEBRA_LIB)
	@echo

WRAPPERS_DIR = wrappers

dgemm: $(WRAPPERS_DIR)/dgemm.c $(WRAPPERS_DIR)/wrappers.h $(WRAPPERS_DIR)/c_lib.h
	$(CC) -W -Wall -std=c99 -pedantic -O3 $(LINEAR_ALGEBRA_INC) -c $<
	@echo

dsyev: $(WRAPPERS_DIR)/dsyev.c $(WRAPPERS_DIR)/wrappers.h $(WRAPPERS_DIR)/c_lib.h
	$(CC) -W -Wall -std=c99 -pedantic -O3 $(LINEAR_ALGEBRA_INC) -c $<
	@echo

syminv: $(WRAPPERS_DIR)/syminv.c $(WRAPPERS_DIR)/wrappers.h $(WRAPPERS_DIR)/c_lib.h
	$(CC) -W -Wall -std=c99 -pedantic -O3 $(LINEAR_ALGEBRA_INC) -c $<
	@echo

clean:
	rm -f *.o *.out

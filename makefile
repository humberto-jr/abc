# How to:
#
# First compile the external user defined PES into
# object files *.o (-c flag). Then, ABC is built
# for many scenarios as follows:
#
# 1) Intel ifort + Intel MKL (default).
# make
#
# Option (1) requires that the env variable MKLROOT
# is properly defined. Check the content of
# $MKLROOT first.
#
# If instead GNU gfortran is intended:
#
# 2) GNU gfortran + Intel MKL.
# make FC=gfortran
#
# Make sure both the PES and ABC are compiled with
# the same kind of compiler in order to avoid missing
# libraries during the link.
#
# 3) Intel ifort + Intel MKL + MAGMA.
# make use_magma LINEAR_ALGEBRA=MAGMA
#
# If MAGMA is not installed in standard locations,
# you can provide the path to its components:
#
# 4) Intel ifort + Intel MKL + MAGMA.
# make use_magma LINEAR_ALGEBRA=MAGMA MAGMAROOT=[path/to/magma-2.5.2] CUDAROOT=[path/to/cuda/10.1]
#
# Wrappers for MAGMA are written in C and are built
# with GNU gcc. If Intel icc is needed, then:
#
# 5) Intel ifort + Intel icc + Intel MKL + MAGMA.
# make use_magma CC=icc LINEAR_ALGEBRA=MAGMA MAGMAROOT=[path/to/magma-2.5.2] CUDAROOT=[path/to/cuda/10.1]
#
# 6) GNU gfortran + Intel icc + Intel MKL + MAGMA.
# make use_magma FC=gfortran CC=icc LINEAR_ALGEBRA=MAGMA MAGMAROOT=[path/to/magma-2.5.2] CUDAROOT=[path/to/cuda/10.1]
#
# Humberto Jr
# Jun, 2020

SHELL = /bin/sh

CC = gcc
FC = ifort
LINEAR_ALGEBRA = MKL

#
# Intel MKL:
#

MKL_LIB=$(MKLROOT)/lib/intel64
MKL_INC=$(MKLROOT)/mkl/include

ifeq ($(FC), ifort)
	LINEAR_ALGEBRA_LIB = -mkl=parallel -L$(MKL_LIB) -I$(MKL_INC) -Wl,--start-group $(MKL_LIB)/libmkl_intel_lp64.a $(MKL_LIB)/libmkl_intel_thread.a $(MKL_LIB)/libmkl_core.a $(MKL_LIB)/libmkl_blas95_lp64.a $(MKL_LIB)/libmkl_lapack95_lp64.a -Wl,--end-group -fopenmp -lpthread
endif

ifeq ($(FC), ifx)
	LINEAR_ALGEBRA_LIB = -qmkl=parallel -L$(MKL_LIB) -I$(MKL_INC) -Wl,--start-group $(MKL_LIB)/libmkl_intel_lp64.a $(MKL_LIB)/libmkl_intel_thread.a $(MKL_LIB)/libmkl_core.a $(MKL_LIB)/libmkl_blas95_lp64.a $(MKL_LIB)/libmkl_lapack95_lp64.a -Wl,--end-group -fopenmp -lpthread
endif

ifeq ($(FC), gfortran)
	LINEAR_ALGEBRA_LIB = -Wl,--start-group $(MKL_LIB)/libmkl_gf_lp64.a $(MKL_LIB)/libmkl_sequential.a $(MKL_LIB)/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
endif

#
# MAGMA library:
#

MAGMAROOT = /usr/local/magma
CUDAROOT = /usr/local/cuda

ifeq ($(LINEAR_ALGEBRA), MAGMA)
	LINEAR_ALGEBRA_INC += -I$(CUDAROOT)/include -I$(MAGMAROOT)/include -DADD_
	LINEAR_ALGEBRA_LIB += -L$(MAGMAROOT)/lib -L$(CUDAROOT)/lib64 -lmagma -lm
endif

#
# Rules:
#

all: abc link
use_magma: abc dgemm dsyev syminv dsyr link

abc:
	$(FC) -O3 *.f90 -c
	@echo

link:
	$(FC) *.o -o abc.out $(LINEAR_ALGEBRA_LIB)
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

dsyr: $(WRAPPERS_DIR)/dsyr.c $(WRAPPERS_DIR)/wrappers.h $(WRAPPERS_DIR)/c_lib.h
	$(CC) -W -Wall -std=c99 -pedantic -O3 $(LINEAR_ALGEBRA_INC) -c $<
	@echo

clean:
	rm -f *.o *.out

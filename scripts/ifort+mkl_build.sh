#!/bin/bash

if [ ! -e $1 ]
then
	echo
	echo "$0, error: $1 not found"
	echo
	exit 666
fi

if [ ! -e $MKLROOT ]
then
	echo
	echo "$0, error: MKLROOT undefined"
	echo
	exit 666
else
	MKL_LIB=$MKLROOT/lib/intel64
	MKL_INC=$MKLROOT/mkl/include
fi

#
# Serial:
#

#ifort -O3 *.f90 -o abc.out $1 -L$MKL_LIB -I$MKL_INC -Wl,--start-group $MKL_LIB/libmkl_intel_lp64.a $MKL_LIB/libmkl_sequential.a $MKL_LIB/libmkl_core.a $MKL_LIB/libmkl_blas95_lp64.a $MKL_LIB/libmkl_lapack95_lp64.a -Wl,--end-group -lpthread

#
# Multi-threaded:
#

ifort -O3 *.f90 -o abc.out $1 -mkl=parallel -L$MKL_LIB -I$MKL_INC -Wl,--start-group $MKL_LIB/libmkl_intel_lp64.a $MKL_LIB/libmkl_intel_thread.a $MKL_LIB/libmkl_core.a $MKL_LIB/libmkl_blas95_lp64.a $MKL_LIB/libmkl_lapack95_lp64.a -Wl,--end-group -fopenmp -lpthread

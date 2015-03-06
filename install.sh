#!/bin/bash
MKL_LIB=/home/programmes/intel/mkl/11.0.4.183/mkl/lib/intel64
MKL_INC=/home/programmes/intel/mkl/11.0.4.183/mkl/include
#
# Serial:
#
#ifort -O3 abc.f fun.f lin.f pot.f -o abc.out -L$MKL_LIB -I$MKL_INC -Wl,--start-group $MKL_LIB/libmkl_intel_lp64.a $MKL_LIB/libmkl_sequential.a $MKL_LIB/libmkl_core.a $MKL_LIB/libmkl_blas95_lp64.a $MKL_LIB/libmkl_scalapack_lp64.a -Wl,--end-group -lpthread
#
# Multi-threaded:
#
ifort -O3 abc.f fun.f lin.f pot.f -o abc.out -L$MKL_LIB -I$MKL_INC -Wl,--start-group $MKL_LIB/libmkl_intel_lp64.a $MKL_LIB/libmkl_intel_thread.a $MKL_LIB/libmkl_core.a $MKL_LIB/libmkl_blas95_lp64.a $MKL_LIB/libmkl_scalapack_lp64.a -Wl,--end-group -openmp -lpthread

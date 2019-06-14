#!/bin/bash

ScratchPath='/hydra/ptmp/tjyderri/octopus-svn-15667/trunk'

module purge
module load intel/14.0
module load mkl/11.1
module load mpi.ibm/1.3.0
module load gsl
module load fftw/3.3.4

./configure CC=mpiicc CFLAGS='-O3 -xHost -openmp' FC=mpiifort FCFLAGS='-O3 -xHost -openmp' LDFLAGS='-L/u/system/SLES11/soft/intel1/14.0/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm -Xlinker -rpath=/u/system/SLES11/soft/intel1/14.0/mkl/lib/intel64:/u/system/SLES11/soft/gsl/1.15/lib' --prefix=$ScratchPath --enable-mpi --disable-libnbc --disable-gdlib --with-gsl-prefix=/u/system/SLES11/soft/gsl/1.15 --with-libxc-prefix=/hydra/u/system/SLES11/soft/libxc/2.0.3/intel-14.0/mpi.ibm-1.3 --enable-openmp --with-fftw-prefix='/u/system/SLES11/soft/fftw/3.3.4/intel-14.0/mpi.ibm-1.3' --with-fftw-lib='-L/u/system/SLES11/soft/fftw/3.3.4/intel-14.0/mpi.ibm-1.3/lib -lfftw3' --with-blas='-L/u/system/SLES11/soft/intel1/14.0/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm' --with-lapack='-L/u/system/SLES11/soft/intel1/14.0/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm' --with-blacs='-L/u/system/SLES11/soft/intel1/14.0/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm' --with-scalapack='-L/u/system/SLES11/soft/intel1/14.0/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm'

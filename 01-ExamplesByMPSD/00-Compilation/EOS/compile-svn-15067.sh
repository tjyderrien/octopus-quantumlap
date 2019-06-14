#!/bin/bash

echo "/!\ Use autoreconf -i before running this script."
# autoreconf -i

echo "Compilation of Octopus for EOS"

echo "Loading modules..."
module purge
module load intel/15.0
module load mkl/11.1
module load impi/5.1.1
module load gsl/1.16
module load fftw/3.3.4
module load nfft/3.2.4

echo "Starting compilation... using line given in /u/system/SLES11/soft/octopus/5.0.0/intel-15.0/impi-5.1.1/.build.doc/config.log"

# ./configure CC=mpiicc CFLAGS='-O3 -xCORE-AVX2 -fma -ip -openmp -fzero-initialized-in-bss -init=arrays -init=zero' FC=mpiifort FCFLAGS='-O3 -xCORE-AVX2 -fma -ip -openmp -fzero-initialized-in-bss -init=arrays -init=zero' LDFLAGS='-L/afs/@cell/common/soft/intel/ics2015/15.0/compiler/lib/intel64 -liomp5 -L/afs/@cell/common/soft/intel/ics2015/15.0/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm -Xlinker -rpath=/afs/@cell/common/soft/intel/ics2015/15.0/mkl/lib/intel64:/afs/@cell/.cs/gsl/1.16/@sys/lib:/u/system/SLES11/soft/pfft/1.0.8a/intel-15.0/impi-5.1.1/lib:/afs/rzg/common/soft/netcdf/4.3.2/amd64_sles11/intel/14.0/serial/lib' --prefix=$HOME/octopus-svn-15067-nfft --enable-mpi --enable-openmp --enable-avx --disable-gdlib --with-gsl-prefix=/afs/@cell/.cs/gsl/1.16/@sys --with-libxc-prefix=/u/system/SLES11/soft/libxc/2.1.2/intel-15.0/impi-5.1.1 --with-fftw-prefix='/afs/rzg/.cs/fftw/fftw-3.3.4/@sys/intel-15.0/impi-5.0/' --with-blas='-L/afs/@cell/common/soft/intel/ics2015/15.0/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm' --with-lapack='-L/afs/@cell/common/soft/intel/ics2015/15.0/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm' --with-blacs='-L/afs/@cell/common/soft/intel/ics2015/15.0/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm' --with-scalapack='-L/afs/@cell/common/soft/intel/ics2015/15.0/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm' --with-pfft-prefix=/u/system/SLES11/soft/pfft/1.0.8a/intel-15.0/impi-5.1.1 --with-netcdf-prefix=/afs/rzg/common/soft/netcdf/4.3.2/amd64_sles11/intel/14.0/serial --with-nfft=/scratch/tjyderri/nfft/3.2.4

make -j
# make install

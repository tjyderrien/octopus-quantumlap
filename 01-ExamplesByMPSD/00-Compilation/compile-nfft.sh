#!/bin/bash

module purge
module load fftw/3.3.4
# module load intel/15.0
# module load mkl/11.1

FFTW_HOME="/afs/rzg/.cs/fftw/fftw-3.3.4/@sys/intel-15.0/impi-5.0"
#export CFLAGS="$(CFLAGS) -O3 -D GAUSSIAN"
# export LDFLAGS="$(LDFLAGS) -L/afs/@cell/common/soft/intel/ics2015/15.0/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm -Xlinker -rpath=/afs/@cell/common/soft/intel/ics2015/15.0/mkl/lib/intel64:/afs/@cell/.cs/gsl/1.16/@sys/lib:/u/system/SLES11/soft/pfft/1.0.8a/intel-15.0/impi-5.1.1/lib:/afs/rzg/common/soft/netcdf/4.3.2/amd64_sles11/intel/14.0/serial/lib"

# ./configure LDFLAGS='-L/afs/@cell/common/soft/intel/ics2015/15.0/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm -Xlinker -rpath=/afs/@cell/common/soft/intel/ics2015/15.0/mkl/lib/intel64:/afs/@cell/.cs/gsl/1.16/@sys/lib:/u/system/SLES11/soft/pfft/1.0.8a/intel-15.0/impi-5.1.1/lib:/afs/rzg/common/soft/netcdf/4.3.2/amd64_sles11/intel/14.0/serial/lib' 

# Thibault's attempt. 
# ./configure -lm -lmkl_core -lmkl_intel_thread -lm  --with-fftw3=${FFTW_HOME} --with-fftw3-libdir=${FFTW_HOME}/lib --with-fftw3-includedir=${FFTW_HOME}/include -L/afs/@cell/common/soft/intel/ics2015/15.0/mkl/lib/intel64

# Heiko's 
module load intel/15.0 mkl/11.1 impi/5.1.1; ./configure --enable-all --prefix=$HOME CC=icc LDLFAGS='-L/afs/@cell/common/soft/intel/ics2015/15.0/compiler/lib/intel64 -liomp5 -L/afs/@cell/common/soft/intel/ics2015/15.0/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm' --with-fftw3=/afs/rzg/.cs/fftw/fftw-3.3.4/@sys/intel-15.0/impi-5.0/

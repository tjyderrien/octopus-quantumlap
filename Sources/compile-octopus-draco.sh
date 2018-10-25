#!/bin/bash

. /etc/profile.d/modules.sh

module purge  
module load intel/17.0
module load impi/2017.4
module load mkl/2017.4
module load libxc/3.0
module load gsl/2.3  
module load fftw-serial/3.3.8
module load hdf5-serial/1.8.20
module load netcdf-serial/4.4.1
module list

export MKL="-L$MKL_HOME/lib/intel64  -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm"

export CC=mpiicc 
export CFLAGS="-g -O3  -xHost -I$FFTW_HOME/include" #-openmp
export CXX=mpiicpc
export FC=mpiifort 
export FCFLAGS="-g -O3  -xHost -I$FFTW_HOME/include" #-openmp

export LDFLAGS="-I$FFTW_HOME/include -L$MKL_HOME/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core  -lmkl_sequential -lmkl_blacs_intelmpi_lp64  -lpthread -lm -Xlinker -rpath=$MKL_HOME/lib/intel64:$GSL_HOME/lib:$NETCDF_HOME/lib" 

#export LD_LIBRARY_PATH=$ISF_HOME/lib:$LD_LIBRARY_PATH 
make mostlyclean
make clean
make distclean

autoreconf -i
#NOTE: OpenMP was disabled, because not useful in our case.
./configure --enable-mpi  --disable-gdlib --disable-openmp \
   --enable-avx \
   --with-gsl-prefix=$GSL_HOME \
   --with-libxc-prefix=$LIBXC_HOME \
   --with-blas="$MKL" \
   --with-lapack="$MKL" \
   --with-blacs="$MKL" \
   --with-scalapack="$MKL" \
   --with-fftw-prefix="$FFTW_HOME" \
   --with-netcdf-prefix=$NETCDF_HOME \
   --prefix="/u/tjyderri/octopus.git/"

 #  --with-isf-prefix=$ISF_HOME \
 #  --with-isf-include=$ISF_HOME/include \

# make clean 
make -j
make install
#MPIEXEC=srun
#make check-short >& check-short.log

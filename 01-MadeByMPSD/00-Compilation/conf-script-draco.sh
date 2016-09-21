#!/bin/bash

. /etc/profile.d/modules.sh

module purge  
module load intel/16.0
module load impi/5.1.3
module load mkl/11.3
module load libxc/2.2.3
module load gsl/2.1  
module load fftw/3.3.4
module load netcdf-serial/4.4.0
module list

export MKL="-L$MKL_HOME/lib/intel64  -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm"

export CC=mpiicc 
export CFLAGS="-g -O3  -xHost -openmp -I$FFTW_HOME/include"
export CXX=mpiicpc
export FC=mpiifort 
export FCFLAGS="-g -O3  -xHost -openmp -I$FFTW_HOME/include" 

export LDFLAGS="-I$FFTW_HOME/include -L$MKL_HOME/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core  -lmkl_sequential -lmkl_blacs_intelmpi_lp64  -lpthread -lm -Xlinker -rpath=$MKL_HOME/lib/intel64:$GSL_HOME/lib:$NETCDF_HOME/lib:" 

#export LD_LIBRARY_PATH=$ISF_HOME/lib:$LD_LIBRARY_PATH 

autoreconf -i

./configure --enable-mpi  --disable-gdlib --enable-openmp \
   --enable-avx \
   --with-gsl-prefix=$GSL_HOME \
   --with-libxc-prefix=$LIBXC_HOME \
   --with-blas="$MKL" \
   --with-lapack="$MKL" \
   --with-blacs="$MKL" \
   --with-scalapack="$MKL" \
   --with-fftw-prefix="$FFTW_HOME" \
   --with-netcdf-prefix=$NETCDF_HOME \
   --prefix="/u/tnico/octopus-trunk/"

 #  --with-isf-prefix=$ISF_HOME \
 #  --with-isf-include=$ISF_HOME/include \

make clean 
make -j
make install


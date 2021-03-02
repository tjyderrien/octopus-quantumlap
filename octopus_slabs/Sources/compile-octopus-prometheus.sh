#!/bin/bash

module purge  
# module load icc/2017.7.259-GCC-6.3.0-2.27
# module load plgrid/tools/intel/17.0.5
module load plgrid/tools/intel/18.0.0 
# module load ifort/2017.7.259-GCC-6.3.0-2.27
# module load impi/2017.4.239-iccifort-2017.7.259-GCC-6.3.0-2.27
module load plgrid/tools/impi/2018

# module load imkl/2017.4.239-iimpi-2017c 
module load plgrid/libs/mkl/2018.0.0

# module load libxc/4.2.3-intel-2017a
module load plgrid/libs/libxc/4.2.3

# module load GSL/2.5-intel-2017a
module load plgrid/libs/gsl/2.4

# module load FFTW/3.3.6-intel-2017a 
module load plgrid/libs/fftw/3.3.7

# module load HDF5/1.10.2-intel-2017a
module load plgrid/libs/hdf5/1.10.4

# module load netCDF/4.4.1.1-intel-2017a
module load plgrid/libs/netcdf/4.7.0

module list

#Currently Loaded Modules:
#  1) plgrid/libs/libxc/4.2.3    4) plgrid/libs/fftw/3.3.7      7) plgrid/tools/impi/2018
#  2) plgrid/libs/gsl/2.3        5) plgrid/tools/gcc/6.4.0      8) plgrid/libs/hdf5/1.10.4
#  3) plgrid/libs/mkl/2018.0.0   6) plgrid/tools/intel/18.0.0   9) plgrid/libs/netcdf/4.7.0

 


# export LIBXC_HOME=/net/people/prace/pr2c0043/libxc
export LIBXC_HOME=/net/software/local/libxc/4.2.3

export MKL="-L$MKL_HOME/lib/intel64  -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm"

export CC=mpiicc 
export CFLAGS="-g -O3  -xHost -I$FFTW_HOME/include" #-openmp
export CXX=mpiicpc
export FC=mpiifort 
export FCFLAGS="-g -O3  -xHost -I$FFTW_HOME/include" #-openmp

export LDFLAGS="-I$FFTW_HOME/include -L$MKL_HOME/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core  -lmkl_sequential -lmkl_blacs_intelmpi_lp64  -lpthread -lm -Xlinker -rpath=$MKL_HOME/lib/intel64:$GSL_HOME/lib:$NETCDF_HOME/lib:$FFTW_HOME/lib" 

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
   --prefix="/net/people/prace/pr2c0043/octopus.git/" 

 #  --with-isf-prefix=$ISF_HOME \
 #  --with-isf-include=$ISF_HOME/include \

# make clean 
make -j 4
make install
#MPIEXEC=srun
#make check-short >& check-short.log

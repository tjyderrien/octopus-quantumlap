#!/bin/bash

set -e

module purge  

## Build for Intel Fortran 2017. #Warning: Outdated from Octopus commit d3529c86 (Jan 31st, 2020).
ImportLibsFortran2017() {

module load Autotools/20180311
# module load icc/2017.7.259-GCC-6.3.0-2.27
module load icc/2017.1.132-GCC-6.3.0-2.27

# module load ifort/2017.7.259-GCC-6.3.0-2.27
module load ifort/2017.1.132-GCC-6.3.0-2.27

module load iccifort/2017.1.132-GCC-6.3.0-2.27

# module load impi/2017.4.239-iccifort-2017.7.259-GCC-6.3.0-2.27
module load impi/2017.1.132-iccifort-2017.1.132-GCC-6.3.0-2.27
module load iimpi/2017a

module load imkl/2017.1.132-iimpi-2017a
# module load imkl/2017.4.239-iimpi-2017c 

export MKL_HOME=/apps/all/imkl/2017.1.132-iimpi-2017a/mkl/lib/intel64:/apps/all/imkl/2017.1.132-iimpi-2017a #May 21st, 8:12am <derrien@fzu.cz> 
module load libxc/4.2.3-intel-2017a
module load GSL/2.5-intel-2017a
export GSL_HOME=/apps/all/GSL/2.5-intel-2017a #May 21st, 8:20am <derrien@fzu.cz>
module load FFTW/3.3.6-intel-2017a #(old CentSO6 module)
export FFTW_HOME=/apps/all/FFTW/3.3.6-intel-2017a #May 21st, 8:19am <derrien@fzu.cz>
module load HDF5/1.10.2-intel-2017a
# module load netCDF/4.4.1.1-intel-2017a
# module load netCDF/4.6.3-intel-2017a
# export NETCDF_HOME=/apps/all/netCDF/4.6.3-intel-2017a #May 21st, 8:21am <derrien@fzu.cz>
}

## Build with Intel Fortran from 2020
ImportLibsFortran2020() {
module purge
module load Autotools/20180311
module load iccifort/2020.0.166
module load iimpi/2020a
module load imkl/2020.0.166-iimpi-2020a
export MKL_HOME=/apps/all/imkl/2020.0.166-iimpi-2020a #May 21st, 8:12am <derrien@fzu.cz> 
module load libxc/4.3.4-intel-2020a
module load GSL/2.5-intel-2020a
export GLS_HOME=/apps/all/GSL/2.5-intel-2020a
module load FFTW/3.3.8-intel-2020a
export FFTW_HOME=/apps/all/FFTW/3.3.8-intel-2020a
module load HDF5/1.10.6-intel-2020a
}

ImportLibsFortran2021() {
module purge
module load Autotools/20200321-GCCcore-10.2.0
module load iccifort/2020.4.304
module load iimpi/2020b
module load imkl/2020.4.304-iimpi-2020b
export MKL_HOME=/apps/all/imkl/2020.4.304-iimpi-2020b
module load libxc/4.3.4-intel-2020a #sends iccifort back to 2020.1.217 and iimpi.2020a #WARNING: no longer valid since March 2022. 
module load GSL/2.6-iccifort-2020.1.217
export GLS_HOME=/apps/all/GSL/2.6-iccifort-2020.1.217
module load FFTW/3.3.8-intel-2020b
export FFTW_HOME=/apps/all/FFTW/3.3.8-intel-2020b
module load HDF5/1.10.7-iimpi-2020b

# module load METIS/5.1.0-intel-2020a
module load SuiteSparse/5.8.1-intel-2020b-METIS-5.1.0 #since March 2022
# module load ParMETIS/4.0.3-intel-2020a

# module load CUDA/10.2.89
# export CUDA_FLAGS="--enable-cuda --enable-nvtx --with-cuda-prefix=$CUDA_HOME"

export CUDA_FLAGS=""

}

ImportLibsFortran2021

## Valid for both
# module load libyaml/0.2.2
# export YAML_HOME=/apps/all/libyaml/0.2.2

module list

export MKL="-L$MKL_HOME/lib/intel64  -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm"

export CC=mpiicc
export CFLAGS="-g -O3  -xHost -I$FFTW_HOME/include" #-openmp
export CXX=mpiicpc
export FC=mpiifort 
export FCFLAGS="-g -O3  -xHost -I$FFTW_HOME/include" #-openmp

# export LDFLAGS="-I$FFTW_HOME/include -L$MKL_HOME/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core  -lmkl_sequential -lmkl_blacs_intelmpi_lp64  -lpthread -lm -Xlinker -rpath=$MKL_HOME/lib/intel64:$GSL_HOME/lib:$NETCDF_HOME/lib:$FFTW_HOME/lib" 

# May 21st, 8:12am <derrien@fzu.cz> #NETCDF fails anyway. 
export LDFLAGS="-I$FFTW_HOME/include -L$MKL_HOME/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core  -lmkl_sequential -lmkl_blacs_intelmpi_lp64  -lpthread -lm -Xlinker -rpath=$MKL_HOME/lib/intel64:$GSL_HOME/lib:$FFTW_HOME/lib:$LIBRARY_PATH"

#export LD_LIBRARY_PATH=$ISF_HOME/lib:$LD_LIBRARY_PATH 
make mostlyclean
make clean
make distclean

autoreconf -i
#NOTE: OpenMP was disabled, because not useful in our case.
./configure $CUDA_FLAGS \
   FCFLAGS_FFTW="-I$MKL_HOME/include/fftw" \
   FCCPP="cpp -ffreestanding" \
   --enable-mpi  --disable-gdlib --disable-openmp \
   --enable-avx \
   --with-gsl-prefix=$GSL_HOME \
   --with-libxc-prefix=$LIBXC_HOME \
   --with-blas="$MKL" \
   --with-lapack="$MKL" \
   --with-blacs="$MKL" \
   --with-scalapack="$MKL" \
   --with-yaml-prefix="$YAML_HOME/lib" \
   --prefix="/home/tjyderrien/octopus.me/"

#   --prefix="/home/gazdovak/octopus.bin/"
#    --with-parmetis-prefix="$PARMETIS_HOME"
#   --with-metis-prefix="$METIS_HOME" \
#   --with-netcdf-prefix=$NETCDF_HOME \ #attempting in manual because this loads $NETCDF_HOME/lib not $NETCDF_HOME/lib64
#    --with-fftw-prefix="$FFTW_HOME" \
 #  --with-isf-prefix=$ISF_HOME \
 #  --with-isf-include=$ISF_HOME/include \
#   --with-yaml-prefix="$YAML_HOME/lib" \

# make clean 
make -j 12
make install
# MPIEXEC=mpirun
#make check-short >& check-short.log

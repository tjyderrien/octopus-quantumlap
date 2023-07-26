#!/bin/bash

module purge  
# module load IntelOneAPI/compiler/2021.1.1
# module load IntelOneAPI/intel_ipp_intel64/2021.1.1
# module load IntelOneAPI/intel_ippcp_intel64/2021.1.1

module load GCC/12.1.0

# module load ifort/2017.7.259-GCC-6.3.0-2.27
# module load impi/2017.4.239-iccifort-2017.7.259-GCC-6.3.0-2.27
# module load plgrid/tools/impi/2018
# module load IntelOneAPI/mpi/2021.1.1
module load OpenMPI/4.1.3

# module load imkl/2017.4.239-iimpi-2017c 
# module load IntelOneAPI/mkl/2021.1.1
# module load IntelOneAPI/tbb/2021.1.1
module load BLAS/3.10.0

module load LAPACK/3.10.1
export LAPACK_HOME="/APPS/centos7/centos79/software/Compiler/GCC-12.1/LAPACK/3.10.1"
# module load plgrid/libs/libxc/4.2.3

module load libxc/6.2.2
# module load libxc/5.2.3-GCC-8.1.0
export LIBXC_HOME="/APPS/centos7/centos79/software/Compiler/GCC-12.1/libxc/6.2.2"

# module load plgrid/libs/gsl/2.4
module load GSL/2.7
export GSL_HOME="/APPS/centos7/centos79/software/Compiler/GCC-12.1/GSL/2.7"

# module load plgrid/libs/fftw/3.3.7
module load FFTW3/3.3.10
export FFTW_HOME="/APPS/centos7/centos79/software/MPI/GCC-12.1-OpenMPI-4.1/FFTW3/3.3.10"
# module load plgrid/libs/hdf5/1.10.4
module load HDF5/1.12.0
# module load plgrid/libs/netcdf/4.7.0
# module load netCDF
module list

#Currently Loaded Modules:
#  1) plgrid/libs/libxc/4.2.3    4) plgrid/libs/fftw/3.3.7      7) plgrid/tools/impi/2018
#  2) plgrid/libs/gsl/2.3        5) plgrid/tools/gcc/6.4.0      8) plgrid/libs/hdf5/1.10.4
#  3) plgrid/libs/mkl/2018.0.0   6) plgrid/tools/intel/18.0.0   9) plgrid/libs/netcdf/4.7.0

 


# export LIBXC_HOME=/net/people/prace/pr2c0043/libxc
# export LIBXC_HOME=/net/software/local/libxc/4.2.3

# export MKL="-L$MKL_HOME/lib/intel64  -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm"

# export CC=mpiicc 
# export CFLAGS="-g -O3  -xHost -I$FFTW_HOME/include" #-openmp
# export CXX=mpiicpc
# export FC=mpiifort 



# export MKL="-L$MKL_HOME/lib/intel64  -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm"

# export CC=mpicc 
# export CFLAGS="-g -O3  -xHost -I$FFTW_HOME/include" #-openmp
# export CXX=mpi++
# export FC=mpifort 

# export FCFLAGS="-g -O3  -xHost -I$FFTW_HOME/include" #-openmp

# export LDFLAGS="-I$FFTW_HOME/include -L$MKL_HOME/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core  -lmkl_sequential -lmkl_blacs_intelmpi_lp64  -lpthread -lm -Xlinker -rpath=$MKL_HOME/lib/intel64:$GSL_HOME/lib:$NETCDF_HOME/lib:$FFTW_HOME/lib" 

make mostlyclean
make clean
make distclean

cd .. 
spack load autoconf
autoreconf -i
cd -

spack load automake
#NOTE: OpenMP was disabled, because not useful in our case.
../configure --enable-mpi  --disable-gdlib --disable-openmp \
    --with-gsl-prefix=$GSL_HOME \
    --with-libxc-prefix=$LIBXC_HOME \
    --with-blas="$BLAS" \
    --with-lapack-prefix="$LAPACK" \
    --with-fftw-prefix="$FFTW_HOME" \
    --prefix="/gpfs/home/tderrien/octopus.gnu.mpiio/"
#   --enable-avx \ 
#   --with-lapack="$MKL" \
#   --with-netcdf-prefix=$NETCDF_HOME \
#  --with-isf-prefix=$ISF_HOME \
#  --with-isf-include=$ISF_HOME/include \
#   --with-blacs="$MKL" \
#   --with-scalapack="$MKL" \


make -j 16
make install

#MPIEXEC=srun
 #make check-short >& check-short.log


#!/bin/bash

echo "LOADING SPACK..."
cd ~
. spack/share/spack/setup-env.sh
cd -

spack load autoconf@2.71%gcc@=12.1.0
# spack load automake@1.13.4%gcc@=12.1.0
# spack load automake@1.16.5
spack load automake@1.16.5%gcc@=12.1.0
spack load m4@1.4.19%gcc@=12.1.0
spack load libtool

module purge  
module load GCC/12.1.0

module load OpenMPI/4.1.3
# spack load openmpi

module load BLAS/3.10.0
# spack load openblas@0.3.23%gcc@=12.1.0
# export BLAS="/gpfs/home/tderrien/spack/opt/spack/linux-centos7-x86_64_v4/gcc-12.1.0/openblas-0.3.23-kex72bdrxntzmj2wya4vgnq443gqgeyv/lib/libopenblas.so"

# OpenBLAS does not work as BLAS, maybe GCC is too new. Should work with GCC 9.3.
# HOWEVER this version replaces Lapack (it passes configure).
# spack load openblas@0.3.23%gcc@=12.1.0
# export LAPACK_HOME="/gpfs/home/tderrien/spack/opt/spack/linux-centos7-x86_64_v4/gcc-12.1.0/openblas-0.3.23-kex72bdrxntzmj2wya4vgnq443gqgeyv/lib/libopenblas.so"

# Edwin's lapack: fails: missing reference to a certain function: undefined reference to 'scnrm2_'
module load LAPACK/3.10.1
export LAPACK_HOME="/APPS/centos7/centos79/software/Compiler/GCC-12.1/LAPACK/3.10.1/liblapack.a"

# Netlib lapack with spack: worked well. 
# spack load netlib-lapack@3.11.0%gcc@=12.1.0
# export LAPACK_HOME="/gpfs/home/tderrien/spack/opt/spack/linux-centos7-x86_64_v4/gcc-12.1.0/netlib-lapack-3.11.0-h23ksrvrfen5cnmtegxfh5gkzaijrjb2/lib64/liblapack.so"

# spack load scalapack
# export SCALAPACK_HOME="/gpfs/home/tderrien/spack/opt/spack/linux-centos7-x86_64_v4/gcc-12.1.0/netlib-scalapack-2.2.0-etgr4pnwxtxqxr4yjmtl7q4pbmiqowbz/lib/libscalapack.so"

# module load libxc/6.2.2
# spack load libxc
module load libxc/5.2.3
export LIBXC_HOME="/APPS/centos7/centos79/software/Compiler/GCC-12.1/libxc/5.2.3"
# export LIBXC_HOME="/HOME/tderrien/spack/opt/spack/linux-centos7-x86_64_v4/gcc-12.1.0/libxc-5.2.3-fz4lueyi4clctlo56kzedbah3esizyfg"

module load GSL/2.7
export GSL_HOME="/APPS/centos7/centos79/software/Compiler/GCC-12.1/GSL/2.7"
# spack load gsl
# export GSL_HOME="/gpfs/home/tderrien/spack/opt/spack/linux-centos7-x86_64_v4/gcc-12.1.0/gsl-2.7.1-ndwqq3cpp37wz7klp5q6x7rh7m7es64v"

# Passes configure, but compilation fails. 
module load FFTW3/3.3.10
export FFTW_HOME="/APPS/centos7/centos79/software/MPI/GCC-12.1-OpenMPI-4.1/FFTW3/3.3.10"
export FFTW_HOME2="/APPS/centos7/centos79/software/MPI/FFTW/GCC-12.1-OpenMPI-4.1/FFTW3/3.3.10"
# spack load fftw
# export FFTW_HOME="/gpfs/home/tderrien/spack/opt/spack/linux-centos7-x86_64_v4/gcc-12.1.0/fftw-3.3.10-bsbgb7g2cu2h4aq67gc5dwllt7qlkm7i"

# module load HDF5/1.12.0
# module load plgrid/libs/netcdf/4.7.0
# module load netCDF
module list

echo $LD_LIBRARY_PATH

# export LIBXC_HOME=/net/people/prace/pr2c0043/libxc
# export LIBXC_HOME=/net/software/local/libxc/4.2.3

# export MKL="-L$MKL_HOME/lib/intel64  -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm"

# export CC=mpiicc 
#  export CFLAGS="-g -O3  -xHost -I$FFTW_HOME/include" #-openmp
# export CXX=mpiicpc
# export FC=mpiifort 

export CC=mpicc
export FC=mpif90
export CXX=mpicxx
export CFLAGS="-O3"
export FCFLAGS="$CFLAGS -ffree-line-length-none"
export CXXFLAGS="$CFLAGS"
export LDFLAGS="-Xlinker -rpath=$GSL_HOME/lib"	    

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
# echo "LOADING SPACK PACKAGES..."
# spack load autoconf@2.71%gcc@=12.1.0
# spack load m4@1.4.19%gcc@=12.1.0
# spack load automake@1.16.5%gcc@=12.1.0

echo "AUTORECONF..."
autoreconf -i -f
cd -

echo "STARTING CONFIGURATION..."
#NOTE: OpenMP was disabled, because not useful in our case.
../configure --enable-mpi  --disable-gdlib --disable-openmp \
    --with-libxc-prefix=$LIBXC_HOME \
    --with-blas="$BLAS" \
    --with-lapack-prefix="$LAPACK/lib64" \
    --with-gsl-prefix=$GSL_HOME \
    --with-fftw-prefix="$FFTW_HOME" \
    --prefix="/gpfs/home/tderrien/octopus.gnu.spack.mpiio/"

#    --with-lapack="$LAPACK_HOME" \

#    --with-lapack="$LAPACK/liblapack.a -llapack" \
#    --with-lapack="$LAPACK_HOME" \
#    --with-scalapack="$SCALAPACK_HOME" \
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


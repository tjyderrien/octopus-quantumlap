#!/bin/bash

module purge  
module load GCC/12.1.0

module load OpenMPI/4.1.3
module load BLAS/3.10.0

# Edwin's lapack: fails: missing reference to a certain function: undefined reference to 'scnrm2_'. No longer needed after update of BLAS by Edwin. 
module load LAPACK/3.10.1
# export LAPACK_HOME="/APPS/centos7/centos79/software/Compiler/GCC-12.1/LAPACK/3.10.1/liblapack.a"

# module load libxc/6.2.2
module load libxc/5.2.3
# module load libxc/5.2.3-GCC-8.1.0
export LIBXC_HOME="/APPS/centos7/centos79/software/Compiler/GCC-12.1/libxc/5.2.3"

module load GSL/2.7
export GSL_HOME="/APPS/centos7/centos79/software/Compiler/GCC-12.1/GSL/2.7"

# Passes configure, but compilation fails. 
module load FFTW3/3.3.10
export FFTW_HOME="/APPS/centos7/centos79/software/MPI/GCC-12.1-OpenMPI-4.1/FFTW3/3.3.10"

module load HDF5/1.12.0

# module load netCDF
module list

echo $LD_LIBRARY_PATH

export CC=mpicc
export FC=mpif90
export CXX=mpicxx
export CFLAGS="-g -pg -debug"
export FCFLAGS="$CFLAGS -ffree-line-length-none"
export CXXFLAGS="$CFLAGS"
export LDFLAGS="-Xlinker -rpath=$GSL_HOME/lib"	    

# export CC=mpicc 
# export CFLAGS="-g -O3  -xHost -I$FFTW_HOME/include" #-openmp
# export CXX=mpi++
# export FC=mpifort 

# export FCFLAGS="-g -O3  -xHost -I$FFTW_HOME/include" #-openmp

make mostlyclean
make clean
make distclean

cd ..
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
    --prefix="/gpfs/home/tderrien/octopus.gnu.mpiio/"



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


make -j 20
make install

#MPIEXEC=srun
#make check-short >& check-short.log


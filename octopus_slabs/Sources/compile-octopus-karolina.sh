#!/bin/sh

set -e

module purge  

 if [[ ! -f ../configure.ac ]]; then
   echo "Error! Please execute this script in a subfolder of the root directory."
   exit 1
 fi

ImportLibsAocc2021()  {
module purge

echo "Enable shell support for spack..."
. /home/tjyderrien/spack/share/spack/setup-env.sh

# module load Autotools/20200321-GCCcore-10.2.0
# module load libtool/2.4.6-GCCcore-10.2.0
# module load M4/1.4.19

echo "Loading AMD compilers + MPI suite..."
# module load AOCC/3.1.0-GCCcore-10.2.0
# module load AOCL/3.0.1-AOCC-3.1.0
# module load OpenMPI/4.1.1-GCC-10.2.0-AOCL-3.0.1-AOCC-3.1.0
module load OAMD/2021a
# spack load openmpi@4.1.1

echo "Loading LibXC..."
# module load libxc/4.3.4-iccifort-2020.4.304
# export LIBXC_HOME=/apps/all/libxc/4.3.4-iccifort-2020.4.304/
# spack load libxc@5.1.5
spack load libxc@4.3.4
module load libxc-4.3.4-aocc-3.1.0-bcpjlku
# export LIBXC_HOME=/home/tjyderrien/spack/opt/spack/linux-centos7-zen2/aocc-3.1.0/libxc-5.1.5-plerszrmfxohojimlwegi3pivw65ok2g/
export LIBXC_HOME=/home/tjyderrien/spack/opt/spack/linux-centos7-zen2/aocc-3.1.0/libxc-4.3.4-bcpjlkupgvfmjxrgdx7pfzg7bvpfbmwy/

echo "Loading GSL..."
# module load GSL/2.7-GCC-10.3.0
# export GSL_HOME=/apps/all/GSL/2.7-GCC-10.3.0/
spack load gsl@2.7
module load gsl-2.7-aocc-3.1.0-syclsmo
export GSL_HOME=/home/tjyderrien/spack/opt/spack/linux-centos7-zen2/aocc-3.1.0/gsl-2.7-syclsmoeyhdg4xtmyp37ko5afymfaxgl/

# echo "FlexiBLAS..."
# module load FlexiBLAS/3.0.4-GCC-10.3.0

# echo "Loading OpenBLAS..."
# spack load openblas
# module load openblas-0.3.18-aocc-3.1.0-ckpuajx

# module load FFTW/3.3.10
# export FFTW_HOME=/apps/all/FFTW/3.3.9-gompi-2021a/
# spack load fftw@3.3.10
# export FFTW_HOME=/home/tjyderrien/spack/opt/spack/linux-centos7-zen2/aocc-3.1.0/fftw-3.3.10-7jdy5byzctb6lwpvyjcuqndh56pelpku/

}

ImportIntelOnKarolina() {

module purge

# module load intel-compilers/2021.2.0

## MKL Intel2020b (compatible with LibXC 4)
module load imkl/2020.4.304-iimpi-2020b
#export MKL_HOME=$MKLROOT:$MKLROOT/lib/intel64:


## MKL old
# export MKL_HOME="/apps/all/imkl/2021.2.0-iimpi-2021a/mkl/latest/"
export MKL_HOME="$MKLROOT"
export MKL="-I$MKLROOT/include -L$MKL_HOME/lib/intel64  -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl"
export FFTW_HOME="${MKL_HOME}/include/fftw/"

# From EasyBuild
export BLAS_LIB_DIR="/apps/all/imkl/2020.4.304-iimpi-2020b/mkl/lib/intel64" 
export LIBBLAS="-Wl,-Bstatic -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -Wl,-Bdynamic"
export BLACS_HOME="$MKLROOT/lib/intel64/libmkl_blacs_intelmpi_lp64.a"
export SCALAPACK_HOME="$MKLROOT/lib/intel64/libmkl_scalapack_lp64.a"


## LibXC: let's keep libxc 4 for the version I need. 
module load libxc/4.3.4-intel-2020b
export LIBXC_HOME=/apps/all/libxc/4.3.4-intel-2020b/
# module load libxc/5.1.5-intel-compilers-2021.2.0
# export LIBXC_HOME="/apps/all/libxc/5.1.5-intel-compilers-2021.2.0/"

module load GSL/2.6-iccifort-2020.4.304
export GSL_HOME="/apps/all/GSL/2.6-iccifort-2020.4.304/"

}


# ImportLibsAocc2021
ImportIntelOnKarolina

## Valid for both
# module load libyaml/0.2.5-GCCcore-10.2.0
# export YAML_HOME=/apps/all/libyaml/0.2.5-GCCcore-10.2.0
# spack load libyaml@0.2.5
# export YAML_HOME=/home/tjyderrien/spack/opt/spack/linux-centos7-zen2/aocc-3.1.0/libyaml-0.2.5-jn4av4foia3gigrk4klwt7lzzxgayevb/


module list

export CC=mpicc
export FC=mpiifort
export CXX=mpicxx

## AOCL
# export CFLAGS="-g -O3 -m64 -march=znver1 -mfma -fvectorize -mfma -mtune=znver1 -mavx2 -m3dnow -fuse-ld=lld -I$AOCL_ROOT/include" 
#-openmp
# export FCFLAGS="-g -O3 -m64 -march=znver1 -mfma -fvectorize -mfma -mtune=znver1 -mavx2 -m3dnow -fuse-ld=lld -I$AOCL_ROOT/include" 
FCFLAGS="-O2 -march=core-avx2" 
# try: -march=znver2 or haswell #https://prace-ri.eu/wp-content/uploads/Best-Practice-Guide_AMD.pdf
#-openmp -ffast-math

# Intel
export CFLAGS="-g -O3  -xHost -I$MKL_HOME/include -I$FFTW_HOME/include"
# export FCFLAGS="$CFLAGS" #-g -O3  -xHost -I$FFTW_HOME/include"

## Common (intel|amd)
export CXXFLAGS="$CFLAGS"

# export LDFLAGS="-I$FFTW_HOME/include -L$MKL_HOME/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core  -lmkl_sequential -lmkl_blacs_intelmpi_lp64  -lpthread -lm -Xlinker -rpath=$MKL_HOME/lib/intel64:$GSL_HOME/lib:$NETCDF_HOME/lib:$FFTW_HOME/lib" 


## AOCL
# May 21st, 8:12am <derrien@fzu.cz> #NETCDF fails anyway. 
# export LDFLAGS="-I$FFTW_HOME/include -L$AOCL_ROOT/lib -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core  -lmkl_sequential -lmkl_blacs_intelmpi_lp64  -lpthread -lm -Xlinker -rpath=$MKL_HOME/lib/intel64:$GSL_HOME/lib:$FFTW_HOME/lib:$LIBRARY_PATH"

# Nov 10th 2021 <derrien@fzu.cz>
#export LDFLAGS="-I$AOCL_ROOT/include -L$AOCL_ROOT/lib -lalm -lm -Xlinker -rpath=$AOCL_ROOT/lib:$GSL_HOME/lib:$SPACK_LD_LIBRARY_PATH:$LD_LIBRARY_PATH:$PATH"

#export LD_LIBRARY_PATH=$ISF_HOME/lib:$LD_LIBRARY_PATH 
# make mostlyclean || echo "mostlyclean: not needed."
# make clean || echo "clean: not needed."
# make distclean || echo "distclean: not needed"

## Intel
export LDFLAGS="-I$FFTW_HOME/include -L$MKL_HOME/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core  -lmkl_sequential -lmkl_blacs_intelmpi_lp64  -lpthread -lm -Xlinker -rpath=$MKL_HOME/lib/intel64:$GSL_HOME/lib:$FFTW_HOME/lib:$LIBRARY_PATH"

echo "Autoreconf..."
pushd .. && autoreconf -i --force && popd
echo "Configure..."
../configure $CUDA_FLAGS \
   FCFLAGS_FFTW="-I$FFTW_HOME" \
   FCCPP="cpp -ffreestanding" \
   --enable-mpi  --disable-gdlib --disable-openmp \
   --with-gsl-prefix="$GSL_HOME" \
   --with-libxc-prefix="$LIBXC_HOME" \
   --with-blas="$BLAS_LIB_DIR $LIBBLAS" \
   --with-blacs="$BLACS_HOME" \
   --with-scalapack="$SCALAPACK_HOME" \
   --prefix='/home/tjyderrien/octopus.intel.mini'

## Intel
#   --with-fftw="$MKL" \
#   --with-lapack="$MKL" \

## AOCL
#   --prefix='/home/tjyderrien/octopus.aocc'
#   --enable-avx \ 
#   --with-fftw="$AOLC_ROOT" \
#   --with-amdblis="$AOCL_ROOT" \
#   --with-lapack="$AOCL_ROOT" \
#   --with-scalapack="$AOCL_ROOT" \

#   --with-blas="$AOCL_ROOT" \ 
#   --with-scalapack="$AOCL_ROOT" \
#   --with-yaml-prefix="$YAML_HOME/lib", \
#   FCCPP="cpp -ffreestanding" \
#   FCFLAGS_FFTW="-I$FFTW_HOME/include/fftw" \  #this inline crashes configure and config.sub




#FFTW: otherwise, use -I$AOCL_ROOT/include/fftw
#   --host=x86_64-unknown-linux-gnu --build=x86_64-unknown-linux-gnu \
#   --with-fftw-prefix="$FFTW_HOME" \
#   --with-metis-prefix="$METIS_HOME" \
#   --with-parmetis-prefix="$PARMETIS_HOME" \
#   --enable-avx \ #maybe this is because of AMD ?

#   --with-netcdf-prefix=$NETCDF_HOME \ #attempting in manual because this loads $NETCDF_HOME/lib not $NETCDF_HOME/lib64

 #  --with-isf-prefix=$ISF_HOME \
 #  --with-isf-include=$ISF_HOME/include \

# make clean 
make -j 64
make install
# MPIEXEC=mpirun
#make check-short >& check-short.log

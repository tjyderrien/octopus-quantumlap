#!/bin/bash

set -e

module purge  


ImportLibsAocc2021()  {
module purge

# Enable shell support for spack
. /home/tjyderrien/spack/share/spack/setup-env.sh

module load Autotools/20180311
module load AOCC/3.1.0-GCCcore-10.2.0
module load AOCL/3.0.1-AOCC-3.1.0

module load OpenMPI/4.1.1-GCC-10.2.0-AOCL-3.0.1-AOCC-3.1.0
# spack load openmpi@4.1.1

# module load libxc/4.3.4-iccifort-2020.4.304
# export LIBXC_HOME=/apps/all/libxc/4.3.4-iccifort-2020.4.304/
spack load libxc@5.1.5
export LIBXC_HOME=/home/tjyderrien/spack/opt/spack/linux-centos7-zen2/aocc-3.1.0/libxc-5.1.5-plerszrmfxohojimlwegi3pivw65ok2g/

# module load GSL/2.7-GCC-10.3.0
# export GSL_HOME=/apps/all/GSL/2.7-GCC-10.3.0/
spack load gsl@2.7
export GSL_HOME=/home/tjyderrien/spack/opt/spack/linux-centos7-zen2/aocc-3.1.0/gsl-2.7-syclsmoeyhdg4xtmyp37ko5afymfaxgl/

# module load FFTW/3.3.10
# export FFTW_HOME=/apps/all/FFTW/3.3.9-gompi-2021a/
spack load fftw@3.3.10
export FFTW_HOME=/home/tjyderrien/spack/opt/spack/linux-centos7-zen2/aocc-3.1.0/fftw-3.3.10-7jdy5byzctb6lwpvyjcuqndh56pelpku/

}

ImportLibsFortran2021() {
module purge
module load Autotools/20180311

module load iccifort/2020.4.304
module load iimpi/2020b

module load imkl/2020.4.304-iimpi-2020b
export MKL_HOME=/apps/all/imkl/2020.4.304-iimpi-2020b

module load libxc/4.3.4-iccifort-2020.4.304
export LIBXC_HOME=/apps/all/libxc/4.3.4-iccifort-2020.4.304/

#TODO Maybe install / ask for GSL/2.6-iccifort-2020.4.304 ? 
module load GSL/2.6-iccifort-2020.4.304
export GSL_HOME=/apps/all/GSL/2.6-iccifort-2020.4.304

# module load FFTW/3.3.8-gompi-2020a-amd
# export FFTW_HOME=/apps/all/FFTW/3.3.8-gompi-2020a-amd
module load FFTW/3.3.8-iccifort-2020.4.304
export FFTW_HOME=/apps/all/FFTW/3.3.8-iccifort-2020.4.304

# module load HDF5/1.12.0-iimpi-2020a #not available
# module load HDF5/1.10.7-iimpi-2020b 
# module load METIS/5.1.0-intel-2020a #not available
# module load METIS/5.1.0-GCCcore-10.2.0
# export METIS_HOME=/apps/all/METIS/5.1.0-GCCcore-10.2.0
# module load ParMETIS/4.0.3-intel-202i0a
# module load ParMETIS/4.0.3-intel-2020b
# export PARMETIS_HOME=/apps/all/ParMETIS/4.0.3-intel-2020b

# module load CUDA/10.2.89
# export CUDA_FLAGS="--enable-cuda --enable-nvtx --with-cuda-prefix=$CUDA_HOME"
export CUDA_FLAGS=""
}

# ImportLibsFortran2021 #This is not for AMD Karolina
ImportLibsAocc2021

## Valid for both
# module load libyaml/0.2.5-GCCcore-10.2.0
# export YAML_HOME=/apps/all/libyaml/0.2.5-GCCcore-10.2.0
spack load libyaml@0.2.5
export YAML_HOME=/home/tjyderrien/spack/opt/spack/linux-centos7-zen2/aocc-3.1.0/libyaml-0.2.5-jn4av4foia3gigrk4klwt7lzzxgayevb/


module list

# export MKL="-L$MKL_HOME/lib/intel64  -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm"

export CC=mpicc
export FC=mpif90
export CXX=mpicxx
export CFLAGS="-g -O3 -march=znver1 -mfma -fvectorize -mfma -mtune=znver1 -mavx2 -m3dnow -fuse-ld=lld -I$FFTW_HOME/include" 
#-openmp
export FCFLAGS="-g -O3 -march=znver1 -mfma -fvectorize -mfma -mtune=znver1 -mavx2 -m3dnow -fuse-ld=lld -I$FFTW_HOME/include" 
# try: -march=znver2 or haswell #https://prace-ri.eu/wp-content/uploads/Best-Practice-Guide_AMD.pdf
#-openmp -ffast-math
export CXXFLAGS="$CFLAGS"

# export LDFLAGS="-I$FFTW_HOME/include -L$MKL_HOME/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core  -lmkl_sequential -lmkl_blacs_intelmpi_lp64  -lpthread -lm -Xlinker -rpath=$MKL_HOME/lib/intel64:$GSL_HOME/lib:$NETCDF_HOME/lib:$FFTW_HOME/lib" 

# May 21st, 8:12am <derrien@fzu.cz> #NETCDF fails anyway. 
# export LDFLAGS="-I$FFTW_HOME/include -L$AOCL_ROOT/lib -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core  -lmkl_sequential -lmkl_blacs_intelmpi_lp64  -lpthread -lm -Xlinker -rpath=$MKL_HOME/lib/intel64:$GSL_HOME/lib:$FFTW_HOME/lib:$LIBRARY_PATH"
# Nov 10th 2021 <derrien@fzu.cz>
export LDFLAGS="-I$FFTW_HOME/include -L$AOCL_ROOT/lib -lalm -lm -Xlinker -rpath=$AOCL_ROOT/lib:$GSL_HOME/lib:$FFTW_HOME/lib:$LD_LIBRARY_PATH"

#export LD_LIBRARY_PATH=$ISF_HOME/lib:$LD_LIBRARY_PATH 
# make mostlyclean || echo "mostlyclean: not needed."
# make clean || echo "clean: not needed."
# make distclean || echo "distclean: not needed"

echo "Autoreconf..."
autoreconf -i
./configure $CUDA_FLAGS \
   FCFLAGS_FFTW="-I$FFTW_HOME/include/fftw" \ 
   FCCPP="cpp -ffreestanding" \
   --target=x86_64-unknown-linux-gnu \
   --enable-mpi  --disable-gdlib --disable-openmp \
   --enable-avx \
   --with-gsl-prefix="$GSL_HOME" \
   --with-libxc-prefix="$LIBXC_HOME" \
   --with-blas="$AOCL_ROOT" \
   --with-lapack="$AOCL_ROOT" \
   --with-blacs="$AOCL_ROOT" \
   --with-scalapack="$AOCL_ROOT" \
   --with-yaml-prefix="$YAML_HOME/lib", \
   --prefix='/home/tjyderrien/octopus.aocc'

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

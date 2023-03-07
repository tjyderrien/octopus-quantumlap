#!/bin/bash

module purge  
# module load Autotools/20210128-GCCcore-10.3.0
module load Autotools/20200321-GCCcore-10.2.0

ImportNavGNULibs() {
	module load libxc/5.1.3-GCC-10.2.0
	export LIBXC_HOME="/apps/software/libxc/5.1.3-GCC-10.2.0/"

	module load OpenBLAS/0.3.12-GCC-10.2.0
	export OPENBLAS_HOME="/apps/software/OpenBLAS/0.3.12-GCC-10.2.0/" #octopus seems to refuse this library. How to optimize with GNU then? 

	module load ScaLAPACK/2.1.0-gompic-2020b
	export SCALAPACK_HOME="/apps/software/ScaLAPACK/2.1.0-gompi-2020b/"
}

ImportNavIntelLibs() {

module load iccifort/2019.5.281
module load iimpi/2019b

source /apps/software/imkl/2019.5.281-iimpi-2019b/bin/compilervars.sh intel64 

module load imkl/2019.5.281-iimpi-2019b
export MKL_HOME="/apps/software/imkl/2019.5.281-iimpi-2019b/"

module load libxc/4.3.4-iccifort-2019.5.281
export LIBXC_HOME="/apps/software/libxc/4.3.4-iccifort-2019.5.281"

module load GSL/2.6-iccifort-2019.5.281
export GSL_HOME="/apps/software/GSL/2.6-iccifort-2019.5.281/"
}


ImportNavIntelLibs 

module load libyaml/0.2.2-GCCcore-8.3.0
export YAML_HOME=/apps/software/libyaml/0.2.2-GCCcore-8.3.0/

module list

export CC=mpiicc
export CXX=mpiicpc
export FC=mpiifort 
# export CFLAGS="-g -O3  -xHost -I$FFTW_HOME/include  -DMKL_ILP64  -I${MKLROOT}/include" #-openmp
# Draco advised line
export CFLAGS="-O3" # -qopt-zmm-usage=high -fma -ip"
export FCFLAGS=${CFLAGS}
export CXXFLAGS=${CFLAGS}

# Advised by Octopus website
# export MKL="-L$MKL_HOME/lib/intel64  -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm"
# See https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html: warning: octopus addresses C-API using 32 bits integer, not 64. Then one should use LP not ILP (https://www.smcm.iqfr.csic.es/docs/intel/mkl/mkl_userguide/MKL_UG_linking_your_application/Support_for_ILP64_Programming.htm)
export MKL=" -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl"
export FFTW_HOME=${MKLROOT}/include/fftw

# Draco advised line
# export LDFLAGS="-Xlinker -rpath=$MKL_HOME/lib/intel64:$GSL_HOME/lib:$NETCDF_HOME/lib:$ELPA_HOME/lib:$METIS_HOME/lib:$PARMETIS_HOME/lib"

#export LD_LIBRARY_PATH=$ISF_HOME/lib:$LD_LIBRARY_PATH 

# Custom line (works on Navigator)
# export LDFLAGS="-I$FFTW_HOME/include -L$MKL_HOME/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core  -lmkl_sequential -lmkl_blacs_intelmpi_lp64  -lpthread -lm -Xlinker -rpath=$MKL_HOME/lib/intel64:$GSL_HOME/lib:$NETCDF_HOME/lib:$FFTW_HOME/lib" 

# Dec 1st, 2021 <derrien@fzu.cz> #Replacing by another compilation line
export LDFLAGS="-I$FFTW_HOME $MKL -Xlinker -rpath=$MKLROOT/lib/intel64:$GSL_HOME/lib:$LIBRARY_PATH"

make mostlyclean
make clean
make distclean

pushd .. && autoreconf -i --force && popd
#NOTE: OpenMP was disabled, because not useful in our case.
../configure $CUDA_FLAGS \
   FCCPP="cpp -ffreestanding" \
   FCFLAGS_FFTW="-I$MKLROOT/include/fftw" \
   --enable-mpi  --disable-gdlib --disable-openmp \
   --with-gsl-prefix=$GSL_HOME \
   --with-libxc-prefix=$LIBXC_HOME \
   --with-blas="$MKL" \
   --with-lapack="$MKL" \
   --with-blacs="$MKL" \
   --with-scalapack="$MKL" \
   --with-yaml-prefix="$YAML_HOME/lib" \
   --prefix="/veracruz/home/p/pr2c0043/octopus.mpiio.July/"

#   --prefix="/veracruz/home/p/pr2c0055/octopus.mpiio/"
#   FCFLAGS_FFTW="-I$MKL_HOME/include/fftw" \ 
#    --with-parmetis-prefix="$PARMETIS_HOME"
#   --with-metis-prefix="$METIS_HOME" \
#   --with-netcdf-prefix=$NETCDF_HOME \ #attempting in manual because this loads $NETCDF_HOME/lib not $NETCDF_HOME/lib64
#    --with-fftw-prefix="$FFTW_HOME" \
 #  --with-isf-prefix=$ISF_HOME \
 #  --with-isf-include=$ISF_HOME/include \

# make clean 
make -j 48
make install
# MPIEXEC=mpirun
#make check-short >& check-short.log

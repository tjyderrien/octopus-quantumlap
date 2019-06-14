module purge  
module load intel/15.0
module load mkl/11.2
module load impi/5.1.1
module load gsl/1.16
module load fftw/3.3.4 
module load nfft/3.2.4
module list

export INSTALL_HERE=/scratch/tjyderri/octopus-svn-15246
export MKL="-L$MKL_HOME/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm"
export PFFT_HOME=/u/system/SLES11/soft/pfft/1.0.8a/intel-15.0/impi-5.1.1
#export NFFT_HOME=/scratch/appel/nfft/3.2.4
export LIBXC_HOME=/u/system/SLES11/soft/libxc/2.1.2/intel-15.0/impi-5.1.1
export FFTW4PFFT_HOME=/u/mjr/soft/rzg_software_stack/SLES11/fftw4pfft/fftw-3.3.4/
export NETCDF_HOME=/afs/rzg/common/soft/netcdf/4.3.2/amd64_sles11/intel/14.0/serial
export CC=mpiicc 
export CFLAGS="-O3 -xCORE-AVX2 -fma -ip -openmp -fzero-initialized-in-bss -I$NFFT_HOME/include -I$FFTW_HOME/include"
export FC=mpiifort 
export FCFLAGS="-O3 -xCORE-AVX2 -fma -ip -openmp -fzero-initialized-in-bss -init=arrays -init=zero -I$NFFT_HOME/include -I$FFTW_HOME/include"
export LDFLAGS="-I$FFTW_HOME/include -L$LIBXC_HOME/lib -L$MKL_HOME/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm -Xlinker -rpath=$MKL_HOME/lib/intel64:$GSL_HOME/lib:$NFFT_HOME/lib" 

export LDFLAGS="-L$MKL_HOME/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm -Xlinker -rpath=$MKL_HOME/lib/intel64:$GSL_HOME/lib:$PFFT_HOME/lib:$NETCDF_HOME/lib:$NFFT_HOME/lib:$LDFLAGS"


autoreconf -i

./configure --enable-mpi  --disable-gdlib --enable-openmp \
 --prefix=${INSTALL_HERE} \
 --with-gsl-prefix=$GSL_HOME \
 --with-libxc-prefix=$LIBXC_HOME \
 --with-blas="$MKL" \
 --with-lapack="$MKL" \
 --with-blacs="$MKL" \
 --with-scalapack="$MKL" \
 --with-fftw-prefix=$FFTW_HOME \
 --with-nfft=$NFFT_HOME \
 --with-pfft-prefix=$PFFT_HOME \
 --with-metis-prefix=/u/system/SLES11/soft/metis/5.1.0/intel-15.0 \
  --with-parmetis-prefix=/u/system/SLES11/soft/parmetis/4.0.3/intel-15.0/impi-5.1

make -j
make install 

if [[ ! -f ../configure.ac ]]; then
   echo "Error! Please execute this script in a subfolder of the root directory."
   exit 1
 fi
 
# compiler=intel
compiler=gnu
cuda=no
 
echo Using $compiler compiler.
[[ $cuda == yes ]] && echo Building with CUDA support.
 
module purge
if [[ $compiler == intel ]]; then
  # newest intel version
  module load intel/19.1.2 impi/2019.8 mkl/2020.2
  export CC=mpiicc
  export FC=mpiifort
  export CXX=mpiicpc
  export CFLAGS="-O3 -xCORE-AVX512 -qopt-zmm-usage=high -fma -ip"
  export FCFLAGS="$CFLAGS"
  export CXXFLAGS="$CFLAGS"
  export MKL="-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl"
elif [[ $compiler == gnu ]]; then
  module load gcc-9.3.0-gcc-11.1.0-edodhoo
  module load openmpi-4.1.1-gcc-9.3.0-tw2zexx
  module load openblas-0.3.18-gcc-9.3.0-3n6bxik
  module load fftw-3.3.10-gcc-9.3.0-7uc4kdi
  export FFTW_HOME="/opt-ext/spack/opt/spack/linux-centos7-broadwell/gcc-9.3.0/fftw-3.3.10-7uc4kdimd2h4ethl7zcur4iesmmlbbfz/"
  export CC=mpicc
  export FC=mpif90
  export CXX=mpicxx
  export CFLAGS="-O3"
  export FCFLAGS="$CFLAGS -ffree-line-length-none"
  export CXXFLAGS="$CFLAGS"
#  export MKL="-lgomp -lpthread -lm -ldl" #-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -lgomp -lpthread -lm -ldl"
  export BLAS_HOME="" #/opt-ext/spack/opt/spack/linux-centos7-broadwell/gcc-9.3.0/openblas-0.3.18-3n6bxikyxyv7plf7fljswyckps73i5de/" #fails if used. Works if empty. 
else
  echo "Compiler $compiler unknown."
  exit 1
fi

module load gsl-2.7-gcc-9.3.0-lpiv5vs
export GSL_HOME="/opt-ext/spack/opt/spack/linux-centos7-broadwell/gcc-9.3.0/gsl-2.7-lpiv5vskktfyqyvi5ennrfzsv47ny6zg/"

module load libxc-5.1.5-gcc-9.3.0-awi5re7
export LIBXC_HOME="/opt-ext/spack/opt/spack/linux-centos7-broadwell/gcc-9.3.0/libxc-5.1.5-awi5re7qx6ybyvltonphzpfxdwpsqvpl/"
# module load metis parmetis
export LDFLAGS="-Xlinker -rpath=$GSL_HOME/lib" #FFTW_HOME


if [[ $cuda == yes ]]; then
  module load cuda/10.2
  CUDA_FLAGS="--enable-cuda --enable-nvtx --with-cuda-prefix=$CUDA_HOME"
  LDFLAGS="$LDFLAGS:$CUDA_HOME/lib64"
else
  CUDA_FLAGS=""
fi
module list

export INSTALLDIR=/opt-ext/octopus.develop

if [[ "$1" != "noconfig" ]]; then
  pushd .. && autoreconf -i --force && popd  
  ../configure $CUDA_FLAGS \
    FCCPP="cpp -ffreestanding" \
    --prefix=$INSTALLDIR \
    --enable-mpi --enable-openmp \
    --disable-gdlib \
    --with-gsl-prefix="$GSL_HOME" \
    --with-libxc-prefix="$LIBXC_HOME" \
    || exit 1
fi
#    FCFLAGS_FFTW="-I$MKLROOT/include/fftw" \
#    --with-netcdf-prefix="$NETCDF_HOME" \
#    --with-metis-prefix="$METIS_HOME" \
#    --with-parmetis-prefix="$PARMETIS_HOME" \
#    --with-blas="$BLAS_HOME" \
#    --with-lapack="$BLAS_HOME" \
#    --with-blacs="$BLAS_HOME" \
#    --with-scalapack="$BLAS_HOME" \


echo "\n\nBuilding octopus...\n"

if [[ "$2" == "nodep" ]]; then
  make NODEP=1 -j20 && make NODEP=1 install || exit 1
else
  make -j20 && make install || exit 1
fi

mkdir -p $INSTALLDIR/.build.doc/
cp -f config.log $INSTALLDIR/.build.doc/
cp -f $0 $INSTALLDIR/.build.doc/

echo "... done"

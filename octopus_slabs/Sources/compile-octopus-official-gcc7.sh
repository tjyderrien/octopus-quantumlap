#!/bin/sh

module purge
# module load GCC/12.1.0
module load GCC/7.2.0

# module load OpenMPI/4.1.3
module load OpenMPI/3.1.6
module load BLAS/3.10.0

# Edwin's lapack: fails: missing reference to a certain function: undefined reference to 'scnrm2_'. No longer needed after update of BLAS by Edwin. 
module load LAPACK/3.10.1
# export LAPACK_HOME="/APPS/centos7/centos79/software/Compiler/GCC-12.1/LAPACK/3.10.1/liblapack.a"

# module load libxc/6.2.2
# module load libxc/5.2.3
module load libxc/4.3.4
# module load libxc/5.2.3-GCC-8.1.0
# export LIBXC_HOME="/APPS/centos7/centos79/software/Compiler/GCC-12.1/libxc/5.2.3"
# export LIBXC_HOME="/APPS/centos7/centos79/software/Compiler/GCC-12.1/libxc/4.3.4"
export LIBXC_HOME="/APPS/centos7/centos79/software/Compiler/GCC-7.2/libxc/4.3.4"
module load GSL/2.7
# export GSL_HOME="/APPS/centos7/centos79/software/Compiler/GCC-12.1/GSL/2.7"
export GSL_HOME="/APPS/centos7/centos79/software/Compiler/GCC-7.2/GSL/2.7"

# Passes configure, but compilation fails. 
module load FFTW3/3.3.10
# export FFTW_HOME="/APPS/centos7/centos79/software/MPI/GCC-12.1-OpenMPI-4.1/FFTW3/3.3.10"
export FFTW_HOME="/APPS/centos7/centos79/software/MPI/GCC-7.2-OpenMPI-3.1/FFTW3/3.3.10"
# module load HDF5/1.12.0
module load HDF5/1.10.2

# module load netCDF
module list

make mostlyclean
make clean
make distclean

cd ..
autoreconf -fi
cd -

export CC="mpicc"
MARCH_FLAG="-march=${GCC_ARCH:-native}"
OPTIMISATION_LEVEL="-O3"
export CFLAGS="$MARCH_FLAG $OPTIMISATION_LEVEL -g -fno-var-tracking-assignments"
export CXX="mpicxx"
export CXXFLAGS="$MARCH_FLAG $OPTIMISATION_LEVEL -g -fno-var-tracking-assignments"
export FC="mpif90"
export FCFLAGS="$MARCH_FLAG $OPTIMISATION_LEVEL -g -fno-var-tracking-assignments -ffree-line-length-none -Wno-argument-mismatch -Wno-invalid-pch"

#HG: ugly hack to include rpath while linking
#    becomes necessary for spack >= 0.19, as it does not set LD_LIBRARY_PATH anymore
export LDFLAGS=`echo ${LIBRARY_PATH:+:$LIBRARY_PATH} | sed -e 's/:/ -Wl,-rpath=/g'`

# HG: help configure to find the ELPA F90 modules and libraries
try_mpsd_elpa_version=`expr match "$MPSD_ELPA_ROOT" '.*/elpa-\([0-9.]\+\)'`
if [ -n "$try_mpsd_elpa_version" ] ; then
  if [ -r $MPSD_ELPA_ROOT/include/elpa_openmp-$try_mpsd_elpa_version/modules/elpa.mod ]; then
     export FCFLAGS_ELPA="-I$MPSD_ELPA_ROOT/include/elpa_openmp-$try_mpsd_elpa_version/modules"
     export LIBS_ELPA="-lelpa_openmp"
  elif [ -r $MPSD_ELPA_ROOT/include/elpa-$try_mpsd_elpa_version/modules/elpa.mod ] ; then
     export FCFLAGS_ELPA="-I$MPSD_ELPA_ROOT/include/elpa-$try_mpsd_elpa_version/modules"
  fi
fi
unset try_mpsd_elpa_version

# HG: always keep options in the order listed by ``configure --help``
../configure \
      --enable-mpi \
      --disable-openmp \
      --with-libxc-prefix="$LIBXC_HOME" \
      --with-blas="-L$BLAS/lib -lblas" \
      --with-gsl-prefix="$GSL_HOME" \
      --with-fftw-prefix="$FFTW_HOME" \
      --prefix="$HOME/octopus.parallel_restart.gnu7.libxc4.mpiio/" \
      "$@" | tee 00-configure.log 2>&1

#      --with-fftw-prefix="$FFTW_HOME" \ 
#      --with-libvdwxc-prefix="$MPSD_LIBVDWXC_ROOT" \
#      --with-pfft-prefix="$MPSD_PFFT_ROOT" \
#      --with-nfft="$MPSD_NFFT_ROOT" \
#      --with-pnfft-prefix="$MPSD_PNFFT_ROOT" \
#      --with-berkeleygw-prefix="$MPSD_BERKELEYGW_ROOT" \
#     --with-sparskit="$MPSD_SPARSKIT_ROOT/lib/libskit.a" \
#      --with-nlopt-prefix="$MPSD_NLOPT_ROOT" \
#      --with-blacs="-L$MPSD_NETLIB_SCALAPACK_ROOT/lib -lscalapack" \
#      --with-scalapack="-L$MPSD_NETLIB_SCALAPACK_ROOT/lib -lscalapack" \
#      --with-elpa-prefix="$MPSD_ELPA_ROOT" \
#      --with-cgal="$MPSD_CGAL_ROOT" \
#      --with-boost="$MPSD_BOOST_ROOT" \
#      --with-metis-prefix="$MPSD_METIS_ROOT" \
#      --with-parmetis-prefix="$MPSD_PARMETIS_ROOT" \

echo "-------------------------------------------------------------------------------" >&2
echo "configure output has been saved to 00-configure.log" >&2
if [ "x${GCC_ARCH-}" = x ] ; then
   echo "Microarchitecture optimization: native (set \$GCC_ARCH to override)" >&2
else
   echo "Microarchitecture optimization: $GCC_ARCH (from override \$GCC_ARCH)" >&2
fi
echo "-------------------------------------------------------------------------------" >&2

make -j 48
# make check-short
make install

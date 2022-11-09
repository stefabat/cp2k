# Build cp2k with cmake

This document regroups information about the cp2k cmake system. CMake is used to
detect cp2k dependencies and configure the compilation process. Dependencies
should be installed independently either with a distribution package manager,
easybuild, or spack to name a few.

It is easier to build and install all manually build dependencies in a single
directory ideally where cp2k will also be installed. Cmake will have less
difficulties to find the FindPACKAGE.cmake files and dependent libraries. CMake
will also use environment variables such as ORNL_FFTW3_ROOT, etc. Usually a
standard prefix is used in HPC environments. If known just add it in the
cmake/cp2k_utils.cmake file.

The cmake build system requires a minimum set of dependencies :

- a c, C++, and fortran compiler (gcc, intel oneapi, AMD or nvidia SDK, xlf, etc...)
- an MPI implementation
- DBCSR.
- openmp
- any flavor of BLAS, LAPACK, SCALAPACK.
- cmake of course

Major vendors implementations of BLAS, LAPACK, and scalapack are supported. The
build system was tested with MKL, cray libsci, openblas, flexiblas but it should
also work with blis, or ATLAS. Corresponding findPACKAGE.cmake are included but
they still need testing.

Options turned on by default are CP2K_USE_LIBXSMM, CP2K_USE_FFTW3,
CP2K_USE_LIBXC, CP2K_USE_COSMA, CP2K_USE_LIBINT2. Additionally MPI, DBCSR,
OPENMP, SCALAPACK, and BLAS/LAPACK are mandatory and can not be turned off. the
arguement `-DCP2K_USE_OPTION=ON, OFF` can be added to the cmake command line
turn `ON` or `OFF` a specific option. The list of currently supported optional
dependencies is

- CP2K_USE_SIRIUS = OFF : add SIRIUS support to cp2k

- CP2K_USE_FFTW3 = ON : add support of fftw3 (on by default)

- CP2K_USE_ELPA = OFF : add elpa support (off by default) WARNING : Expect the
  detection to fail at that stage

- CP2K_USE_PEXSI = OFF

- CP2K_USE_SUPERLU = OFF : detection should work but needs improvement

- CP2K_USE_COSMA = ON : Add cosma dropin replacement for sclapack pdgemnm

- CP2K_USE_LIBINT2 = ON : add libint2 support (detection works ok, module files
  may not be found at compilation time though)

- CP2K_USE_VORI = OFF : detection is fine compilation might fail at linking time
  (investigating why)

- CP2K_USE_QUIP = OFF

- CP2K_USE_SPGLIB = ON : everything alright

- CP2K_USE_LIBXC = ON : everything is fine, use pkgconfig by default (ideally
  the library should be built with cmake, if so we can get rid off the
  FindLibxc.cmake)

- CP2K_USE_SPLA = OFF : enable spla off-loading capabilities (use cmake modules
  to detect it)

- CP2K_USE_METIS = OFF :

- CP2K_USE_LIBXSMM = ON : use libxsmm library for small matrices operations.
  detection based on pkg-config

- CP2K_USE_ACCEL = NONE, CUDA, HIP : enable gpu support

- CP2K_BLAS_VENDOR = MKL, SCI, OpenBLAS, FlexiBLAS, Armpl, auto : default is
  auto. cmake will search for the most common blas / lapack implementations. If
  possible indicate which implementation you are using.

- CP2K_SCALAPACK_VENDOR - MKL, SCI, GENERIC : similar to the option previous
  option but for scalapack

- CP2K_BLAS_THREADING = sequential, openmp, etc... : leave the default value (or
  use it at your own peril)

- CP2K_BLAS_INTERFACE = 32 bits, 64 bits : size of the integers for the matrices
  and vectors sizes. default 32 bits

- CP2K_DEV_OPTIONS = OFF : enable developer options. the main purpose is for
  debugging

  - CP2K_USE_GRID_GPU = ON : turn on of gpu support for collocate integrate
  - CP2K_USE_PW_GPU = ON, turn on or off gpu fft support
  - CP2K_USE_DBM_GPU = ON turn on or off dbm gpu support

It is also possible to compile CP2K with GPU support namely CUDA or HIP. To do
so, add `-DCP2K_USE_ACCEL=CUDA,HIP -DCP2K_WITH_GPU=gpu_arch` to the cmake
command line.

While compiling CP2K with CUDA support should not pose problems (finding
libcublas and libcufft might fail though with the nvidia hpc sdk), we should
expect issues when compiling the hip support.

ROCM 5.0.x is known to have a bug in the cmake configuration files. It is
possible to go around this but at the expense of time. The build system was not
tested with ROCM 5.1.x but this version shows performance regression and should
be avoided. The Jiting capabilities of ROCM 5.2.x do not work properly which
affects DBCSR. It is highly recommended to update ROCM to the latest version to
avoid all these issues. CP2K can be built with ROCM 5.2.x but GPU support in
dbcsr should be tunred off otherwise a crash should be expected.

## Threading with blas and lapack

CP2K expect by default a single threaded version of blas and lapack. The option
`-DCP2K_BLAS_THREADING` can change this behavior. Be careful when tweaking this
specific option as many implementations of blas / lapack are easier threaded or
(exclusive) sequential but not both. I think the only exception to this is MKL.
Also note that CP2K dependencies will most likely have the same issue (COSMA
with cray-libsci for instance)

## typical examples of cmake use

The following list gives several examples of cmake command lines. Just add
`-DCP2K_USE_SIRIUS=ON` to add support of SIRIUS in cp2k

`shell cmake -DCP2K_INSTALL_PREFIX=/myprefix ..`

then

`shell make`

- MKL

the command line is

````shell cmake -DCP2K_INSTALL_PREFIX=/myprefix -DCP2K_BLAS_VENDOR=MKL
-DCP2K_SCALAPACK_VENDOR=MKL ..```

- Cray environments (with cray-libsci)

```shell
MPICC=cc MPICXX=CC cmake -DCP2K_INSTALL_PREFIX=/myprefix
-DCP2K_BLAS_VENDOR=SCI -DCP2K_SCALAPACK_VENDOR=SCI .. ```

## CUDA / HIP

Let us consider the case where openblas and netlib scalapack are installed
(openmpi or mpich)

```shell
cmake -DCP2K_INSTALL_PREFIX=/myprefix -DCP2K_BLAS_VENDOR=openblas
-DCP2K_SCALAPACK_VENDOR=GENERIC -DCP2K_USE_ACCEL=CUDA -DCP2K_WITH_GPU=A100 ..```

if HIP is needed than

```shell
cmake -DCP2K_INSTALL_PREFIX=/myprefix -DCP2K_BLAS_VENDOR=openblas
-DCP2K_SCALAPACK_VENDOR=GENERIC -DCP2K_USE_ACCEL=HIP -DCP2K_WITH_GPU=Mi250 ..```

## troubleshooting

This build system is relatevily stable and was tested on Cray, IBM, and redhat
like distributions. However it is not perfect and problems will show up, that's
why the two build systems will be available. We encourage the user to test the
build system just reporting the output of 'cmake ..' is already beneficial.

The best way to report these problems is to open an issue including the cmake
command line, error message, and operating systems.

What is known to fail sometimes

- Nvidia hpc sdk : The location of the cuda maths libraries has changed
recently. While CUDA support will be detected, the cuda maths libraries may not.

- HIP : CMAKE support of ROCM is still under development and is known to fail
from time to time. Update to ROCM 5.3.x or above to solve the issue.

- BLAS / LAPACK / SCALAPACK : use the options 'CP2K_BLAS_VENDOR' and
'CP2K_SCALPACK_VENDOR' if you know that 'MKL' or 'SCI' (cray libsci) are
present. '-DCP2k_BLAS_VENDOR=OpenBLAS' will also help cmake to find OpenBLAS if
it is used. Detecting the scalapack library might also fail if the user
environment is not properly set up.
````

# When provided to CMake using the -C argument, the command line must
# also have a -DTPL_DIR:PATH=[third party library path] argument to 
# identify the location of the installed raja, umpire, camp libraries.

set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-10.3.1/bin/g++" CACHE PATH "")
set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-10.3.1/bin/gcc" CACHE PATH "")
set(CMAKE_Fortran_COMPILER "/usr/tce/packages/gcc/gcc-10.3.1/bin/gfortran" CACHE PATH "")
set(ENABLE_MPI On CACHE BOOL "")
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")
set(ENABLE_OPENMP Off CACHE BOOL "")
set(ENABLE_HDF5 ON CACHE BOOL "")
set(ENABLE_CUDA Off CACHE BOOL "")
set(ENABLE_HIP Off CACHE BOOL "")
set(ENABLE_RAJA On CACHE BOOL "")
set(ENABLE_UMPIRE On CACHE BOOL "")
set(CMAKE_EXPORT_COMPILE_COMMANDS On CACHE BOOL "")
set(CMAKE_CXX_FLAGS "-O3 -std=c++20" CACHE STRING "")
set(BLT_CXX_STD "c++20" CACHE STRING "")
set(BLT_FORTRAN_FLAGS "" CACHE STRING "")
set(ENABLE_TESTS On CACHE BOOL "") 
set(ENABLE_SAMRAI_TESTS On CACHE BOOL "")

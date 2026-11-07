# When provided to CMake using the -C argument, the command line must
# also have a -DTPL_DIR:PATH=[third party library path] argument to 
# identify the location of the installed raja, umpire, camp libraries.

# ml cce/18.0.0
# ml rocm/6.1.2
# export MPICH_GPU_SUPPORT_ENABLED=1
# export LD_LIBRARY_PATH=/opt/cray/pe/lib64:${LD_LIBRARY_PATH}
# export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}
# export LD_LIBRARY_PATH=${CRAY_MPICH_ROOTDIR}/gtl/lib:${LD_LIBRARY_PATH}

set(ENABLE_MPI On CACHE BOOL "")
set(CMAKE_CXX_COMPILER "/usr/tce/packages/cce-tce/cce-18.0.0/bin/crayCC" CACHE PATH "")
set(CMAKE_C_COMPILER "/usr/tce/packages/cce-tce/cce-18.0.0/bin/craycc" CACHE PATH "")
set(CMAKE_Fortran_COMPILER "/usr/tce/packages/cce-tce/cce-18.0.0/bin/crayftn" CACHE PATH "")
#set(CMAKE_Fortran_COMPILER "/usr/tce/packages/rocmcc-tce/rocmcc-6.1.2/bin/amdflang" CACHE PATH "")
set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "")
set(ENABLE_OPENMP Off CACHE BOOL "")
set(ENABLE_HDF5 Off CACHE BOOL "")
set(ENABLE_CUDA Off CACHE BOOL "")
set(ENABLE_HIP On CACHE BOOL "")
set(ROCM_PATH "/opt/rocm-6.1.2" CACHE PATH "")
set(HIP_ROOT_DIR "/opt/rocm-6.1.2/hip" CACHE PATH "")
# set(HIP_HIPCC_FLAGS "-D__HIP_ROCclr__ -D__HIP_PLATFORM_AMD__ -DCAMP_USE_PLATFORM_DEFAULT_STREAM -D__HIP_ARCH_GFX90A__=1 --rocm-path=/opt/rocm-6.1.2 -std=c++20 -x hip --offload-arch=gfx942" CACHE STRING "")
set(CMAKE_HIP_ARCHITECTURES "gfx942" CACHE STRING "")
set(GPU_TARGETS "gfx942" CACHE STRING "")
set(AMDGPU_TARGETS "gfx942" CACHE STRING "")
set(ENABLE_RAJA On CACHE BOOL "")
set(ENABLE_UMPIRE On CACHE BOOL "")
set(CMAKE_EXPORT_COMPILE_COMMANDS On CACHE BOOL "")
set(CMAKE_CXX_FLAGS "-std=c++20" CACHE STRING "")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -DNDEBUG" CACHE STRING "")
set(BLT_CXX_STD "c++20" CACHE STRING "")
set(BLT_FORTRAN_FLAGS "" CACHE STRING "")
set(ENABLE_TESTS On CACHE BOOL "") 
set(ENABLE_SAMRAI_TESTS On CACHE BOOL "")
set(MIN_TEST_PROCS "2" CACHE STRING "")
set(MPI_CXX_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.30-cce-18.0.0/bin/mpicxx" CACHE PATH "")
set(MPI_C_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.30-cce-18.0.0/bin/mpicc" CACHE PATH "")
set(MPI_Fortran_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.30-cce-18.0.0/bin/mpifort" CACHE PATH "")

set(HIP_HIPCC_FLAGS "-D__HIP_ROCclr__ -D__HIP_PLATFORM_AMD__ -DCAMP_USE_PLATFORM_DEFAULT_STREAM -std=c++20 -x hip --offload-arch=gfx942" CACHE STRING "")

set(MPI_CXX_COMPILER_INCLUDE_DIRS "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.30-cce-18.0.0/include" CACHE STRING "")
set(MPI_CXX_LINK_FLAGS "-Wl,-rpath,/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.30-cce-18.0.0/ib -Wl,-rpath,/opt/cray/libfabric/2.1/lib64:/opt/cray/pe/pmi/6.1.14/lib:/opt/cray/pe/pals/1.2.12/lib:/opt/cray/pe/cce/18.0.0/cce/x86_64/lib:/opt/rocm-6.1.2/llvm/lib -L/opt/cray/pe/mpich/8.1.30/gtl/lib -lmpi_gtl_hsa" CACHE STRING "")
set(MPI_C_LINK_FLAGS "-Wl,-rpath,/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.30-cce-18.0.0/lib -Wl,-rpath,/opt/cray/libfabric/2.1/lib64:/opt/cray/pe/pmi/6.1.14/lib:/opt/cray/pe/pals/1.2.12/lib:/opt/cray/pe/cce/18.0.0/cce/x86_64/lib:/opt/rocm-6.1.2/llvm/lib" CACHE STRING "")
set(MPI_Fortran_LINK_FLAGS "-Wl,-rpath,/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.30-cce-18.0.0/lib -Wl,-rpath,/opt/cray/libfabric/2.1/lib64:/opt/cray/pe/pmi/6.1.14/lib:/opt/cray/pe/pals/1.2.12/lib:/opt/cray/pe/cce/18.0.0/cce/x86_64/lib:/opt/rocm-6.1.2/llvm/lib" CACHE STRING "")

set(MPI_mpi_cray_LIBRARY "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.30-cce-18.0.0/lib/libmpi_cray.so" CACHE PATH "")
set(MPI_mpifort_cray_LIBRARY "/usr/tce/packages/cray-mpich-tce/cray-mpich-8.1.30-cce-18.0.0/lib/libmpifort_cray.so" CACHE PATH "")

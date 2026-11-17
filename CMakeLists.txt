cmake_minimum_required(VERSION 3.8)

project(SAMRAI C CXX Fortran)

set(SAMRAI_VERSION_MAJOR 4)
set(SAMRAI_VERSION_MINOR 5)
set(SAMRAI_VERSION_PATCHLEVEL 0)
set(SAMRAI_VERSION
  "${SAMRAI_VERSION_MAJOR}.${SAMRAI_VERSION_MINOR}.${SAMRAI_VERSION_PATCHLEVEL}")


#------------------------------------------------------------------------------
# Initialize BLT build system
#------------------------------------------------------------------------------
if (DEFINED BLT_SOURCE_DIR)

    # Support a shared BLT outside of the repository if given a BLT_SOURCE_DIR
    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
        message(FATAL_ERROR "Cannot find SetupBLT.cmake in [${BLT_SOURCE_DIR}]")
    endif()

else()

    # Use internal 'blt' submodule path if BLT_SOURCE_DIR not provided
    set(BLT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/blt" CACHE PATH "")
    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
         message(FATAL_ERROR "The BLT submodule is not present. \
                  If in the git repository, please run the following commands:\n \
                  git submodule init\n \
                  git submodule update")
    endif()

endif()

# BLT Options need to be set in the CACHE
set(ENABLE_FORTRAN On CACHE Bool "Enable Fortran")
set(ENABLE_OPENMP Off CACHE Bool "Enable OpenMP")
set(ENABLE_CUDA Off CACHE Bool "Enable CUDA")
set(ENABLE_HIP Off CACHE Bool "Enable HIP")
set(ENABLE_MPI On CACHE Bool "Enable MPI.")
set(ENABLE_TESTS On CACHE Bool "Enable Tests.")
set(ENABLE_DOCS Off CACHE Bool "Enable Docs.")
set(BLT_CXX_STD "c++20" CACHE STRING "Version of C++ standard")
set(gtest_disable_pthreads ON CACHE BOOL "")

# Set allow list for C++ standard.  C++ 20 is the required minimum.
set(ALLOWED_CXX_STD "c++20")
list(FIND ALLOWED_CXX_STD "${BLT_CXX_STD}" FOUND_INDEX)
if (FOUND_INDEX EQUAL -1)
    message(FATAL_ERROR "BLT_CXX_STD is set to '${BLT_CXX_STD}', which is not supported.  The recommended standard is C++ 20.")
endif()

message(STATUS "Enable Test Flag Top: ${ENABLE_TESTS} " )

# Third party libraries
option(ENABLE_UMPIRE "Enable UMPIRE memory management" Off)
option(ENABLE_RAJA "Enable RAJA" Off)
option(ENABLE_CONDUIT "Enable CONDUIT" Off)
option(ENABLE_CALIPER "Enable CALIPER" Off)
option(ENABLE_HDF5 "Enable HDF5." On)
option(ENABLE_HYPRE "Enable HYPRE" Off)
option(ENABLE_PETSC "Enable PETSc" Off)
option(ENABLE_SILO "Enable Silo" Off)
option(ENABLE_SUNDIALS "Enable SUNDIALS" Off)

# SAMRAI options
option(ENABLE_SAMRAI_TESTS "Enable SAMRAI Test Programs" On)
option(ENABLE_PERF_TESTS "Enable Performance Tests." Off)
set(MIN_TEST_PROCS 1 CACHE INT "Minimum number of procs for tests.")
set(NUM_PERF_PROCS 8 CACHE INT "Number of processors for performance tests.")
option(ENABLE_CHECK_ASSERTIONS "Enable assertion checking." On)
option(ENABLE_BOX_COUNTING "Turns on box telemetry." Off)
option(ENABLE_DEPRECATED "Build with deprecated features." On)
option(ENABLE_TIMERS "Enable SAMRAI timers." On)
option(ENABLE_TOOLS "Turns on building of executable tools." Off)
option(DEBUG_INITIALIZE_UNDEFINED "Initialize new memory to undefined values" Off)
set(MAXDIM 3 CACHE INT "Maximum allowed dimension")
set(SAMRAI_MAX_COMM_BUFFERS 10 CACHE INT "Max number of chunks allowed to split large buffers")
option(ENABLE_NVTX_REGIONS "Enable NVTX regions with push/pop for NVVP." Off)
set(CUDA_ARCH "sm_70" CACHE STRING "Compute architecture to pass to CUDA builds")
set(CMAKE_CUDA_FLAGS "" CACHE STRING "")
set(CMAKE_INSTALL_LIBDIR lib)
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH Off CACHE Bool "Rpath uses Link path")
set(SAMRAI_RAJA_WORKGROUP_THREADS 512 CACHE INT "Number of workgroup threads")
set(ENABLE_SAMRAI_DEVICE_ALLOC CACHE BOOL "Use Device allocator on allocations for GPU")

include(GNUInstallDirs)

include(${BLT_SOURCE_DIR}/SetupBLT.cmake)

include_directories(${CMAKE_BINARY_DIR}/include)

include(cmake/CMakeBasics.cmake)

set(SAMRAI_LIBRARIES SAMRAI_algs SAMRAI_appu SAMRAI_geom SAMRAI_hier
    SAMRAI_math SAMRAI_mesh SAMRAI_pdat SAMRAI_solv SAMRAI_tbox
    SAMRAI_xfer)

add_subdirectory(source)
if (ENABLE_TOOLS)
  add_subdirectory(tools)
endif()

if (ENABLE_DOCS)
  add_subdirectory(docs)
endif()

export(TARGETS ${SAMRAI_LIBRARIES}
  FILE ${CMAKE_CURRENT_BINARY_DIR}/SAMRAITargets.cmake)

set(SAMRAI_MODULE "share/SAMRAI/cmake")
file(RELATIVE_PATH REL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}" "${CMAKE_INSTALL_PREFIX}/include")

set(CONF_INCLUDE_DIRS "${CMAKE_INSTALL_PREFIX}/${REL_INCLUDE_DIR}")
configure_file(SAMRAIConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/SAMRAIConfig.cmake" @ONLY)


configure_file(SAMRAIConfigVersion.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/SAMRAIConfigVersion.cmake" @ONLY)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/SAMRAIConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/SAMRAIConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${SAMRAI_MODULE}
  )

install(EXPORT SAMRAITargets DESTINATION
  ${CMAKE_INSTALL_PREFIX}/${SAMRAI_MODULE})

install(DIRECTORY
  ${CMAKE_CURRENT_BINARY_DIR}/include/SAMRAI
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

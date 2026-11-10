# Configure macros for calling fortran
include(FortranCInterface)

FortranCInterface_HEADER(
  ${CMAKE_BINARY_DIR}/include/SAMRAI/FC.h
  MACRO_NAMESPACE "CMAKE_FORTRAN_")

install(FILES ${CMAKE_BINARY_DIR}/include/SAMRAI/FC.h
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/SAMRAI)

include(CheckIncludeFiles)
include(CheckTypeSize)
include(CheckFunctionExists)
include(CheckCXXSourceCompiles)

if (${ENABLE_BOX_COUNTING})
  set (BOX_TELEMETRY On)
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set (DEBUG_INITIALIZE_UNDEFINED On)
  set (DEBUG_CHECK_ASSERTIONS On)
  set (DEBUG_CHECK_DIM_ASSERTIONS On)
endif()

#HAVE_CMATH
check_include_files("math.h" HAVE_CMATH)

#HAVE_CMATH_ISNAN
check_function_exists(std::isnan "cmath" HAVE_CMATH_ISNAN)

#HAVE_CTIME
check_include_files("time.h" HAVE_CTIME)

#HAVE_EXCEPTION_HANDLING

#HAVE_INLINE_ISNAND
check_function_exists(__inline_isnand "math.h" HAVE_INLINE_ISNAND)

#HAVE_ISNAN
check_cxx_source_compiles("int test = std::isnan(0.0)" HAVE_ISNAN)

#HAVE_ISNAND

#HAVE_MALLINFO
check_function_exists(mallinfo HAVE_MALLINFO)

#HAVE_MALLINFO2
check_function_exists(mallinfo2 HAVE_MALLINFO2)

#HAVE_MALLOC_H
check_include_files(malloc.h HAVE_MALLOC_H)

#SAMRAI_HAVE_SYS_TIMES_H
check_include_files(sys/times.h SAMRAI_HAVE_SYS_TIMES_H)

#SAMRAI_HAVE_UNISTD_H
check_include_files(unistd.h SAMRAI_HAVE_UNISTD_H)

#IOMANIP_HEADER_FILE
set(IOSTREAM_HEADER_FILE "<iomanip>")

#IOSTREAM_HEADER_FILE
set(IOSTREAM_HEADER_FILE "<iostream>")

#LACKS_SSTREAM
check_include_files(sstream HAVE_SSTREAM_H)

#LACKS_TEMPLATE_COMPLEX

#OPT_BUILD


#OSTRINGSTREAM_TYPE_IS_BROKEN

#OSTRSTREAM_TYPE_IS_BROKEN

#STL_SSTREAM_HEADER_FILE

if (${ENABLE_TIMERS})
  set(ENABLE_SAMRAI_TIMERS On)
endif ()

set(SAMRAI_MAXIMUM_DIMENSION ${MAXDIM})

set(SAMRAI_MAX_COMM_BUFFERS ${SAMRAI_MAX_COMM_BUFFERS})

configure_file(${PROJECT_SOURCE_DIR}/config/SAMRAI_config.h.cmake.in ${CMAKE_BINARY_DIR}/include/SAMRAI/SAMRAI_config.h)

install(FILES ${CMAKE_BINARY_DIR}/include/SAMRAI/SAMRAI_config.h
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/SAMRAI)


set (memory_utilities_sources
  main.cpp)

set (memory_utilities_depends SAMRAI_tbox)

if (ENABLE_OPENMP)
  set(memory_utilities_depends ${memory_utilities_depends} openmp)
endif ()

if (ENABLE_CUDA)
  set(memory_utilities_depends ${memory_utilities_depends} cuda)
endif ()

if (ENABLE_HIP)
  set(memory_utilities_depends ${memory_utilities_depends} blt::hip)
endif ()

blt_add_executable(
  NAME memory_utilities
  SOURCES ${memory_utilities_sources}
  DEPENDS_ON ${memory_utilities_depends})

set_target_properties(memory_utilities PROPERTIES LINKER_LANGUAGE CXX)

if(ENABLE_MPI)
  set(TASKS 1)
else()
  set(TASKS 0)
endif()

blt_add_test(
  NAME memory_utilities
  COMMAND memory_utilities
  NUM_MPI_TASKS ${TASKS})

if(ENABLE_MPI)
  blt_add_test(
    NAME memory_utilities_2procs
    COMMAND memory_utilities
    NUM_MPI_TASKS 2)
endif()

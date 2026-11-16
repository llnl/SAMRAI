set(boxutilities_depends SAMRAI_hier SAMRAI_tbox)

if (ENABLE_OPENMP)
  set(boxutilities_depends ${boxutilities_depends} openmp)
endif ()

if (ENABLE_CUDA)
  set(boxutilities_depends ${boxutilities_depends} cuda)
endif ()

if (ENABLE_HIP)
  set(boxutilities_depends ${boxutilities_depends} blt::hip)
endif ()

blt_add_executable(
  NAME boxutilities
  SOURCES main.cpp
  DEPENDS_ON ${boxutilities_depends})

set_target_properties(boxutilities PROPERTIES LINKER_LANGUAGE CXX)

if (ENABLE_MPI)
  set(TASKS 1)
else ()
  set(TASKS 0)
endif ()

blt_add_test(
  NAME boxutilities
  COMMAND boxutilities
  NUM_MPI_TASKS ${TASKS})

set (mblktree_sources
  main-mbtree.cpp)

set (mblktree_depends_on
    SAMRAI_hier
    SAMRAI_geom
    SAMRAI_mesh
    SAMRAI_tbox
    testlib)

if (ENABLE_OPENMP)
  set(mblktree_depends_on ${mblktree_depends_on} openmp)
endif ()

if (ENABLE_CUDA)
  set(mblktree_depends_on ${mblktree_depends_on} cuda)
endif ()

if (ENABLE_HIP)
  set(mblktree_depends_on ${mblktree_depends_on} blt::hip)
endif ()

blt_add_executable(
  NAME mblktree
  SOURCES ${mblktree_sources}
  DEPENDS_ON
    ${mblktree_depends_on})

target_compile_definitions(mblktree PUBLIC TESTING=1)

set_target_properties(mblktree PROPERTIES LINKER_LANGUAGE CXX)

file (GLOB test_inputs ${CMAKE_CURRENT_SOURCE_DIR}/test_inputs/*.input)
file (GLOB test_baselines ${CMAKE_CURRENT_SOURCE_DIR}/test_inputs/*baseline*)

execute_process(COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/test_inputs)
execute_process(COMMAND ln -s ${test_baselines} ${CMAKE_CURRENT_BINARY_DIR}/test_inputs)

samrai_add_tests(
  NAME mblktree
  EXECUTABLE mblktree
  INPUTS ${test_inputs}
  PARALLEL FALSE)

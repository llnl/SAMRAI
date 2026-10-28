set (testlib_headers
  AutoTester.h
  BaselineData.h
  DerivedVisOwnerData.h
  GeomUtilsAMR.h
  get-input-filename.h
  MblkHyperbolicLevelIntegrator.h
  MblkHyperbolicPatchStrategy.h
  MeshGenerationStrategy.h
  ShrunkenLevelGenerator.h
  SinusoidalFrontGenerator.h
  SkeletonBoundaryUtilities2.h
  SkeletonBoundaryUtilities3.h
  SphericalShellGenerator.h)

set (testlib_sources
  AutoTester.cpp
  BaselineData.cpp
  DerivedVisOwnerData.cpp
  MblkHyperbolicLevelIntegrator.cpp
  MblkHyperbolicPatchStrategy.cpp
  SkeletonBoundaryUtilities2.cpp
  SkeletonBoundaryUtilities3.cpp
  SkeletonCellDoubleConservativeLinearRefine.cpp
  SkeletonCellDoubleWeightedAverage.cpp
  MeshGenerationStrategy.cpp
  SinusoidalFrontGenerator.cpp
  SphericalShellGenerator.cpp
  ShrunkenLevelGenerator.cpp
  get-input-filename.cpp
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/skelbdry2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/skelbdry3d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/upfluxsum2d.f
  ${CMAKE_CURRENT_BINARY_DIR}/fortran/upfluxsum3d.f)

set (testlib_depends_on
  SAMRAI_algs
  SAMRAI_appu
  SAMRAI_hier
  SAMRAI_pdat
  SAMRAI_tbox
  )

# TODO CMake should resolve this dependency for us...
if (ENABLE_OPENMP)
  set(testlib_depends_on ${testlib_depends_on} openmp)
endif ()

if (ENABLE_CUDA)
  set(testlib_depends_on ${testlib_depends_on} cuda)
endif ()

if (ENABLE_HIP)
  set(testlib_depends_on ${testlib_depends_on} blt::hip)
endif ()

process_m4(NAME fortran/skelbdry2d)
process_m4(NAME fortran/skelbdry3d)
process_m4(NAME fortran/upfluxsum2d)
process_m4(NAME fortran/upfluxsum3d)

blt_add_library( NAME testlib
  SOURCES ${testlib_sources}
  HEADERS ${testlib_headers}
  DEPENDS_ON ${testlib_depends_on})

set_target_properties(testlib PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories( testlib
  PUBLIC ${PROJECT_SOURCE_DIR}/source/test/testlib)
# the scope of the following does not extend beyond target testlib
target_compile_definitions(testlib PRIVATE TESTING=1)

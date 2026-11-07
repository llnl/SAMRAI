set (hier_headers
  AssumedPartition.h
  AssumedPartitionBox.h
  BaseConnectorAlgorithm.h
  BaseGridGeometry.h
  BlockId.h
  BlueprintUtils.h
  BlueprintUtilsStrategy.h
  BoundaryBox.h
  BoundaryBoxUtils.h
  BoundaryLookupTable.h
  Box.h
  BoxContainer.h
  BoxContainerSingleBlockIterator.h
  BoxContainerSingleOwnerIterator.h
  BoxContainerUtils.h
  BoxGeometry.h
  BoxId.h
  BoxLevel.h
  BoxLevelConnectorUtils.h
  BoxLevelHandle.h
  BoxLevelStatistics.h
  BoxNeighborhoodCollection.h
  BoxOverlap.h
  BoxTree.h
  BoxUtilities.h
  CoarseFineBoundary.h
  CoarsenOperator.h
  ComponentSelector.h
  Connector.h
  ConnectorStatistics.h
  FlattenedHierarchy.h
  ForAll.h
  GlobalId.h
  HierarchyNeighbors.h
  Index.h
  IntVector.h
  LocalId.h
  MappingConnector.h
  MappingConnectorAlgorithm.h
  MultiblockBoxTree.h
  OverlapConnectorAlgorithm.h
  Patch.h
  PatchBoundaries.h
  PatchData.h
  PatchDataFactory.h
  PatchDataRestartManager.h
  PatchDescriptor.h
  PatchFactory.h
  PatchGeometry.h
  PatchHierarchy.h
  PatchLevel.h
  PatchLevelFactory.h
  PeriodicId.h
  PeriodicShiftCatalog.h
  PersistentOverlapConnectors.h
  ProcessorMapping.h
  RealBoxConstIterator.h
  RefineOperator.h
  SequentialLocalIdGenerator.h
  SingularityFinder.h
  TimeInterpolateOperator.h
  TransferOperatorRegistry.h
  Transformation.h
  UncoveredBoxIterator.h
  Variable.h
  VariableContext.h
  VariableDatabase.h)

set (hier_sources
  AssumedPartition.cpp
  AssumedPartitionBox.cpp
  BaseConnectorAlgorithm.cpp
  BaseGridGeometry.cpp
  BlockId.cpp
  BlueprintUtils.cpp
  BoundaryBox.cpp
  BoundaryBoxUtils.cpp
  BoundaryLookupTable.cpp
  Box.cpp
  BoxContainer.cpp
  BoxContainerSingleBlockIterator.cpp
  BoxContainerSingleOwnerIterator.cpp
  BoxContainerUtils.cpp
  BoxGeometry.cpp
  BoxId.cpp
  BoxLevel.cpp
  BoxLevelConnectorUtils.cpp
  BoxLevelHandle.cpp
  BoxLevelStatistics.cpp
  BoxNeighborhoodCollection.cpp
  BoxOverlap.cpp
  BoxTree.cpp
  BoxUtilities.cpp
  CoarseFineBoundary.cpp
  CoarsenOperator.cpp
  ComponentSelector.cpp
  Connector.cpp
  ConnectorStatistics.cpp
  FlattenedHierarchy.cpp
  GlobalId.cpp
  HierarchyNeighbors.cpp
  Index.cpp
  IntVector.cpp
  LocalId.cpp
  MappingConnector.cpp
  MappingConnectorAlgorithm.cpp
  MultiblockBoxTree.cpp
  OverlapConnectorAlgorithm.cpp
  Patch.cpp
  PatchBoundaries.cpp
  PatchData.cpp
  PatchDataFactory.cpp
  PatchDataRestartManager.cpp
  PatchDescriptor.cpp
  PatchFactory.cpp
  PatchGeometry.cpp
  PatchHierarchy.cpp
  PatchLevel.cpp
  PatchLevelFactory.cpp
  PeriodicId.cpp
  PeriodicShiftCatalog.cpp
  PersistentOverlapConnectors.cpp
  ProcessorMapping.cpp
  RealBoxConstIterator.cpp
  RefineOperator.cpp
  SingularityFinder.cpp
  TimeInterpolateOperator.cpp
  TransferOperatorRegistry.cpp
  Transformation.cpp
  UncoveredBoxIterator.cpp
  Variable.cpp
  VariableContext.cpp
  VariableDatabase.cpp)

set (hier_depends SAMRAI_tbox)

if (ENABLE_MPI)
  set (hier_depends ${hier_depends} mpi)
endif ()

if (ENABLE_CUDA)
  set (hier_depends ${hier_depends} cuda)
endif ()

if (ENABLE_HIP)
  set (hier_depends ${hier_depends} blt::hip)
endif ()

blt_add_library(
  NAME SAMRAI_hier
  SOURCES ${hier_sources}
  HEADERS ${hier_headers}
  DEPENDS_ON ${hier_depends})

target_include_directories(
  SAMRAI_hier
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>
  $<INSTALL_INTERFACE:include>)

install(TARGETS SAMRAI_hier
  EXPORT SAMRAITargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(FILES ${hier_headers}
  DESTINATION include/SAMRAI/hier)

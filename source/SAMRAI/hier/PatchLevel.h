/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A collection of patches at one level of the AMR hierarchy
 *
 ************************************************************************/

#ifndef included_hier_PatchLevel
#define included_hier_PatchLevel

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BoxContainerSingleBlockIterator.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/PatchFactory.h"
#include "SAMRAI/hier/ProcessorMapping.h"
#include "SAMRAI/tbox/StagedKernelFusers.h"
#include "SAMRAI/tbox/Utilities.h"

#include <map>
#include <vector>
#include <memory>

namespace SAMRAI {
namespace hier {

class BaseGridGeometry;

/*!
 * @brief Container class for patches defined at a single level of the
 * AMR hierarchy.
 *
 * The patches in a patch level are distributed across the processors
 * of a parallel machine, so not all patches reside on the local processor.
 * (However, each patch is assigned to one and only one processor.)
 *
 * To iterate over the local patches in a patch level, use the patch
 * level iterator class (PatchLevel::Iterator).
 *
 * @see BasePatchLevel
 * @see Patch
 * @see PatchDescriptor
 * @see PatchFactory
 * @see PatchLevelFactory
 * @see PatchLevel::Iterator
 */

class PatchLevel
{
public:
   /*!
    * @brief Default constructor.  PatchLevel must be initialized before it can
    * be used.
    */
   explicit PatchLevel(
      const tbox::Dimension& dim);

   /*!
    * @brief Construct a new patch level given a BoxLevel.
    *
    * This constructor makes a COPY of the supplied BoxLevel.  If the caller
    * intends to modify the supplied BoxLevel for other purposes after creating
    * this PatchLevel, then this constructor must be used rather than the
    * constructor taking a std::shared_ptr<BoxLevel>.
    *
    * The BoxLevel provides refinement ratio information, establishing
    * the ratio between the index space of the new level and some reference
    * level (typically level zero) in some patch hierarchy.
    *
    * The ratio information provided by the BoxLevel is also used
    * by the grid geometry instance to initialize geometry information
    * of both the level and the patches on that level.
    *
    * @param[in]  box_level
    * @param[in]  grid_geometry
    * @param[in]  descriptor The PatchDescriptor used to allocate patch data
    *             on the local processor
    * @param[in]  factory Optional PatchFactory.  If none specified, a default
    *             (standard) patch factory will be used.
    * @param[in]  defer_boundary_box_creation Flag to indicate suppressing
    *             construction of the boundary boxes.
    *
    * @pre grid_geometry
    * @pre descriptor
    * @pre box_level.getDim() == grid_geometry->getDim()
    * @pre box_level.getRefinementRatio() != IntVector::getZero(getDim())
    * @pre all components of box_level's refinement ratio must be nonzero and,
    *      all components not equal to 1 must have the same sign
    */
   PatchLevel(
      const BoxLevel& box_level,
      const std::shared_ptr<BaseGridGeometry>& grid_geometry,
      const std::shared_ptr<PatchDescriptor>& descriptor,
      const std::shared_ptr<PatchFactory>& factory =
         std::shared_ptr<PatchFactory>(),
      bool defer_boundary_box_creation = false);

   /*!
    * @brief Construct a new patch level given a BoxLevel.
    *
    * This constructor ACQUIRES the supplied BoxLevel.  If the caller will not
    * modify the supplied BoxLevel for other purposes after creating this
    * PatchLevel, then this constructor may be used rather than the constructor
    * taking a BoxLevel&.  Use of this constructor where permitted is more
    * efficient as it avoids copying an entire BoxLevel.  Note that this
    * constructor locks the supplied BoxLevel so that any attempt by the caller
    * to modify it after calling this constructor will result in an
    * unrecoverable error.
    *
    * The BoxLevel provides refinement ratio information, establishing
    * the ratio between the index space of the new level and some reference
    * level (typically level zero) in some patch hierarchy.
    *
    * The ratio information provided by the BoxLevel is also used
    * by the grid geometry instance to initialize geometry information
    * of both the level and the patches on that level.
    *
    * @param[in]  box_level
    * @param[in]  grid_geometry
    * @param[in]  descriptor The PatchDescriptor used to allocate patch data
    *             on the local processor
    * @param[in]  factory Optional PatchFactory.  If none specified, a default
    *             (standard) patch factory will be used.
    * @param[in]  defer_boundary_box_creation Flag to indicate suppressing
    *             construction of the boundary boxes.
    *
    * @pre grid_geometry
    * @pre descriptor
    * @pre box_level.getDim() == grid_geometry->getDim()
    * @pre box_level.getRefinementRatio() != IntVector::getZero(getDim())
    * @pre all components of box_level's refinement ratio must be nonzero and,
    *      all components not equal to 1 must have the same sign
    */
   PatchLevel(
      const std::shared_ptr<BoxLevel> box_level,
      const std::shared_ptr<BaseGridGeometry>& grid_geometry,
      const std::shared_ptr<PatchDescriptor>& descriptor,
      const std::shared_ptr<PatchFactory>& factory =
         std::shared_ptr<PatchFactory>(),
      bool defer_boundary_box_creation = false);

   /*!
    * @brief Construct a new patch level from the specified PatchLevel
    * restart database.
    *
    * The box, mapping, and ratio to level zero data which are normally
    * passed in during the construction of a new patch level are
    * retrieved from the specified restart database.
    *
    * @param[in]  restart_database
    * @param[in]  grid_geometry
    * @param[in]  descriptor The PatchDescriptor used to allocate patch
    *             data.
    * @param[in]  factory
    * @param[in]  defer_boundary_box_creation Flag to indicate suppressing
    *             construction of the boundary boxes.  @b Default: false
    *
    * @pre restart_database
    * @pre grid_geometry
    * @pre descriptor
    */
   PatchLevel(
      const std::shared_ptr<tbox::Database>& restart_database,
      const std::shared_ptr<BaseGridGeometry>& grid_geometry,
      const std::shared_ptr<PatchDescriptor>& descriptor,
      const std::shared_ptr<PatchFactory>& factory,
      bool defer_boundary_box_creation = false);

   /*!
    * @brief The virtual destructor for patch level deallocates all patches.
    */
   virtual ~PatchLevel() = default;

   /*!
    * @brief Get the level number
    *
    * @return the number of this level in a hierarchy, or the number of
    * a hierarchy level matching the index space of this level. If this
    * level does not align with the index space of a level in the hierarchy,
    * then this value is -1.  When the level is in a hierarchy, the return
    * value of the number of the level in the hierarchy.
    *
    * @see inHierarchy()
    */
   int
   getLevelNumber() const
   {
      return d_level_number;
   }

   /*!
    * @brief Set the number of this level to the level in the hierarchy
    * aligning with the index space of this level.
    *
    * The default value is -1 meaning the level index space does not align
    * with that of any hierarchy level.
    *
    * @param[in]  level
    */
   void
   setLevelNumber(
      const int level)
   {
      d_level_number = level;
      for (Iterator p(begin()); p != end(); ++p) {
         p->setPatchLevelNumber(d_level_number);
      }
   }

   /*!
    * @brief Convenience method to get the next coarser level
    * number in a hierarchy.
    *
    * Used for data interpolation from coarser levels.  If the
    * level is in a hierarchy, then this value is getLevelNumber() - 1.
    *
    * @see inHierarchy()
    *
    * @return The next coarser level in the hierarchy or -1 if the level
    * does not exist in the hierarchy.
    */
   int
   getNextCoarserHierarchyLevelNumber() const
   {
      return d_next_coarser_level_number;
   }

   /*!
    * @brief Convenience method to set the number of of the next coarser
    * level in a hierarchy.
    *
    * For the purposes of data interpolation from coarser levels, set the
    * next coarser level in a hierarchy.  The default of -1 means the
    * level does not relate to any hierarchy.
    *
    * @param[in]  level
    */
   void
   setNextCoarserHierarchyLevelNumber(
      const int level)
   {
      d_next_coarser_level_number = level;
   }

   /*!
    * @brief Determine if this level resides in a hierarchy.
    *
    * @return true if this level resides in a hierarchy, otherwise false.
    */
   bool
   inHierarchy() const
   {
      return d_in_hierarchy;
   }

   /*!
    * @brief Setting to indicate whether this level resides in a hierarchy.
    *
    * @param[in]  in_hierarchy Flag to indicate whether this level resides
    *             in a hierarchy.  @b Default: false
    */
   void
   setLevelInHierarchy(
      bool in_hierarchy)
   {
      d_in_hierarchy = in_hierarchy;
      for (Iterator p(begin()); p != end(); ++p) {
         p->setPatchInHierarchy(d_in_hierarchy);
      }
   }

   /*!
    * @brief Get the number of patches.
    *
    * This is equivalent to calling PatchLevel::getGlobalNumberOfPatches().
    */
   int
   getNumberOfPatches() const
   {
      return getGlobalNumberOfPatches();
   }

   /*!
    * @brief Get the local number of patches
    */
   int
   getLocalNumberOfPatches() const
   {
      return static_cast<int>(d_box_level->getLocalNumberOfBoxes());
   }

   /*!
    * @brief Get the global number of patches.
    */
   int
   getGlobalNumberOfPatches() const
   {
      return d_box_level->getGlobalNumberOfBoxes();
   }

   /*!
    * @brief Get the local number of Cells.
    */
   int
   getLocalNumberOfCells() const
   {
      return static_cast<int>(d_box_level->getLocalNumberOfCells());
   }

   /*!
    * @brief Get the global number of cells
    */
   long int
   getGlobalNumberOfCells() const
   {
      return d_box_level->getGlobalNumberOfCells();
   }

   /*!
    * @brief Get a Patch based on its GlobalId.
    *
    * @param[in]  gid
    *
    * @return A std::shared_ptr to the Patch indicated by the GlobalId.
    */
   const std::shared_ptr<Patch>&
   getPatch(
      const GlobalId& gid) const
   {
      BoxId mbid(gid);
      PatchContainer::const_iterator it = d_patches.find(mbid);
      if (it == d_patches.end()) {
         TBOX_ERROR("PatchLevel::getPatch error: GlobalId "
            << gid << " does not exist locally.\n"
            << "You must specify the GlobalId of a current local patch.");
      }
      return it->second;
   }

   /*!
    * @brief Get a Patch based on its BoxId.
    *
    * @param[in]  mbid
    *
    * @return A std::shared_ptr to the Patch indicated by the BoxId.
    *
    * @pre d_patches.find(mbid) != d_patches.end()
    */
   std::shared_ptr<Patch>
   getPatch(
      const BoxId& mbid) const
   {
      const PatchContainer::const_iterator mi = d_patches.find(mbid);
      if (mi == d_patches.end()) {
         TBOX_ERROR("PatchLevel::getPatch error: BoxId "
            << mbid << " does not exist locally.\n"
            << "You must specify the BoxId of a current local box"
            << " that is not a periodic image.");
      }
      return (*mi).second;
   }

   /*!
    * @brief Get a patch using a random access index.
    *
    * The index specifies the position of the patch as would be
    * encountered when iterating through the patches.
    */
   const std::shared_ptr<Patch>& getPatch(size_t index) const
   {
      if (index >= d_patch_vector.size()) {
         TBOX_ERROR("PatchLevel::getPatch error: index "
            << index << " is too big.\n"
            << "There are only " << d_patch_vector.size() << " patches.");
      }
      return d_patch_vector[index];
   }

   /*!
    * @brief Get the PatchDescriptor
    *
    * @return pointer to the patch descriptor for the hierarchy.
    */
   std::shared_ptr<PatchDescriptor>
   getPatchDescriptor() const
   {
      return d_descriptor;
   }

   /*!
    * @brief Get the PatchFactory
    *
    * @return the factory object used to created patches in the level.
    */
   std::shared_ptr<PatchFactory>
   getPatchFactory() const
   {
      return d_factory;
   }

   /*!
    * @brief Get the grid geometry
    *
    * @return A std::shared_ptr to the grid geometry description.
    */
   std::shared_ptr<BaseGridGeometry>
   getGridGeometry() const
   {
      return d_geometry;
   }

   /*!
    * @brief Update this patch level through refining.
    *
    * The data members of this patch level are updated by refining the
    * information on a given coarse level using the given ratio between
    * the two levels.  The fine level will cover the same physical space as
    * the coarse level and will have the same number of patches with the
    * same mapping of those patches to processors.  However, the index
    * space of the level will be refined by the specified ratio.
    *
    * @par Assumptions
    * If the fine grid geometry is null (default case), then it is assumed
    * that this level is to use the same grid geometry as the given coarse
    * level and the ratio to level zero is set relative to the given coarse
    * level.  Otherwise, we use the given grid geometry (assumed to be a proper
    * refinement of the grid geometry used on the given coarse level) and copy
    * ratio to level zero from given coarse level.  In other words, the function
    * can be used to produce two different results.
    *
    * <ol>
    *   <li> When passed a null grid geometry pointer, the refined patch level
    *        can be used for data exchange operations with the AMR hierarchy
    *        in which the coarse level resides -- both levels are defined with
    *        respect to the index space of the grid geometry object which they
    *        share.  Thus, the refined patch level can be used in data
    *        exchanges with the AMR hierarchy of the coarse level
    *        automatically.
    *   <li> Second, when passed a non-null fine grid geometry pointer, the
    *        level is defined relative to that geometry and the refined patch
    *        level cannot be used in data exchanges with the AMR hierarchy
    *        of the coarse level automatically in general.  This mode is
    *        used to construct a refined copy of an entire patch hierarchy,
    *        typically.
    * </ol>
    *
    * @param[in]  coarse_level
    * @param[in]  refine_ratio
    * @param[in]  fine_grid_geometry @b Default: std::shared_ptr to a null
    *             grid geometry
    * @param[in]  defer_boundary_box_creation @b Default: false
    *
    * @pre coarse_level
    * @pre refine_ratio > IntVector::getZero(getDim())
    * @pre (getDim() == coarse_level->getDim()) &&
    *      (getDim() == refine_ratio.getDim())
    * @pre !fine_grid_geometry || getDim() == fine_grid_geometry->getDim()
    */
   void
   setRefinedPatchLevel(
      const std::shared_ptr<PatchLevel>& coarse_level,
      const IntVector& refine_ratio,
      const std::shared_ptr<BaseGridGeometry>& fine_grid_geometry =
         std::shared_ptr<BaseGridGeometry>(),
      bool defer_boundary_box_creation = false);

   /*!
    * @brief Update this patch through coarsening.
    *
    * The data members of this patch level are updated by coarsening the
    * information on a given fine level using the given ratio between
    * the two levels.  The coarse level will cover the same physical space as
    * the fine level and will have the patches with the same
    * GlobalIndices.  However, the index space of the level will be coarsened
    * by the specified ratio.
    * @par Assumptions
    * If the coarse grid geometry is null (default case), then it is assumed
    * that this level is to use the same grid geometry as the given fine
    * level and the ratio to level zero is set relative to the given fine
    * level.  Otherwise, we use the given grid geometry (assumed to be a proper
    * coarsening of the grid geometry used on the given fine level) and copy
    * ratio to level zero from given fine level.  In other words, the function
    * can be used to produce two different results.
    *
    * <ol>
    *   <li> When passed a null grid geometry pointer, the coarsened
    *        patch level can be used for data exchange operations with the
    *        AMR hierarchy in which the fine level resides -- both levels
    *        are defined with respect to the index space of the grid geometry
    *        object which they share.  Thus, the coarsened patch level can be
    *        used in data exchanges with the AMR hierarchy of the fine level
    *        automatically.
    *   <li> When passed a non-null coarse grid geometry pointer, the level is
    *        defined relative to that geometry and the coarsened patch level
    *        cannot be used in data exchanges with the AMR hierarchy of the
    *        fine level automatically in general.  This mode is used to
    *        construct a coarsened copy of an entire patch hierarchy,
    *        typically.
    * </ol>
    *
    * @param[in]  fine_level
    * @param[in]  coarsen_ratio
    * @param[in]  coarse_grid_geom @b Default: std::shared_ptr to a null
    *             grid geometry
    * @param[in]  defer_boundary_box_creation @b Default: false
    *
    * @pre fine_level
    * @pre coarsen_ratio > IntVector::getZero(getDim())
    * @pre (getDim() == fine_level->getDim()) &&
    *      (getDim() == coarsen_ratio.getDim())
    * @pre !coarse_grid_geom || getDim() == coarse_grid_geom->getDim()
    */
   void
   setCoarsenedPatchLevel(
      const std::shared_ptr<PatchLevel>& fine_level,
      const IntVector& coarsen_ratio,
      const std::shared_ptr<BaseGridGeometry>& coarse_grid_geom =
         std::shared_ptr<BaseGridGeometry>(),
      bool defer_boundary_box_creation = false);

   /*!
    * @brief Create and store the boundary boxes for this level.
    *
    * If boundary boxes have already been constructed, this function
    * does nothing.
    * @note
    * If the level is constructed with boundary box creation deferred,
    * this method must be called before any attempt at filling data at
    * physical boundaries.  This function is called from
    * xfer::RefineSchedule prior to any physical boundary operations.
    */
   void
   setBoundaryBoxes()
   {
      if (!d_boundary_boxes_created) {
         d_geometry->setBoundaryBoxes(*this);
         d_boundary_boxes_created = true;
      }
   }

   /*!
    * @brief Get the physical domain.
    *
    * @return A const reference to the box array that defines
    * the extent of the index space on the level.
    */
   const std::vector<BoxContainer>&
   getPhysicalDomainArray() const
   {
      return d_physical_domain;
   }

   const BoxContainer&
   getPhysicalDomain(
      const BlockId& block_id) const
   {
      return d_physical_domain[block_id.getBlockValue()];
   }

   /*!
    * @brief Get the box defining the patches on the level.
    *
    * The internal state of PatchLevel (where boxes are concern) is
    * dependent on the BoxLevel associated with it, and computed
    * only if getBoxes() is called.  The first call to getBoxes() must be
    * done by all processors as it requires communication.
    *
    * @return a const reference to the box array that defines
    * the patches on the level.
    */
   const BoxContainer&
   getBoxes() const
   {
      if (!d_has_globalized_data) {
         initializeGlobalizedBoxLevel();
      }
      return d_boxes;
   }

   /*!
    * @brief Get boxes for a particular block.
    *
    * The boxes from only the specified block will be stored in the
    * output BoxContainer.
    *
    * @param[out] boxes
    * @param[in] block_id
    */
   void
   getBoxes(
      BoxContainer& boxes,
      const BlockId& block_id) const;

   /*!
    * @brief Get the BoxLevel associated with the PatchLevel.
    *
    * @return a reference to a std::shared_ptr to the BoxLevel
    * associated with the PatchLevel.
    */
   const std::shared_ptr<BoxLevel>&
   getBoxLevel() const
   {
      return d_box_level;
   }

   /*!
    * @brief Get the globalized version of the BoxLevel associated
    * with the PatchLevel.
    * @note
    * The first time this method is used, a global communication is
    * done.  Thus all processors must use this method the first time
    * any processor uses it.
    *
    * @return The globalized version of the BoxLevel associated
    * with the PatchLevel.
    */
   const BoxLevel&
   getGlobalizedBoxLevel() const
   {
      if (!d_has_globalized_data) {
         initializeGlobalizedBoxLevel();
      }
      return d_box_level->getGlobalizedVersion();
   }

   /*!
    * @brief Get the mapping of patches to processors.
    *
    * @return A const reference to the mapping of patches to processors.
    */
   const ProcessorMapping&
   getProcessorMapping() const
   {
      if (!d_has_globalized_data) {
         initializeGlobalizedBoxLevel();
      }
      return d_mapping;
   }

   /*!
    * @brief Get the ratio between the index space of this PatchLevel and
    * the reference level in the AMR hierarchy.
    *
    * @return A const reference to the vector ratio between the index
    * space of this patch level and that of a reference level in AMR
    * hierarchy (that is, level zero).
    */
   const IntVector&
   getRatioToLevelZero() const
   {
      return d_ratio_to_level_zero;
   }

   /*!
    * @brief Get the ratio between this level and the next coarser
    * level in the patch hierarchy.
    *
    * This vector is set with the setRatioToCoarserLevel() function.
    * If the level is not in a hierarchy, a default ratio of zero is returned.
    *
    * @return the vector ratio between this level and the next coarser
    * level in the patch hierarchy.
    */
   const IntVector&
   getRatioToCoarserLevel() const
   {
      return d_ratio_to_coarser_level;
   }

   /*!
    * @brief Set the ratio between this level and the next coarser
    * level in the patch hierarchy.
    *
    * This is required only when level resides in a hierarchy.
    *
    * @param[in] ratio
    */
   void
   setRatioToCoarserLevel(
      const IntVector& ratio)
   {
      d_ratio_to_coarser_level = ratio;
   }

   /*!
    * @brief Get the processor mapping for the patch.
    *
    * @return the processor that owns the specified patch.  The patches
    * are numbered starting at zero.
    *
    * @param[in] box_id Patch's BoxId
    */
   int
   getMappingForPatch(
      const BoxId& box_id) const
   {
      // Note: p is required to be a local index.
      /*
       * This must be for backward compatability, because if p is a local
       * index, the mapping is always to d_box_level->getRank().
       * Here is the old code:
       *
       * return d_box_level->getBoxStrict(p)->getOwnerRank();
       */
      NULL_USE(box_id);
      return d_box_level->getMPI().getRank();
   }

   /*!
    * @brief Get the box for the specified patch
    *
    * @return The box for the specified patch.
    *
    * @param[in] box_id Patch's BoxId
    *
    * @pre box_id.getOwnerRank() == getBoxLevel()->getMPI().getRank()
    */
   const Box&
   getBoxForPatch(
      const BoxId& box_id) const
   {
      TBOX_ASSERT(box_id.getOwnerRank() == d_box_level->getMPI().getRank());
      return getPatch(box_id)->getBox();
   }

   /*!
    * @brief Determine if the patch is adjacent to a non-periodic
    * physical domain boundary.
    *
    * @param[in] box_id Patch's BoxId
    *
    * @return True if patch with given number is adjacent to a non-periodic
    * physical domain boundary.  Otherwise, false.
    *
    * @pre box_id.getOwnerRank() == getBoxLevel()->getMPI().getRank()
    */
   bool
   patchTouchesRegularBoundary(
      const BoxId& box_id) const
   {
      TBOX_ASSERT(box_id.getOwnerRank() == d_box_level->getMPI().getRank());
      return getPatch(box_id)->getPatchGeometry()->getTouchesRegularBoundary();
   }

   /*!
    * @brief Allocate the specified component on all patches.
    *
    * @param[in]  id
    * @param[in]  timestamp @b Default: zero (0.0)
    */
   void
   allocatePatchData(
      const int id,
      const double timestamp = 0.0)
   {
      for (Iterator ip(begin()); ip != end(); ++ip) {
         ip->allocatePatchData(id, timestamp);
      }

#if defined(HAVE_RAJA)
      tbox::StagedKernelFusers::getInstance()->launchAndCleanup();
#endif

   }

   /*!
    * @brief Allocate the specified components on all patches.
    *
    * @param[in]  components The componentSelector indicating
    *             which elements to allocate
    * @param[in]  timestamp @b Default: zero (0.0)
    */
   void
   allocatePatchData(
      const ComponentSelector& components,
      const double timestamp = 0.0)
   {
      for (Iterator ip(begin()); ip != end(); ++ip) {
         ip->allocatePatchData(components, timestamp);
      }

#if defined(HAVE_RAJA)
      tbox::StagedKernelFusers::getInstance()->launchAndCleanup();
#endif

   }

   /*!
    * @brief Determine if the patch data has been allocated.
    *
    * @return True if (1) there are no patches in this patch level or  if
    *         (2) all of the patches have allocated the patch data component,
    *         otherwise false.
    *
    * @param[in] id The patch data identifier.
    */
   bool
   checkAllocated(
      const int id) const
   {
      bool allocated = true;
      for (PatchContainer::const_iterator mi = d_patches.begin();
           mi != d_patches.end(); ++mi) {
         allocated &= (*mi).second->checkAllocated(id);
      }
      return allocated;
   }

   /*!
    * @brief  Deallocate the specified component on all patches.
    *
    * This component will need to be reallocated before its next use.
    *
    * @param[in]  id The patch data identifier
    */
   void
   deallocatePatchData(
      const int id)
   {
      for (Iterator ip(begin()); ip != end(); ++ip) {
         ip->deallocatePatchData(id);
      }

#if defined(HAVE_RAJA)
      tbox::StagedKernelFusers::getInstance()->launchAndCleanup();
#endif
   }

   /*!
    * @brief Deallocate the specified components on all patches.
    *
    * Components will need to be reallocated before their next use.
    *
    * @param[in]  components The ComponentSelector indicating which
    *             components to deallocate.
    */
   void
   deallocatePatchData(
      const ComponentSelector& components)
   {
      for (Iterator ip(begin()); ip != end(); ++ip) {
         ip->deallocatePatchData(components);
      }

#if defined(HAVE_RAJA)
      tbox::StagedKernelFusers::getInstance()->launchAndCleanup();
#endif
   }

   /*!
    * @brief Get the dimension of this object.
    *
    * @return the dimension of this object.
    */
   const tbox::Dimension&
   getDim() const
   {
      return d_dim;
   }

   /*!
    * @brief Set the simulation time for the specified patch component.
    *
    * @param[in]  timestamp
    * @param[in]  id The patch identifier
    */
   void
   setTime(
      const double timestamp,
      const int id)
   {
      for (Iterator ip(begin()); ip != end(); ++ip) {
         ip->setTime(timestamp, id);
      }
   }

   /*!
    * @brief Set the simulation time for the specified patch components.
    *
    * @param[in] timestamp
    * @param[in] components The ComponentSelector indicating on which
    *            components to set the simulation time.
    */
   void
   setTime(
      const double timestamp,
      const ComponentSelector& components)
   {
      for (Iterator ip(begin()); ip != end(); ++ip) {
         ip->setTime(timestamp, components);
      }
   }

   /*!
    * @brief Set the simulation time for all allocated patch components.
    *
    * @param[in]  timestamp
    */
   void
   setTime(
      const double timestamp)
   {
      for (Iterator ip(begin()); ip != end(); ++ip) {
         ip->setTime(timestamp);
      }
   }

   /*!
    * @brief Find an overlap Connector with the given PatchLevel's BoxLevel as
    * its head and minimum Connector width.  If the specified Connector is not
    * found, take the specified action.
    *
    * If multiple Connectors fit the criteria, the one with the
    * smallest ghost cell width (based on the algebraic sum of the
    * components) is selected.
    *
    * @param[in] head Find the overlap Connector with this PatchLevel's
    *      BoxLevel as the head.
    * @param[in] min_connector_width Find the overlap Connector satisfying
    *      this minimum Connector width.
    * @param[in] not_found_action Action to take if Connector is not found.
    * @param[in] exact_width_only If true, the returned Connector will
    *      have exactly the requested connector width. If only a Connector
    *      with a greater width is found, a connector of the requested width
    *      will be generated.
    *
    * @return The Connector which matches the search criterion.
    *
    * @pre getBoxLevel()->isInitialized()
    * @pre head.getBoxLevel()->isInitialized()
    */
   const Connector&
   findConnector(
      const PatchLevel& head,
      const IntVector& min_connector_width,
      ConnectorNotFoundAction not_found_action,
      bool exact_width_only = true) const
   {
      return getBoxLevel()->findConnector(*head.getBoxLevel(),
         min_connector_width,
         not_found_action,
         exact_width_only);
   }

   /*!
    * @brief Find an overlap Connector with its transpose with the given
    * PatchLevel's BoxLevel as its head and minimum Connector widths.  If the
    * specified Connector is not found, take the specified action.
    *
    * If multiple Connectors fit the criteria, the one with the
    * smallest ghost cell width (based on the algebraic sum of the
    * components) is selected.
    *
    * @param[in] head Find the overlap Connector with this PatchLevel's
    *      BoxLevel as the head.
    * @param[in] min_connector_width Find the overlap Connector satisfying
    *      this minimum Connector width.
    * @param[in] transpose_min_connector_width Find the transpose overlap
    *      Connector satisfying this minimum Connector width.
    * @param[in] not_found_action Action to take if Connector is not found.
    * @param[in] exact_width_only If true, the returned Connector will
    *      have exactly the requested connector width. If only a Connector
    *      with a greater width is found, a connector of the requested width
    *      will be generated.
    * @return The Connector which matches the search criterion.
    *
    * @pre getBoxLevel()->isInitialized()
    * @pre head.getBoxLevel()->isInitialized()
    */
   const Connector&
   findConnectorWithTranspose(
      const PatchLevel& head,
      const IntVector& min_connector_width,
      const IntVector& transpose_min_connector_width,
      ConnectorNotFoundAction not_found_action,
      bool exact_width_only = true) const
   {
      return getBoxLevel()->findConnectorWithTranspose(*head.getBoxLevel(),
         min_connector_width,
         transpose_min_connector_width,
         not_found_action,
         exact_width_only);
   }

   /*!
    * @brief Create an overlap Connector, computing relationships by
    * globalizing data.
    *
    * The base will be this PatchLevel's BoxLevel.
    * Find Connector relationships using a (non-scalable) global search.
    *
    * @see Connector
    * @see Connector::initialize()
    *
    * @param[in] head This PatchLevel's BoxLevel will be the head.
    * @param[in] connector_width
    *
    * @return A const reference to the newly created overlap Connector.
    *
    * @pre getBoxLevel()->isInitialized()
    * @pre head.getBoxLevel()->isInitialized()
    */
   const Connector&
   createConnector(
      const PatchLevel& head,
      const IntVector& connector_width) const
   {
      return getBoxLevel()->createConnector(*head.getBoxLevel(),
         connector_width);
   }

   /*!
    * @brief Create an overlap Connector with its transpose, computing
    * relationships by globalizing data.
    *
    * The base will be this PatchLevel's BoxLevel.
    * Find Connector relationships using a (non-scalable) global search.
    *
    * @see Connector
    * @see Connector::initialize()
    *
    * @param[in] head This PatchLevel's BoxLevel will be the head.
    * @param[in] connector_width
    * @param[in] transpose_connector_width
    *
    * @return A const reference to the newly created overlap Connector.
    *
    * @pre getBoxLevel()->isInitialized()
    * @pre head.getBoxLevel()->isInitialized()
    */
   const Connector&
   createConnectorWithTranspose(
      const PatchLevel& head,
      const IntVector& connector_width,
      const IntVector& transpose_connector_width) const
   {
      return getBoxLevel()->createConnectorWithTranspose(*head.getBoxLevel(),
         connector_width,
         transpose_connector_width);
   }

   /*!
    * @brief Cache the supplied overlap Connector and its transpose
    * if it exists.
    *
    * @param[in] connector
    *
    * @pre connector
    * @pre getBoxLevel()->isInitialized()
    * @pre getBoxLevel() == connector->getBase()
    */
   void
   cacheConnector(
      std::shared_ptr<Connector>& connector) const
   {
      return getBoxLevel()->cacheConnector(connector);
   }

   /*!
    * @brief Returns whether the object has overlap Connectors with the given
    * PatchLevel's BoxLevel as the head and minimum Connector width.
    *
    * TODO:  does the following comment mean that this must be called
    * before the call to findConnector?
    *
    * If this returns true, the Connector fitting the specification
    * exists and findConnector() will not throw an assertion.
    *
    * @param[in] head Find the overlap Connector with this PatchLevel's
    *      BoxLevel as the head.
    * @param[in] min_connector_width Find the overlap Connector satisfying
    *      this minimum ghost cell width.
    *
    * @return True if a Connector is found, otherwise false.
    */
   bool
   hasConnector(
      const PatchLevel& head,
      const IntVector& min_connector_width) const
   {
      return getBoxLevel()->hasConnector(*head.getBoxLevel(),
         min_connector_width);
   }

   /*!
    * @brief Use the PatchLevel restart database to set the state of the
    * PatchLevel and to create all patches on the local processor.
    *
    * @par Assertions
    * Assertions will check that database is a non-null std::shared_ptr,
    * that the data being retrieved from the database are of
    * the type expected.  Also checked is the number of patches is positive,
    * and the number of patches and size of processor mapping array are the
    * same, and that the number of patches and the number of boxes on the
    * level are equal.
    *
    * @param[in,out] restart_db
    *
    * @pre restart_db
    */
   void
   getFromRestart(
      const std::shared_ptr<tbox::Database>& restart_db);

   /*!
    * @brief Write data to the restart database.
    *
    * Writes the data from the PatchLevel to the restart database.
    * Also tells all local patches to write out their state to
    * the restart database.
    *
    * @param[in,out]  restart_db
    *
    * @pre restart_db
    */
   void
   putToRestart(
      const std::shared_ptr<tbox::Database>& restart_db) const;

   /*!
    * @brief Print a patch level to varying details.
    *
    * If depth>0, print function will be called for each patch in the level.
    *
    * @param[in]  os The std::ostream in which to print to
    * @param[in]  border @b Default: empty string
    * @param[in]  depth @b Default: zero (0).
    *
    * @return 0.  Always.
    */
   int
   recursivePrint(
      std::ostream& os,
      const std::string& border = std::string(),
      int depth = 0);

private:
   /*
    * Static integer constant describing class's version number.
    */
   static const int HIER_PATCH_LEVEL_VERSION;

   /*
    * @brief Container of distributed patches on level.
    */
   typedef std::map<BoxId, std::shared_ptr<Patch> > PatchContainer;

   /*
    * @brief Vector of local patches on level.
    */
   typedef std::vector<std::shared_ptr<Patch> > PatchVector;

public:
   /*!
    * @brief Iterator for looping through local patches.
    */
   class Iterator
   {
      friend class PatchLevel;

public:
      /*!
       * @brief Copy constructor.
       *
       * @param[in]  other
       */
      Iterator(
         const Iterator& other);

      /*!
       * @brief Assignment operator
       */
      Iterator&
      operator = (
         const Iterator& rhs)
      {
         d_iterator = rhs.d_iterator;
         d_patches = rhs.d_patches;
         return *this;
      }

      /*!
       * @brief Dereference operator.
       */
      const std::shared_ptr<Patch>&
      operator * () const
      {
         return d_iterator->second;
      }

      /*!
       * @brief Delegation operations to the Patch pointer.
       */
      const std::shared_ptr<Patch>&
      operator -> () const
      {
         return d_iterator->second;
      }

      /*!
       * @brief Equality comparison.
       */
      bool
      operator == (
         const Iterator& rhs) const
      {
         return d_iterator == rhs.d_iterator;
      }

      /*!
       * @brief Inequality operator.
       */
      bool
      operator != (
         const Iterator& rhs) const
      {
         return d_iterator != rhs.d_iterator;
      }

      /*!
       * @brief Pre-increment.
       */
      Iterator&
      operator ++ ()
      {
         ++d_iterator;
         return *this;
      }

      /*!
       * @brief Post-increment.
       */
      Iterator
      operator ++ (
         int)
      {
         Iterator tmp_iterator = *this;
         ++d_iterator;
         return tmp_iterator;
      }

private:
      /*
       * Unimplemented default constructor.
       */
      Iterator();

      /*!
       * @brief Construct from a PatchLevel.
       *
       * @param[in]  patch_level
       * @param[in]  begin
       */
      Iterator(
         const PatchLevel* patch_level,
         bool begin);

      /*!
       * @brief The real iterator (this class is basically a wrapper).
       */
      PatchContainer::const_iterator d_iterator;

      /*!
       * @brief For supporting backward-compatible interface.
       */
      const PatchContainer* d_patches;

   };

   /*!
    * @brief Construct an iterator pointing to the first Patch in the
    * PatchLevel.
    */
   Iterator
   begin() const
   {
      return Iterator(this, true);
   }

   /*!
    * @brief Construct an iterator pointing to the last Patch in the
    * PatchLevel.
    */
   Iterator
   end() const
   {
      return Iterator(this, false);
   }

   /*!
    * @brief Typdef PatchLevel::Iterator to standard iterator nomenclature.
    */
   typedef Iterator iterator;

private:
   /**
    * @brief Static initialization to be done at startup.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback();

   /**
    * @brief Static cleanup to be done at shutdown.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

   void
   initializeGlobalizedBoxLevel() const;

   /*!
    * @brief Dimension of the object
    */
   const tbox::Dimension d_dim;

   /*!
    * @brief Number of blocks that can be represented by this level.
    */
   size_t d_number_blocks;

   /*!
    * Primary metadata describing the PatchLevel.
    */
   std::shared_ptr<BoxLevel> d_box_level;

   /*
    * Whether we have a globalized version of d_box_level.
    */
   mutable bool d_has_globalized_data;
   /*
    * Boxes for all level patches.
    *
    * d_boxes is slave to d_box_level and computed only if getBoxes() is called.
    * This means that the first getBoxes() has to be called by all processors,
    * because it requires communication.
    */
   mutable BoxContainer d_boxes;

   /*
    * Patch mapping to processors.
    */
   mutable ProcessorMapping d_mapping;

   /*
    * ratio to reference level
    */
   IntVector d_ratio_to_level_zero;

   /*
    * Grid geometry description.
    */
   std::shared_ptr<BaseGridGeometry> d_geometry;
   /*
    * PatchDescriptor - patch data info shared by all patches in the hierarchy
    */
   std::shared_ptr<PatchDescriptor> d_descriptor;
   /*
    * Factory for creating patches.
    */
   std::shared_ptr<PatchFactory> d_factory;

   /*
    * Local number of patches on the level.
    */
   int d_local_number_patches;

   /*
    * Extent of the index space.
    */
   std::vector<BoxContainer> d_physical_domain;

   /*
    * The ratio to coarser level applies only when the level resides
    * in a hierarchy.  The level number is that of the hierarchy level
    * that aligns with the index space of the level; if level aligns with
    * no such level then the value is -1 (default value).  The next coarser
    * level number is the next coarser level in the hierarchy for the
    * purposes of filling data from coarser levels.   It is -1 by default
    * but is usually a valid level number more often than level number.
    * The boolean is true when the level is in a hierarchy, false otherwise.
    */
   IntVector d_ratio_to_coarser_level;

   /*
    * Level number in the hierarchy.
    */
   int d_level_number;

   /*
    * Aligning with the index space of the next coarser level number.
    */
   int d_next_coarser_level_number;

   /*
    * Flag indicating the level is in a hierarchy.
    */
   bool d_in_hierarchy;

   /*
    * Container for patches.
    */
   PatchContainer d_patches;

   /*!
    * @brief Vector holding the same patches in d_patches, in the same order.
    *
    * This allows random access to the patches.
    */
   PatchVector d_patch_vector;

   /*
    * Flag to indicate boundary boxes are created.
    */
   bool d_boundary_boxes_created;

   /*!
    * @brief Has shutdown handler been initialized.
    *
    * This should be checked and set in every ctor.
    */
   static bool s_initialized;

   /*!
    * @brief Initialize static state
    */
   static bool
   initialize();

   static std::shared_ptr<tbox::Timer> t_level_constructor;
   static std::shared_ptr<tbox::Timer> t_constructor_setup;
   static std::shared_ptr<tbox::Timer> t_constructor_phys_domain;
   static std::shared_ptr<tbox::Timer> t_constructor_touch_boundaries;
   static std::shared_ptr<tbox::Timer> t_constructor_set_geometry;
   static std::shared_ptr<tbox::Timer> t_set_patch_touches;
   static std::shared_ptr<tbox::Timer> t_constructor_compute_shifts;

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;
};

}
}

#endif

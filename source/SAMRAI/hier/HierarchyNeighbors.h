/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Basic neighbor information within a hierarchy
 *
 ************************************************************************/

#ifndef included_hier_HierarchyNeighbors
#define included_hier_HierarchyNeighbors

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/PatchHierarchy.h"

#include <vector>
#include <memory>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Class providing basic information about the neighors of a Patch
 * within a PatchHierarchy.
 *
 * When a Patch exists within a PatchHierarchy, it is often useful to know
 * what other Patches in the hierarchy are its neighbors--meaning
 * those Patches on the same or an adjacent level that either overlap the
 * first Patch in space or touch it along a Patch boundary.  This class
 * provides a basic interface for identifying the neighbors of a given
 * local Patch.
 *
 * @see Connector
 */
class HierarchyNeighbors
{
public:

   /*!
    * @brief Constructor of HierarchyNeighbors
    *
    * HierarchyNeighbors is constructed based on the given PatchHierarchy
    * and a range of levels.  The range of levels may include the entire
    * hierarchy or a subset of its levels.  Any levels outside of the given
    * range will be treated by this class as if they do not exist.
    *
    * The coarsest and finest level numbers may be equal, but when this is so,
    * this class should only be used to find neighbors on the same level.
    *
    * The width argument can be used to specify a cell width such that the
    * neighbors of a patch are defined as all patches closer in proximity to
    * the patch than the given width.  If the default value of 1 is used, then
    * neighbors will be only those patches which overlap or touch at a patch
    * boundary.
    *
    * If no intra-level neighbors information is desired, set the
    * do_same_level_nbrs parameter to false.
    *
    * @param hierarchy        PatchHierarchy that will be used
    * @param coarsest_level   The coarsest level that will be used
    * @param finest_level     The finest level that will be used
    * @param do_same_level_numbers  Find intra-level neigbhors as well as
    *                               coarse-fine neighbors
    * @param width            Cell width to find neighbors
    *
    * @pre coarsest_level >= 0;
    * @pre coarsest_level <= finest_level
    * @pre finest_level < hierarchy->getNumberOfLevels()
    * @pre width >= 1
    */
   HierarchyNeighbors(const PatchHierarchy& hierarchy,
                      int coarsest_level,
                      int finest_level,
                      bool do_same_level_nbrs = true,
                      int width = 1);

   /*!
    * @brief Destructor
    */
   ~HierarchyNeighbors() = default;

   /*!
    * @brief Get the coarsest level number represented by this object
    */
   int getCoarsestLevelNumber() const
   {
      return d_coarsest_level;
   }

   /*!
    * @brief Get the finest level number represented by this object
    */
   int getFinestLevelNumber() const
   {
      return d_finest_level;
   }

   /*!
    * @brief Get neighbors from the next finer level.
    *
    * Given a Box from a local Patch on a given level, this
    * returns a container holding the neighboring boxes on the next finer
    * level.  The container will be empty if the Patch has no neighbors on
    * the finer level.
    *
    * An error will occur if the Box does not represent a local Patch on
    * the given level.
    *
    * @param box      Box representing a Patch
    * @param ln       Level number of the level holding the Patch. Both ln
    *                 and ln+1 must be within the range of levels used to
    *                 construct this object. 
    *
    * @pre ln >= d_coarsest_level && ln < d_finest_level
    */
   const BoxContainer& getFinerLevelNeighbors(
      const Box& box,
      int ln) const
   {
      TBOX_ASSERT(ln >= d_coarsest_level && ln < d_finest_level);

      std::map<BoxId, BoxContainer >::const_iterator itr =
         d_finer_level_nbrs[ln].find(box.getBoxId());

      if (itr == d_finer_level_nbrs[ln].end()) {
         TBOX_ERROR("HierarchyNeighbors::getFinerLevelNeighbors error: Box "
            << box << " does not exist locally on level " << ln << ".\n"
            << "You must specify the Box from a current local patch.");
      }

      return itr->second;
   }

   /*!
    * @brief Get neighbors from the next coarser level.
    *
    * Given a Box from a local Patch on a given level, this
    * returns a container holding the neighboring boxes on the next coarser
    * level.
    * 
    * An error will occur if the Box does not represent a local Patch on
    * the given level.
    *
    * @param box      Box representing a Patch
    * @param ln       Level number of the level holding the Patch. Both ln
    *                 and ln-1 must be within the range of levels used to
    *                 construct this object. 
    *
    * @pre ln <= d_finest_level && ln >= d_coarsest_level
    */
   const BoxContainer& getCoarserLevelNeighbors(
      const Box& box,
      int ln) const
   {
      TBOX_ASSERT(ln <= d_finest_level && ln > d_coarsest_level);

      std::map<BoxId, BoxContainer >::const_iterator itr =
         d_coarser_level_nbrs[ln].find(box.getBoxId());

      if (itr == d_coarser_level_nbrs[ln].end()) {
         TBOX_ERROR("HierarchyNeighbors::getCoarserLevelNeighbors error: Box "
            << box << " does not exist locally on level " << ln << ".\n"
            << "You must specify the Box from a current local patch.");
      }

      return itr->second;
   }

   /*!
    * @brief Get neighbors from the same level.
    *
    * Given a Box which represents a local Patch on a given level, this
    * returns a container holding the neighboring boxes on that level.
    * The container will be empty if the Patch has no neighbors.
    *
    * An error will occur if the Box does not represent a local Patch on
    * the given level.
    *
    * @param box      Box representing a Patch
    * @param ln       Level number of the level holding the Patch
    *
    * @pre ln <= d_finest_level && ln >= d_coarsest_level
    */
   const BoxContainer& getSameLevelNeighbors(const Box& box, int ln) const
   {
      TBOX_ASSERT(ln <= d_finest_level && ln >= d_coarsest_level);

      if (!d_do_same_level_nbrs) {
         TBOX_ERROR("HierarchyNeighbors::getSameLevelNeighbors error:  "
            << "HierarchyNeighbors object was constructed with argument "
            << "to not find same level neighbors." << std::endl);
      }

      std::map<BoxId, BoxContainer >::const_iterator itr =
         d_same_level_nbrs[ln].find(box.getBoxId());


      if (itr == d_same_level_nbrs[ln].end()) {
         TBOX_ERROR("HierarchyNeighbors::getSameLevelNeighbors error: Box "
            << box << " does not exist locally on level " << ln << ".\n"
            << "You must specify the BoxId of a current local patch."
            << std::endl);
      }

      return itr->second;
   }

private:

   /*!
    * Level numbers for the range of levels represented in this object.
    */
   int d_coarsest_level;
   int d_finest_level;

   /*!
    * Flag to tell whether this object has computed intra-level neighbors.
    */
   bool d_do_same_level_nbrs;

   /*!
    * @brief Containers for the neighbor boxes
    *
    * The vectors are indexed by level number, and for each level, the BoxId
    * of a Patch is mapped to a BoxContainer holding the relevant neighbors of
    * that Patch.
    */
   std::vector< std::map<BoxId, BoxContainer> > d_same_level_nbrs;
   std::vector< std::map<BoxId, BoxContainer> > d_coarser_level_nbrs;
   std::vector< std::map<BoxId, BoxContainer> > d_finer_level_nbrs;

   /*!
    * @brief Cell width that defines proximity to search for neighbors.
    */
   IntVector d_neighbor_width;
};

}
}

#endif

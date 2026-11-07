/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   A flattened representation of a hierarchy
 *
 ************************************************************************/

#ifndef included_hier_FlattenedHierarchy
#define included_hier_FlattenedHierarchy

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/HierarchyNeighbors.h"
#include "SAMRAI/hier/PatchHierarchy.h"

#include <string>
#include <vector>
#include <memory>

namespace SAMRAI {
namespace hier {

/*!
 * @brief A flattened representation of a hierarchy.
 *
 * Class FlattenedHierarchy provides a way to view a PatchHierarchy or a
 * subset of levels of PatchHierarchy in a "flattened" representation.  In
 * this representation, only the finest existing mesh is visible at any given
 * location in the problem space.
 *
 * For any Patch in a PatchHierarchy, an associated container of visible boxes
 * is created.  If the Patch has no Patches from a finer level covering any
 * of its index space, then the container will simply hold a box that is
 * spatially equivalent to the Patch itself.  If the Patch's index space is
 * entirely covered by Patches on a finer level, the container will be empty.
 * If the Patch's index space is partially covered by a finer level, then
 * the container will hold boxes that compose the uncovered parts of the Patch.
 */
class FlattenedHierarchy
{
public:

   /*!
    * @brief Constructor of FlattenedHierarchy
    *
    * A FlattenedHierarchy is constructed based on the given PatchHierarchy
    * and a range of levels.  The range of levels may include the entire
    * hierarchy or a subset of its levels.  Any levels outside of the given
    * range will be treated by this class as if they do not exist.
    *
    * If the coarsest and finest level numbers are equal, the flattened
    * representation will be equivalent to the boxes of the level.
    *
    * @param hierarchy        PatchHierarchy that will be represented
    * @param coarsest_level   The coarsest level that will be used
    * @param finest_level     The finest level that will be used
    *
    * @pre coarsest_level >= 0;
    * @pre coarsest_level <= finest_level
    * @pre finest_level < hierarchy->getNumberOfLevels()
    */
   FlattenedHierarchy(const PatchHierarchy& hierarchy,
                      int coarsest_level,
                      int finest_level);

   /*!
    * @brief Destructor
    */
   ~FlattenedHierarchy() = default;

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
    * @brief Get the visible boxes for a given Box.
    *
    * The visible boxes represent the parts of the Box that cover space that
    * is not covered by finer levels in the hierarchy.  If the entire
    * index space of the Box is covered by a finer level, the returned
    * BoxContainer will be empty.
    *
    * Each Box in the returned BoxContainer will have a BoxId that is unique
    * across all visible boxes at all levels of resolution that can be returned
    * by a FlattenedHierarchy object.  This differs from the usage of BoxId
    * it PatchLevel or BoxLevel, where a BoxId is guaranteed to be unique only
    * on one level, not across all levels.
    *
    * @param box  The Box for a local patch on the PatchLevel having
    *             level number ln.
    * @param ln   Level number of the box's level of resolution
    *
    * @pre ln <= d_finest_level && ln >= d_coarsest_level 
    */
   const BoxContainer& getVisibleBoxes(const Box& box, int ln) const
   {
      TBOX_ASSERT(ln <= d_finest_level && ln >= d_coarsest_level);

      std::map<BoxId, BoxContainer >::const_iterator itr =
         d_visible_boxes[ln].find(box.getBoxId());

      if (itr == d_visible_boxes[ln].end()) {
         TBOX_ERROR("FlattenedHierarchy::getVisibleBoxes error: Box "
            << box << " does not exist locally on level " << ln << ".\n"
            << "You must specify the Box of a current local patch.");
      }

      return itr->second;
   }

   /*
    * @brief Get the PatchHierachy associated with this FlattenedHierarchy.
    */
   const PatchHierarchy& getPatchHierarchy() const
   {
      return *d_patch_hierarchy;
   }

   /*!
    * @brief Returns an iterator to the first uncovered Box in this hierarchy.
    */
   UncoveredBoxIterator
   beginUncovered()
   {
      return UncoveredBoxIterator(this, true);
   }

   /*!
    * @brief Returns an iterator to the last uncovered Box in this hierarchy.
    */
   UncoveredBoxIterator
   endUncovered()
   {
      return UncoveredBoxIterator(this, false);
   }

private:

   /*!
    * Level numbers for the range of levels represented in this object.
    */
   int d_coarsest_level;
   int d_finest_level;

   /*!
    * @brief Container for the boxes in the flattened hierarchy representation
    *
    * The vector is indexed by level number, and for each level, the BoxId
    * of a Patch is mapped to a BoxContainer holding the visible parts of
    * that Patch's box.
    */
   std::vector< std::map<BoxId, BoxContainer> > d_visible_boxes;

   /*!
    * @brief Pointer to the PatchHierarchy that was used to create this object.
    */
   const PatchHierarchy* d_patch_hierarchy;

};

}
}

#endif

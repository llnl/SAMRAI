/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Base class for geometry management on patches
 *
 ************************************************************************/
#include "SAMRAI/hier/PatchGeometry.h"

#include "SAMRAI/hier/BoundaryLookupTable.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

PatchGeometry::PatchGeometry(
   const IntVector& ratio_to_level_zero,
   const TwoDimBool& touches_regular_bdry,
   const BlockId& block_id):
   d_dim(ratio_to_level_zero.getDim()),
   d_ratio_to_level_zero(ratio_to_level_zero),
   d_patch_boundaries(ratio_to_level_zero.getDim()),
   d_touches_regular_bdry(ratio_to_level_zero.getDim()),
   d_block_id(block_id)

{
   TBOX_ASSERT_OBJDIM_EQUALITY2(ratio_to_level_zero, touches_regular_bdry);

#ifdef DEBUG_CHECK_ASSERTIONS

   /*
    * All components of ratio must be nonzero.  Additionally, all components
    * of ratio not equal to 1 must have the same sign.
    */
   TBOX_ASSERT(ratio_to_level_zero != 0);
   if (d_dim.getValue() > 1) {
      BlockId::block_t b = block_id.getBlockValue();
      for (unsigned int i = 0; i < d_dim.getValue(); ++i) {
         bool pos0 = d_ratio_to_level_zero(b,i) > 0;
         bool pos1 = d_ratio_to_level_zero(b,(i + 1) % d_dim.getValue()) > 0;
         TBOX_ASSERT(pos0 == pos1
            || (d_ratio_to_level_zero(b,i) == 1)
            || (d_ratio_to_level_zero(b,(i + 1) % d_dim.getValue()) == 1));
      }
   }
#endif

   d_has_regular_boundary = false;
   d_has_periodic_boundary = false;

   for (int axis = 0; axis < d_dim.getValue(); ++axis) {
      for (int dir = 0; dir < 2; ++dir) {
         d_touches_regular_bdry(axis, dir) = touches_regular_bdry(axis, dir);

         if (d_touches_regular_bdry(axis, dir)) {
            d_has_regular_boundary = true;
         }
      }
   }
}

Box
PatchGeometry::getBoundaryFillBox(
   const BoundaryBox& bbox,
   const Box& patch_box,
   const IntVector& gcw) const
{

   TBOX_ASSERT_OBJDIM_EQUALITY3(bbox, patch_box, gcw);

#ifdef DEBUG_CHECK_ASSERTIONS
   for (int i = 0; i < d_dim.getValue(); ++i) {
      TBOX_ASSERT(gcw(i) >= 0);
   }
#endif
   Box tmp_box(patch_box);
   tmp_box.grow(gcw);
   Box fill_box(bbox.getBox() * tmp_box);

   int bdry_type = bbox.getBoundaryType();
   int location_index = bbox.getLocationIndex();

   // Get the singleton class lookup table
   const BoundaryLookupTable* blut;
   blut = BoundaryLookupTable::getLookupTable(d_dim);

#ifdef DEBUG_CHECK_ASSERTIONS
   const std::vector<int>& location_index_max = blut->getMaxLocationIndices();
#endif
   TBOX_ASSERT(bdry_type > 0);
   TBOX_ASSERT(bdry_type <= d_dim.getValue());
   TBOX_ASSERT(location_index >= 0);

   if (!fill_box.empty()) {

      // Loop over codimension (a.k.a. boundary type)
      for (tbox::Dimension::dir_t codim = 1; codim <= d_dim.getValue(); ++codim) {

         // When we get a match on the boundary type
         if (bdry_type == codim) {

            TBOX_ASSERT(location_index < location_index_max[codim - 1]);

            // Get the directions involved in this boundary type from the
            // lookup table.
            const std::vector<tbox::Dimension::dir_t>& dir =
               blut->getDirections(location_index, codim);

            // For each direction, identify this as an upper or lower boundary.
            for (tbox::Dimension::dir_t i = 0; i < codim; ++i) {
               if (blut->isUpper(location_index, codim, i)) {
                  fill_box.growUpper(dir[i], gcw(dir[i]) - 1);
               } else {
                  fill_box.growLower(dir[i], gcw(dir[i]) - 1);
               }
            }

            // We've found boundary type, so break out of the loop.
            break;
         }
      }
   }

   return fill_box;
}

void
PatchGeometry::setCodimensionBoundaries(
   const std::vector<BoundaryBox>& bdry_boxes,
   int codim)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   for (int i = 0; i < static_cast<int>(bdry_boxes.size()); ++i) {
      TBOX_ASSERT(bdry_boxes[i].getBoundaryType() == codim);
   }
#endif
   TBOX_ASSERT(codim <= d_dim.getValue());
   TBOX_ASSERT(codim > 0);

   d_patch_boundaries[codim - 1].clear();
   d_patch_boundaries[codim - 1].reserve(bdry_boxes.size());

   for (int b = 0; b < static_cast<int>(bdry_boxes.size()); ++b) {
      d_patch_boundaries[codim - 1].push_back(bdry_boxes[b]);
   }
}

void
PatchGeometry::setBoundaryBoxesOnPatch(
   const std::vector<std::vector<BoundaryBox> >& bdry)
{
   for (int i = 0; i < d_dim.getValue(); ++i) {
      setCodimensionBoundaries(bdry[i], i + 1);
   }
}

void
PatchGeometry::printClassData(
   std::ostream& stream) const
{
   stream << "\nPatchGeometry::printClassData..." << std::endl;
   stream << "Ratio to level zero = " << d_ratio_to_level_zero << std::endl;
   stream << "d_has_regular_boundary = "
          << d_has_regular_boundary << std::endl;
   stream << "Boundary boxes for patch..." << std::endl;
   for (int d = 0; d < d_dim.getValue(); ++d) {
      const int n = static_cast<int>(d_patch_boundaries[d].size());
      stream << "Boundary box array " << d << " has " << n << " boxes"
             << std::endl;
      for (int i = 0; i < n; ++i) {
         stream << "box " << i << " = "
                << d_patch_boundaries[d][i].getBox() << std::endl;
      }
   }
}

PatchGeometry::TwoDimBool::TwoDimBool(
   const tbox::Dimension& dim):
   d_dim(dim)
{
   setAll(false);
}

PatchGeometry::TwoDimBool::TwoDimBool(
   const tbox::Dimension& dim,
   bool v):
   d_dim(dim)
{
   for (int i = 0; i < 2 * d_dim.getValue(); ++i) {
      d_data[i] = v;
   }
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

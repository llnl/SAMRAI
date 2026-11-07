/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface to user routines for refining AMR data.
 *
 ************************************************************************/
#include "SAMRAI/mesh/MultiblockGriddingTagger.h"

#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Utilities.h"

#include <vector>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

/*
 *************************************************************************
 *
 * The default constructor and virtual destructor do nothing
 * particularly interesting.
 *
 *************************************************************************
 */

MultiblockGriddingTagger::MultiblockGriddingTagger():
   xfer::RefinePatchStrategy(),
   xfer::SingularityPatchStrategy()
{
}

MultiblockGriddingTagger::~MultiblockGriddingTagger()
{
}

hier::IntVector
MultiblockGriddingTagger::getRefineOpStencilWidth(const tbox::Dimension& dim) const
{
   return hier::IntVector::getOne(dim);
}

void
MultiblockGriddingTagger::setScratchTagPatchDataIndex(
   int buf_tag_indx)
{

   std::shared_ptr<hier::Variable> check_var;
   bool indx_maps_to_variable =
      hier::VariableDatabase::getDatabase()->mapIndexToVariable(buf_tag_indx,
         check_var);
   if (!indx_maps_to_variable || !check_var) {
      TBOX_ERROR(
         "MultiblockGriddingTagger::setScratchTagPatchDataIndex error...\n"
         << "Given patch data index = " << buf_tag_indx
         << " is not in VariableDatabase."
         << std::endl);
   } else {
      std::shared_ptr<pdat::CellVariable<int> > t_check_var(
         SAMRAI_SHARED_PTR_CAST<pdat::CellVariable<int>, hier::Variable>(check_var));
      TBOX_ASSERT(t_check_var);
   }

   d_buf_tag_indx = buf_tag_indx;
}

void
MultiblockGriddingTagger::setPhysicalBoundaryConditions(
   hier::Patch& patch,
   const double fill_time,
   const hier::IntVector& ghost_width_to_fill)
{
   NULL_USE(fill_time);

   const tbox::Dimension& dim = patch.getDim();

   const std::shared_ptr<pdat::CellData<int> > tag_data(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
         patch.getPatchData(d_buf_tag_indx)));

   TBOX_ASSERT(tag_data);

   hier::IntVector gcw =
      hier::IntVector::min(ghost_width_to_fill,
         tag_data->getGhostCellWidth());

   std::shared_ptr<hier::PatchGeometry> pgeom(patch.getPatchGeometry());

   for (int d = 0; d < dim.getValue(); ++d) {

      const std::vector<hier::BoundaryBox>& bbox =
         pgeom->getCodimensionBoundaries(d + 1);

      for (int b = 0; b < static_cast<int>(bbox.size()); ++b) {
         if (!bbox[b].getIsMultiblockSingularity()) {
            hier::Box fill_box = pgeom->getBoundaryFillBox(bbox[b],
                  patch.getBox(),
                  gcw);

            tag_data->fillAll(0, fill_box);
         }
      }
   }
}

void
MultiblockGriddingTagger::fillSingularityBoundaryConditions(
   hier::Patch& patch,
   const hier::PatchLevel& encon_level,
   std::shared_ptr<const hier::Connector> dst_to_encon,
   const hier::Box& fill_box,
   const hier::BoundaryBox& boundary_box,
   const std::shared_ptr<hier::BaseGridGeometry>& grid_geometry)
{
   NULL_USE(boundary_box);
   NULL_USE(grid_geometry);

   TBOX_ASSERT(!grid_geometry->hasEnhancedConnectivity() || dst_to_encon);
   TBOX_ASSERT_OBJDIM_EQUALITY3(patch, fill_box, boundary_box);

   const tbox::Dimension& dim = fill_box.getDim();

   const hier::BoxId& dst_mb_id = patch.getBox().getBoxId();

   const hier::BlockId& patch_blk_id = patch.getBox().getBlockId();

   const std::shared_ptr<pdat::CellData<int> > tag_data(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
         patch.getPatchData(d_buf_tag_indx)));

   TBOX_ASSERT(tag_data);

   hier::Box sing_fill_box(tag_data->getGhostBox() * fill_box);
   tag_data->fillAll(0, sing_fill_box);

   if (grid_geometry->hasEnhancedConnectivity()) {

      hier::Connector::ConstNeighborhoodIterator ni =
         dst_to_encon->findLocal(dst_mb_id);

      if (ni != dst_to_encon->end()) {

         for (hier::Connector::ConstNeighborIterator ei = dst_to_encon->begin(ni);
              ei != dst_to_encon->end(ni); ++ei) {

            std::shared_ptr<hier::Patch> encon_patch(
               encon_level.getPatch(ei->getBoxId()));

            const hier::BlockId& encon_blk_id = ei->getBlockId();

            hier::Transformation::RotationIdentifier rotation =
               hier::Transformation::NO_ROTATE;
            hier::IntVector offset(dim);

            hier::BaseGridGeometry::ConstNeighborIterator itr =
               grid_geometry->find(patch_blk_id, encon_blk_id);
            if (itr != grid_geometry->end(patch_blk_id)) {
               rotation = (*itr).getRotationIdentifier();
               offset = (*itr).getShift(encon_level.getLevelNumber());
            }

            hier::Transformation transformation(
               rotation, offset, encon_blk_id, patch_blk_id);
            hier::Box encon_patch_box(encon_patch->getBox());
            transformation.transform(encon_patch_box);

            hier::Box encon_fill_box(encon_patch_box * sing_fill_box);
            if (!encon_fill_box.empty()) {

               const hier::Transformation::RotationIdentifier back_rotate =
                  hier::Transformation::getReverseRotationIdentifier(
                     rotation, dim);

               hier::IntVector back_shift(dim);

               hier::Transformation::calculateReverseShift(
                  back_shift, offset, rotation);

               hier::Transformation back_trans(back_rotate, back_shift,
                                               encon_fill_box.getBlockId(),
                                               encon_patch->getBox().getBlockId());

               std::shared_ptr<pdat::CellData<int> > sing_data(
                  SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
                     encon_patch->getPatchData(d_buf_tag_indx)));

               TBOX_ASSERT(sing_data);

               pdat::CellIterator ciend(pdat::CellGeometry::end(encon_fill_box));
               for (pdat::CellIterator ci(pdat::CellGeometry::begin(encon_fill_box));
                    ci != ciend; ++ci) {
                  pdat::CellIndex src_index(*ci);
                  pdat::CellGeometry::transform(src_index, back_trans);

                  int sing_val = (*sing_data)(src_index);
                  if (sing_val != 0 && (*tag_data)(*ci) == 0) {
                     (*tag_data)(*ci) = sing_val;
                  }
               }
            }
         }
      }
   }
}

void
MultiblockGriddingTagger::preprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::Box& fine_box,
   const hier::IntVector& ratio)
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_box);
   NULL_USE(ratio);
}

void
MultiblockGriddingTagger::postprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::Box& fine_box,
   const hier::IntVector& ratio)
{
   TBOX_ASSERT_OBJDIM_EQUALITY4(fine, coarse, fine_box, ratio);

   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_box);
   NULL_USE(ratio);
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

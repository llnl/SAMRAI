/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   AMR communication tests for cell-centered patch data
 *
 ************************************************************************/

#include "CellMultiblockTest.h"

#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellVariable.h"
//#include "SAMRAI/pdat/CellDoubleConstantRefine.h"
#include "SAMRAI/pdat/CellConstantRefine.h"
#include "SAMRAI/tbox/NVTXUtilities.h"

#include "MultiblockTester.h"

#include <vector>

using namespace SAMRAI;

using CELL_MBLK_KERNEL_TYPE = double;

CellMultiblockTest::CellMultiblockTest(
   const std::string& object_name,
   const tbox::Dimension& dim,
   std::shared_ptr<tbox::Database> main_input_db,
   const std::string& refine_option):
   PatchMultiblockTestStrategy(dim),
   d_dim(dim)
{

   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(main_input_db);
   TBOX_ASSERT(!refine_option.empty());

   d_object_name = object_name;

   d_refine_option = refine_option;

   d_finest_level_number = main_input_db->
      getDatabase("PatchHierarchy")->
      getInteger("max_levels") - 1;

   //int num_blocks = main_input_db->getDatabase("BlockGridGeometry")->
   //   getInteger("num_blocks");

   //d_skel_grid_geometry.resizeArray(num_blocks);

   std::string geom_name("BlockGridGeometry");

   if (main_input_db->keyExists(geom_name)) {
      getGridGeometry().reset(
         new geom::GridGeometry(
            dim,
            geom_name,
            main_input_db->getDatabase(geom_name)));

   } else {
      TBOX_ERROR("CellMultiblockTest: could not find entry `"
         << geom_name << "' in input.");
   }

   readTestInput(main_input_db->getDatabase("CellMultiblockTest"));
}

CellMultiblockTest::~CellMultiblockTest()
{

}

void CellMultiblockTest::readTestInput(
   std::shared_ptr<tbox::Database> db)
{
   TBOX_ASSERT(db);

   /*
    * Base class reads variable parameters and boxes to refine.
    */

   readVariableInput(db->getDatabase("VariableData"));
   readRefinementInput(db->getDatabase("RefinementData"));
}

void CellMultiblockTest::registerVariables(
   MultiblockTester* commtest)
{
   TBOX_ASSERT(commtest != 0);

   int nvars = static_cast<int>(d_variable_src_name.size());

   d_variables.resize(nvars);

   for (int i = 0; i < nvars; ++i) {
      d_variables[i].reset(
         new pdat::CellVariable<CELL_MBLK_KERNEL_TYPE>(d_dim,
            d_variable_src_name[i],
            d_variable_depth[i]));

      commtest->registerVariable(d_variables[i],
         d_variables[i],
         d_variable_src_ghosts[i],
         d_variable_dst_ghosts[i],
         getGridGeometry(),
         d_variable_refine_op[i]);

   }

}

void CellMultiblockTest::initializeDataOnPatch(
   hier::Patch& patch,
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   const int level_number,
   const hier::BlockId& block_id,
   char src_or_dst)
{
   NULL_USE(hierarchy);
   NULL_USE(src_or_dst);

   if ((d_refine_option == "INTERIOR_FROM_SAME_LEVEL")
       || ((d_refine_option == "INTERIOR_FROM_COARSER_LEVEL")
           && (level_number < d_finest_level_number))) {

      for (int i = 0; i < static_cast<int>(d_variables.size()); ++i) {

         std::shared_ptr<pdat::CellData<CELL_MBLK_KERNEL_TYPE> > cell_data(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<CELL_MBLK_KERNEL_TYPE>, hier::PatchData>(
               patch.getPatchData(d_variables[i], getDataContext())));
         TBOX_ASSERT(cell_data);

         hier::Box dbox = cell_data->getGhostBox();
         cell_data->fillAll((CELL_MBLK_KERNEL_TYPE)block_id.getBlockValue());
      }
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
   }
}

void CellMultiblockTest::tagCellsToRefine(
   hier::Patch& patch,
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number,
   int tag_index)
{
   NULL_USE(hierarchy);

   /*
    * Base class sets tags in box array for each level.
    */
   tagCellsInInputBoxes(patch, level_number, tag_index);

}

void CellMultiblockTest::setPhysicalBoundaryConditions(
   hier::Patch& patch,
   const double time,
   const hier::IntVector& gcw_to_fill) const
{
   NULL_USE(time);

   std::shared_ptr<hier::PatchGeometry> pgeom(patch.getPatchGeometry());

   const std::vector<hier::BoundaryBox>& node_bdry =
      pgeom->getCodimensionBoundaries(d_dim.getValue());
   const int num_node_bdry_boxes = static_cast<int>(node_bdry.size());

   std::vector<hier::BoundaryBox> empty_vector(0, hier::BoundaryBox(d_dim));
   const std::vector<hier::BoundaryBox>& edge_bdry =
      d_dim > tbox::Dimension(1) ?
      pgeom->getCodimensionBoundaries(d_dim.getValue() - 1) : empty_vector;
   const int num_edge_bdry_boxes = static_cast<int>(edge_bdry.size());

   const std::vector<hier::BoundaryBox>& face_bdry =
      d_dim == tbox::Dimension(3) ?
      pgeom->getCodimensionBoundaries(d_dim.getValue() - 2) : empty_vector;
   const int num_face_bdry_boxes = static_cast<int>(face_bdry.size());

   for (int i = 0; i < static_cast<int>(d_variables.size()); ++i) {

      std::shared_ptr<pdat::CellData<CELL_MBLK_KERNEL_TYPE> > cell_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<CELL_MBLK_KERNEL_TYPE>, hier::PatchData>(
            patch.getPatchData(d_variables[i], getDataContext())));
      TBOX_ASSERT(cell_data);

      /*
       * Set node boundary data.
       */
      for (int ni = 0; ni < num_node_bdry_boxes; ++ni) {

         hier::Box fill_box = pgeom->getBoundaryFillBox(node_bdry[ni],
               patch.getBox(),
               gcw_to_fill);

         if (!node_bdry[ni].getIsMultiblockSingularity()) {
            cell_data->fillAll((CELL_MBLK_KERNEL_TYPE)(node_bdry[ni].getLocationIndex() + 100),
               fill_box);
         }
      }

      if (d_dim > tbox::Dimension(1)) {
         /*
          * Set edge boundary data.
          */
         for (int ei = 0; ei < num_edge_bdry_boxes; ++ei) {

            hier::Box fill_box = pgeom->getBoundaryFillBox(edge_bdry[ei],
                  patch.getBox(),
                  gcw_to_fill);

            if (!edge_bdry[ei].getIsMultiblockSingularity()) {
               cell_data->fillAll((CELL_MBLK_KERNEL_TYPE)(edge_bdry[ei].getLocationIndex()
                                           + 100),
                  fill_box);
            }
         }
      }

      if (d_dim == tbox::Dimension(3)) {
         /*
          * Set face boundary data.
          */
         for (int fi = 0; fi < num_face_bdry_boxes; ++fi) {

            hier::Box fill_box = pgeom->getBoundaryFillBox(face_bdry[fi],
                  patch.getBox(),
                  gcw_to_fill);

            if (!face_bdry[fi].getIsMultiblockSingularity()) {
               cell_data->fillAll((CELL_MBLK_KERNEL_TYPE)(face_bdry[fi].getLocationIndex()
                                           + 100),
                  fill_box);
            }
         }
      }

   }
#if defined(HAVE_RAJA)
   tbox::parallel_synchronize();
#endif

}

void CellMultiblockTest::fillSingularityBoundaryConditions(
   hier::Patch& patch,
   const hier::PatchLevel& encon_level,
   std::shared_ptr<const hier::Connector> dst_to_encon,
   const hier::Box& fill_box,
   const hier::BoundaryBox& bbox,
   const std::shared_ptr<hier::BaseGridGeometry>& grid_geometry)
{
   const tbox::Dimension& dim = fill_box.getDim();

   const hier::BoxId& dst_mb_id = patch.getBox().getBoxId();

   const hier::BlockId& patch_blk_id = patch.getBox().getBlockId();

   for (int i = 0; i < static_cast<int>(d_variables.size()); ++i) {

      std::shared_ptr<pdat::CellData<CELL_MBLK_KERNEL_TYPE> > cell_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<CELL_MBLK_KERNEL_TYPE>, hier::PatchData>(
            patch.getPatchData(d_variables[i], getDataContext())));
      TBOX_ASSERT(cell_data);

      hier::Box sing_fill_box(cell_data->getGhostBox() * fill_box);
      cell_data->fillAll(0.0, sing_fill_box);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif


      int depth = cell_data->getDepth();
      int num_encon_used = 0;

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

               hier::Transformation transformation(rotation, offset,
                                                   encon_blk_id,
                                                   patch_blk_id);
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
                                                  patch_blk_id,
                                                  encon_blk_id);

                  std::shared_ptr<pdat::CellData<CELL_MBLK_KERNEL_TYPE> > sing_data(
                     SAMRAI_SHARED_PTR_CAST<pdat::CellData<CELL_MBLK_KERNEL_TYPE>, hier::PatchData>(
                        encon_patch->getPatchData(
                           d_variables[i], getDataContext())));
                  TBOX_ASSERT(sing_data);

                  pdat::CellIterator ciend(pdat::CellGeometry::end(encon_fill_box));
                  for (pdat::CellIterator ci(pdat::CellGeometry::begin(encon_fill_box));
                       ci != ciend; ++ci) {
                     pdat::CellIndex src_index(*ci);
                     pdat::CellGeometry::transform(src_index, back_trans);
                     for (int d = 0; d < depth; ++d) {
                        (*cell_data)(*ci, d) += (*sing_data)(src_index, d);
                     }
                  }

                  ++num_encon_used;
               }
            }
         }
      }

      if (num_encon_used) {
         pdat::CellIterator ciend(pdat::CellGeometry::end(sing_fill_box));
         for (pdat::CellIterator ci(pdat::CellGeometry::begin(sing_fill_box));
              ci != ciend; ++ci) {
            for (int d = 0; d < depth; ++d) {
               (*cell_data)(*ci, d) /= num_encon_used;
            }
         }
      } else {
         cell_data->fillAll(
            (CELL_MBLK_KERNEL_TYPE)bbox.getLocationIndex() + 200.0, fill_box);
      }
   }
#if defined(HAVE_RAJA)
   tbox::parallel_synchronize();
#endif
}

/*
 *************************************************************************
 *
 * Verify results of communication operations.  This test must be
 * consistent with data initialization and boundary operations above.
 *
 *************************************************************************
 */
bool CellMultiblockTest::verifyResults(
   const hier::Patch& patch,
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   const int level_number,
   const hier::BlockId& block_id)
{

   tbox::plog << "\nEntering CellMultiblockTest::verifyResults..." << std::endl;
   tbox::plog << "level_number = " << level_number << std::endl;
   tbox::plog << "Patch box = " << patch.getBox() << std::endl;

   hier::IntVector tgcw(d_dim, 0);
   for (int i = 0; i < static_cast<int>(d_variables.size()); ++i) {
      tgcw.max(patch.getPatchData(d_variables[i], getDataContext())->
         getGhostCellWidth());
   }
   hier::Box pbox = patch.getBox();

   std::shared_ptr<pdat::CellData<CELL_MBLK_KERNEL_TYPE> > solution(
      new pdat::CellData<CELL_MBLK_KERNEL_TYPE>(pbox, 1, tgcw));

   hier::Box tbox(pbox);
   tbox.grow(tgcw);

   std::shared_ptr<hier::BaseGridGeometry> grid_geom(
      hierarchy->getGridGeometry());
   hier::BoxContainer singularity(
      grid_geom->getSingularityBoxContainer(block_id));

   hier::IntVector ratio =
      hierarchy->getPatchLevel(level_number)->getRatioToLevelZero();

   singularity.refine(ratio);

   bool test_failed = false;

   for (int i = 0; i < static_cast<int>(d_variables.size()); ++i) {

      CELL_MBLK_KERNEL_TYPE correct = (CELL_MBLK_KERNEL_TYPE)block_id.getBlockValue();

      std::shared_ptr<pdat::CellData<CELL_MBLK_KERNEL_TYPE> > cell_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<CELL_MBLK_KERNEL_TYPE>, hier::PatchData>(
            patch.getPatchData(d_variables[i], getDataContext())));
      TBOX_ASSERT(cell_data);
      int depth = cell_data->getDepth();

      pdat::CellIterator ciend(pdat::CellGeometry::end(pbox));
      for (pdat::CellIterator ci(pdat::CellGeometry::begin(pbox));
           ci != ciend; ++ci) {
         for (int d = 0; d < depth; ++d) {
            CELL_MBLK_KERNEL_TYPE result = (*cell_data)(*ci, d);

            if (!tbox::MathUtilities<CELL_MBLK_KERNEL_TYPE>::equalEps(correct, result)) {
               tbox::perr << "Test FAILED: ...."
                          << " : cell index = " << *ci << std::endl;
               tbox::perr << "    Variable = " << d_variable_src_name[i]
                          << " : depth index = " << d << std::endl;
               tbox::perr << "    result = " << result
                          << " : correct = " << correct << std::endl;
               test_failed = true;
            }
         }
      }

      hier::Box gbox = cell_data->getGhostBox();

      for (hier::BaseGridGeometry::ConstNeighborIterator ne =
              grid_geom->begin(block_id);
           ne != grid_geom->end(block_id); ++ne) {

         const hier::BaseGridGeometry::Neighbor& nbr = *ne;
         correct = nbr.getBlockId().getBlockValue();

         if (nbr.isSingularity()) {
            continue;
         }

         hier::BoxContainer neighbor_ghost(nbr.getTransformedDomain());
         neighbor_ghost.refine(ratio);
         neighbor_ghost.intersectBoxes(gbox);

         for (hier::BoxContainer::iterator ng = neighbor_ghost.begin();
              ng != neighbor_ghost.end(); ++ng) {

            pdat::CellIterator ciend(pdat::CellGeometry::end(*ng));
            for (pdat::CellIterator ci(pdat::CellGeometry::begin(*ng));
                 ci != ciend; ++ci) {
               for (int d = 0; d < depth; ++d) {
                  CELL_MBLK_KERNEL_TYPE result = (*cell_data)(*ci, d);

                  if (!tbox::MathUtilities<CELL_MBLK_KERNEL_TYPE>::equalEps(correct,
                         result)) {
                     tbox::perr << "Test FAILED: ...."
                                << " : cell index = " << *ci << std::endl;
                     tbox::perr << "    Variable = " << d_variable_src_name[i]
                                << " : depth index = " << d << std::endl;
                     tbox::perr << "    result = " << result
                                << " : correct = " << correct << std::endl;
                     test_failed = true;
                  }
               }
            }
         }
      }

      std::shared_ptr<hier::PatchGeometry> pgeom = patch.getPatchGeometry();
      TBOX_ASSERT(pgeom);

      for (int b = 0; b < d_dim.getValue(); ++b) {
         const std::vector<hier::BoundaryBox>& bdry =
            pgeom->getCodimensionBoundaries(b + 1);

         for (int k = 0; k < static_cast<int>(bdry.size()); ++k) {
            hier::Box fill_box = pgeom->getBoundaryFillBox(bdry[k],
                  patch.getBox(),
                  tgcw);
            fill_box = fill_box * gbox;

            if (bdry[k].getIsMultiblockSingularity()) {
               correct = 0.0;

               int num_sing_neighbors = 0;
               for (hier::BaseGridGeometry::ConstNeighborIterator ns(
                       grid_geom->begin(block_id));
                    ns != grid_geom->end(block_id); ++ns) {
                  const hier::BaseGridGeometry::Neighbor& nbr = *ns;
                  if (nbr.isSingularity()) {
                     hier::BoxContainer neighbor_ghost(
                        nbr.getTransformedDomain());
                     neighbor_ghost.refine(ratio);
                     neighbor_ghost.intersectBoxes(fill_box);
                     if (neighbor_ghost.size()) {
                        ++num_sing_neighbors;
                        correct += nbr.getBlockId().getBlockValue();
                     }
                  }
               }

               if (num_sing_neighbors == 0) {

                  correct = (CELL_MBLK_KERNEL_TYPE)bdry[k].getLocationIndex() + 200.0;

               } else {

                  correct /= (CELL_MBLK_KERNEL_TYPE)num_sing_neighbors;

               }

            } else {
               correct = (CELL_MBLK_KERNEL_TYPE)(bdry[k].getLocationIndex() + 100);
            }

            pdat::CellIterator ciend(pdat::CellGeometry::end(fill_box));
            for (pdat::CellIterator ci(pdat::CellGeometry::begin(fill_box));
                 ci != ciend; ++ci) {
               for (int d = 0; d < depth; ++d) {
                  CELL_MBLK_KERNEL_TYPE result = (*cell_data)(*ci, d);

                  if (!tbox::MathUtilities<CELL_MBLK_KERNEL_TYPE>::equalEps(correct,
                         result)) {
                     tbox::perr << "Test FAILED: ...."
                                << " : cell index = " << *ci << std::endl;
                     tbox::perr << "    Variable = " << d_variable_src_name[i]
                                << " : depth index = " << d << std::endl;
                     tbox::perr << "    result = " << result
                                << " : correct = " << correct << std::endl;
                     test_failed = true;
                  }
               }
            }

         }
      }

   }

   if (!test_failed) {
      tbox::plog << "CellMultiblockTest Successful!" << std::endl;
   } else {
      tbox::perr << "Multiblock CellMultiblockTest FAILED: \n" << std::endl;
   }

   solution.reset();   // just to be anal...

   tbox::plog << "\nExiting CellMultiblockTest::verifyResults..." << std::endl;
   tbox::plog << "level_number = " << level_number << std::endl;
   tbox::plog << "Patch box = " << patch.getBox() << std::endl << std::endl;

   return !test_failed;
}

/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   AMR communication tests for face-centered patch data
 *
 ************************************************************************/

#include "FaceMultiblockTest.h"

#include "SAMRAI/xfer/BoxGeometryVariableFillPattern.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/FaceDoubleConstantRefine.h"
#include "SAMRAI/pdat/FaceConstantRefine.h"
#include "SAMRAI/pdat/FaceVariable.h"

#include "MultiblockTester.h"

#include <vector>

using namespace SAMRAI;

using FACE_MBLK_KERNEL_TYPE = double;

FaceMultiblockTest::FaceMultiblockTest(
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

   std::string geom_name("BlockGridGeometry");

   if (main_input_db->keyExists(geom_name)) {
      getGridGeometry().reset(
         new geom::GridGeometry(
            dim,
            geom_name,
            main_input_db->getDatabase(geom_name)));

   } else {
      TBOX_ERROR("FaceMultiblockTest: could not find entry `"
         << geom_name << "' in input.");
   }

   readTestInput(main_input_db->getDatabase("FaceMultiblockTest"));
}

FaceMultiblockTest::~FaceMultiblockTest()
{
}

void FaceMultiblockTest::readTestInput(
   std::shared_ptr<tbox::Database> db)
{
   TBOX_ASSERT(db);

   /*
    * Base class reads variable parameters and boxes to refine.
    */

   readVariableInput(db->getDatabase("VariableData"));
   readRefinementInput(db->getDatabase("RefinementData"));
}

void FaceMultiblockTest::registerVariables(
   MultiblockTester* commtest)
{
   TBOX_ASSERT(commtest != 0);

   int nvars = static_cast<int>(d_variable_src_name.size());

   d_variables.resize(nvars);

   for (int i = 0; i < nvars; ++i) {
      d_variables[i].reset(
         new pdat::FaceVariable<FACE_MBLK_KERNEL_TYPE>(d_dim,
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

void FaceMultiblockTest::initializeDataOnPatch(
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

         std::shared_ptr<pdat::FaceData<FACE_MBLK_KERNEL_TYPE> > face_data(
            SAMRAI_SHARED_PTR_CAST<pdat::FaceData<FACE_MBLK_KERNEL_TYPE>, hier::PatchData>(
               patch.getPatchData(d_variables[i], getDataContext())));
         TBOX_ASSERT(face_data);

         hier::Box dbox = face_data->getGhostBox();

         face_data->fillAll((FACE_MBLK_KERNEL_TYPE)block_id.getBlockValue());

      }
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
   }
}

void FaceMultiblockTest::tagCellsToRefine(
   hier::Patch& patch,
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   const int level_number,
   int tag_index)
{
   NULL_USE(hierarchy);

   /*
    * Base class sets tags in box array for each level.
    */
   tagCellsInInputBoxes(patch, level_number, tag_index);

}

void FaceMultiblockTest::setPhysicalBoundaryConditions(
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

      std::shared_ptr<pdat::FaceData<FACE_MBLK_KERNEL_TYPE> > face_data(
         SAMRAI_SHARED_PTR_CAST<pdat::FaceData<FACE_MBLK_KERNEL_TYPE>, hier::PatchData>(
            patch.getPatchData(d_variables[i], getDataContext())));
      TBOX_ASSERT(face_data);

      /*
       * Set node boundary data.
       */
      for (int nb = 0; nb < num_node_bdry_boxes; ++nb) {

         hier::Box fill_box = pgeom->getBoundaryFillBox(node_bdry[nb],
               patch.getBox(),
               gcw_to_fill);

         for (tbox::Dimension::dir_t axis = 0; axis < d_dim.getValue(); ++axis) {
            hier::Box patch_face_box =
               pdat::FaceGeometry::toFaceBox(patch.getBox(), axis);
            if (!node_bdry[nb].getIsMultiblockSingularity()) {
               pdat::FaceIterator niend(pdat::FaceGeometry::end(fill_box, axis));
               for (pdat::FaceIterator ni(pdat::FaceGeometry::begin(fill_box, axis));
                    ni != niend; ++ni) {
                  if (!patch_face_box.contains(*ni)) {
                     for (int d = 0; d < face_data->getDepth(); ++d) {
                        (*face_data)(*ni, d) =
                           (FACE_MBLK_KERNEL_TYPE)(node_bdry[nb].getLocationIndex() + 100);
                     }
                  }
               }
            }
         }
      }

      if (d_dim > tbox::Dimension(1)) {
         /*
          * Set edge boundary data.
          */
         for (int eb = 0; eb < num_edge_bdry_boxes; ++eb) {

            hier::Box fill_box = pgeom->getBoundaryFillBox(edge_bdry[eb],
                  patch.getBox(),
                  gcw_to_fill);

            for (tbox::Dimension::dir_t axis = 0; axis < d_dim.getValue(); ++axis) {
               hier::Box patch_face_box =
                  pdat::FaceGeometry::toFaceBox(patch.getBox(), axis);
               hier::Index plower(patch_face_box.lower());
               hier::Index pupper(patch_face_box.upper());

               if (!edge_bdry[eb].getIsMultiblockSingularity()) {
                  pdat::FaceIterator niend(pdat::FaceGeometry::end(fill_box, axis));
                  for (pdat::FaceIterator ni(pdat::FaceGeometry::begin(fill_box, axis));
                       ni != niend; ++ni) {
                     if (!patch_face_box.contains(*ni)) {
                        bool use_index = true;
                        for (tbox::Dimension::dir_t n = 0; n < d_dim.getValue(); ++n) {
                           if (axis == n &&
                               edge_bdry[eb].getBox().numberCells(n) == 1) {
                              if ((*ni)(0) == plower(0) || (*ni)(0) ==
                                  pupper(0)) {
                                 use_index = false;
                                 break;
                              }
                           }
                        }

                        if (use_index) {
                           for (int d = 0; d < face_data->getDepth(); ++d) {
                              (*face_data)(*ni, d) =
                                 (FACE_MBLK_KERNEL_TYPE)(edge_bdry[eb].getLocationIndex()
                                          + 100);
                           }
                        }
                     }
                  }
               }
            }
         }
      }

      if (d_dim == tbox::Dimension(3)) {
         /*
          * Set face boundary data.
          */
         for (int fb = 0; fb < num_face_bdry_boxes; ++fb) {

            hier::Box fill_box = pgeom->getBoundaryFillBox(face_bdry[fb],
                  patch.getBox(),
                  gcw_to_fill);

            for (tbox::Dimension::dir_t axis = 0; axis < d_dim.getValue(); ++axis) {
               hier::Box patch_face_box =
                  pdat::FaceGeometry::toFaceBox(patch.getBox(), axis);
               hier::Index plower(patch_face_box.lower());
               hier::Index pupper(patch_face_box.upper());

               if (!face_bdry[fb].getIsMultiblockSingularity()) {
                  pdat::FaceIterator niend(pdat::FaceGeometry::end(fill_box, axis));
                  for (pdat::FaceIterator ni(pdat::FaceGeometry::begin(fill_box, axis));
                       ni != niend; ++ni) {
                     if (!patch_face_box.contains(*ni)) {
                        bool use_index = true;
                        for (tbox::Dimension::dir_t n = 0; n < d_dim.getValue(); ++n) {
                           if (axis == n &&
                               face_bdry[fb].getBox().numberCells(n) == 1) {
                              if ((*ni)(0) == plower(0) || (*ni)(0) ==
                                  pupper(0)) {
                                 use_index = false;
                                 break;
                              }
                           }
                        }

                        if (use_index) {
                           for (int d = 0; d < face_data->getDepth(); ++d) {
                              (*face_data)(*ni, d) =
                                 (FACE_MBLK_KERNEL_TYPE)(face_bdry[fb].getLocationIndex()
                                          + 100);
                           }
                        }
                     }
                  }
               }
            }
         }
      }

   }

}

void FaceMultiblockTest::fillSingularityBoundaryConditions(
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

      std::shared_ptr<pdat::FaceData<FACE_MBLK_KERNEL_TYPE> > face_data(
         SAMRAI_SHARED_PTR_CAST<pdat::FaceData<FACE_MBLK_KERNEL_TYPE>, hier::PatchData>(
            patch.getPatchData(d_variables[i], getDataContext())));
      TBOX_ASSERT(face_data);

      hier::Box sing_fill_box(face_data->getGhostBox() * fill_box);

      int depth = face_data->getDepth();

      for (tbox::Dimension::dir_t axis = 0; axis < d_dim.getValue(); ++axis) {
         hier::Box pbox =
            pdat::FaceGeometry::toFaceBox(patch.getBox(), axis);

         hier::Index plower(pbox.lower());
         hier::Index pupper(pbox.upper());

         pdat::FaceIterator niend(pdat::FaceGeometry::end(sing_fill_box, axis));
         for (pdat::FaceIterator ni(pdat::FaceGeometry::begin(sing_fill_box, axis));
              ni != niend; ++ni) {
            bool use_index = true;
            for (tbox::Dimension::dir_t n = 0; n < d_dim.getValue(); ++n) {
               if (axis == n && bbox.getBox().numberCells(n) == 1) {
                  if ((*ni)(0) == plower(0) || (*ni)(0) == pupper(0)) {
                     use_index = false;
                     break;
                  }
               }
            }
            if (use_index) {
               for (int d = 0; d < depth; ++d) {
                  (*face_data)(*ni, d) = 0.0;
               }
            }
         }
      }

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

                  std::shared_ptr<pdat::FaceData<FACE_MBLK_KERNEL_TYPE> > sing_data(
                     SAMRAI_SHARED_PTR_CAST<pdat::FaceData<FACE_MBLK_KERNEL_TYPE>, hier::PatchData>(
                        encon_patch->getPatchData(
                           d_variables[i], getDataContext())));
                  TBOX_ASSERT(sing_data);

                  for (tbox::Dimension::dir_t axis = 0; axis < d_dim.getValue(); ++axis) {

                     hier::Box pbox(
                        pdat::FaceGeometry::toFaceBox(patch.getBox(), axis));

                     hier::Index plower(pbox.lower());
                     hier::Index pupper(pbox.upper());

                     pdat::FaceIterator ciend(pdat::FaceGeometry::end(sing_fill_box, axis));
                     for (pdat::FaceIterator ci(pdat::FaceGeometry::begin(sing_fill_box, axis));
                          ci != ciend; ++ci) {
                        bool use_index = true;
                        for (tbox::Dimension::dir_t n = 0; n < d_dim.getValue(); ++n) {
                           if (axis == n && bbox.getBox().numberCells(n) == 1) {
                              if ((*ci)(0) == plower(0) || (*ci)(0) == pupper(0)) {
                                 use_index = false;
                                 break;
                              }
                           }
                        }
                        if (use_index) {

                           pdat::FaceIndex src_index(*ci);
                           pdat::FaceGeometry::transform(src_index, back_trans);

                           for (int d = 0; d < depth; ++d) {
                              (*face_data)(*ci, d) += (*sing_data)(src_index, d);
                           }
                        }
                     }
                  }

                  ++num_encon_used;
               }
            }
         }
      }

      if (num_encon_used) {

         for (tbox::Dimension::dir_t axis = 0; axis < d_dim.getValue(); ++axis) {

            hier::Box pbox =
               pdat::FaceGeometry::toFaceBox(patch.getBox(), axis);

            hier::Index plower(pbox.lower());
            hier::Index pupper(pbox.upper());

            pdat::FaceIterator ciend(pdat::FaceGeometry::end(sing_fill_box, axis));
            for (pdat::FaceIterator ci(pdat::FaceGeometry::begin(sing_fill_box, axis));
                 ci != ciend; ++ci) {
               bool use_index = true;
               for (tbox::Dimension::dir_t n = 0; n < d_dim.getValue(); ++n) {
                  if (axis == n && bbox.getBox().numberCells(n) == 1) {
                     if ((*ci)(0) == plower(0) || (*ci)(0) == pupper(0)) {
                        use_index = false;
                        break;
                     }
                  }
               }
               if (use_index) {
                  for (int d = 0; d < depth; ++d) {
                     (*face_data)(*ci, d) /= num_encon_used;
                  }
               }
            }
         }

      } else {

         /*
          * In cases of reduced connectivity, there are no other blocks
          * from which to acquire data.
          */

         for (tbox::Dimension::dir_t axis = 0; axis < d_dim.getValue(); ++axis) {

            hier::Box pbox =
               pdat::FaceGeometry::toFaceBox(patch.getBox(), axis);

            hier::Index plower(pbox.lower());
            hier::Index pupper(pbox.upper());

            pdat::FaceIterator ciend(pdat::FaceGeometry::end(sing_fill_box, axis));
            for (pdat::FaceIterator ci(pdat::FaceGeometry::begin(sing_fill_box, axis));
                 ci != ciend; ++ci) {
               bool use_index = true;
               for (tbox::Dimension::dir_t n = 0; n < d_dim.getValue(); ++n) {
                  if (axis == n && bbox.getBox().numberCells(n) == 1) {
                     if ((*ci)(0) == plower(0) || (*ci)(0) == pupper(0)) {
                        use_index = false;
                        break;
                     }
                  }
               }
               if (use_index) {
                  for (int d = 0; d < depth; ++d) {
                     (*face_data)(*ci, d) =
                        (FACE_MBLK_KERNEL_TYPE)bbox.getLocationIndex() + 200.0;
                  }
               }
            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Verify results of communication operations.  This test must be
 * consistent with data initialization and boundary operations above.
 *
 *************************************************************************
 */
bool FaceMultiblockTest::verifyResults(
   const hier::Patch& patch,
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number,
   const hier::BlockId& block_id)
{

   tbox::plog << "\nEntering FaceMultiblockTest::verifyResults..." << std::endl;
   tbox::plog << "level_number = " << level_number << std::endl;
   tbox::plog << "Patch box = " << patch.getBox() << std::endl;

   hier::IntVector tgcw(d_dim, 0);
   for (int i = 0; i < static_cast<int>(d_variables.size()); ++i) {
      tgcw.max(patch.getPatchData(d_variables[i], getDataContext())->
         getGhostCellWidth());
   }
   hier::Box pbox = patch.getBox();

   std::shared_ptr<pdat::FaceData<FACE_MBLK_KERNEL_TYPE> > solution(
      new pdat::FaceData<FACE_MBLK_KERNEL_TYPE>(pbox, 1, tgcw));

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

      FACE_MBLK_KERNEL_TYPE correct = (FACE_MBLK_KERNEL_TYPE)block_id.getBlockValue();

      std::shared_ptr<pdat::FaceData<FACE_MBLK_KERNEL_TYPE> > face_data(
         SAMRAI_SHARED_PTR_CAST<pdat::FaceData<FACE_MBLK_KERNEL_TYPE>, hier::PatchData>(
            patch.getPatchData(d_variables[i], getDataContext())));
      TBOX_ASSERT(face_data);
      int depth = face_data->getDepth();

      hier::Box interior_box(pbox);
      interior_box.grow(hier::IntVector(d_dim, -1));

      for (tbox::Dimension::dir_t axis = 0; axis < d_dim.getValue(); ++axis) {
         pdat::FaceIterator ciend(pdat::FaceGeometry::end(interior_box, axis));
         for (pdat::FaceIterator ci(pdat::FaceGeometry::begin(interior_box, axis));
              ci != ciend; ++ci) {
            for (int d = 0; d < depth; ++d) {
               FACE_MBLK_KERNEL_TYPE result = (*face_data)(*ci, d);

               if (!tbox::MathUtilities<FACE_MBLK_KERNEL_TYPE>::equalEps(correct, result)) {
                  tbox::perr << "Test FAILED: ...."
                             << " : face index = " << *ci << std::endl;
                  tbox::perr << "    Variable = " << d_variable_src_name[i]
                             << " : depth index = " << d << std::endl;
                  tbox::perr << "    result = " << result
                             << " : correct = " << correct << std::endl;
                  test_failed = true;
               }
            }
         }
      }

      hier::Box gbox = face_data->getGhostBox();

      for (tbox::Dimension::dir_t axis = 0; axis < d_dim.getValue(); ++axis) {
         hier::Box patch_face_box =
            pdat::FaceGeometry::toFaceBox(pbox, axis);

         hier::BoxContainer tested_neighbors;

         for (hier::BaseGridGeometry::ConstNeighborIterator ne(
                 grid_geom->begin(block_id));
              ne != grid_geom->end(block_id); ++ne) {

            const hier::BaseGridGeometry::Neighbor& nbr = *ne;
            if (nbr.isSingularity()) {
               continue;
            }

            correct = nbr.getBlockId().getBlockValue();

            hier::BoxContainer neighbor_ghost(nbr.getTransformedDomain());
            hier::BoxContainer neighbor_face_ghost;
            for (hier::BoxContainer::iterator nn = neighbor_ghost.begin();
                 nn != neighbor_ghost.end(); ++nn) {
               hier::Box neighbor_ghost_interior(
                  pdat::FaceGeometry::toFaceBox(*nn, axis));
               neighbor_ghost_interior.grow(-hier::IntVector::getOne(d_dim));
               neighbor_face_ghost.pushFront(neighbor_ghost_interior);
            }

            neighbor_face_ghost.refine(ratio);
            neighbor_face_ghost.intersectBoxes(
               pdat::FaceGeometry::toFaceBox(gbox, axis));

            neighbor_face_ghost.removeIntersections(tested_neighbors);

            for (hier::BoxContainer::iterator ng = neighbor_face_ghost.begin();
                 ng != neighbor_face_ghost.end(); ++ng) {

               hier::Box::iterator ciend(ng->end());
               for (hier::Box::iterator ci(ng->begin()); ci != ciend; ++ci) {
                  pdat::FaceIndex fi(*ci, 0, 0);
                  fi.setAxis(axis);
                  if (!patch_face_box.contains(fi)) {
                     for (int d = 0; d < depth; ++d) {
                        FACE_MBLK_KERNEL_TYPE result = (*face_data)(fi, d);

                        if (!tbox::MathUtilities<FACE_MBLK_KERNEL_TYPE>::equalEps(correct,
                               result)) {
                           tbox::perr << "Test FAILED: ...."
                                      << " : face index = " << fi << std::endl;
                           tbox::perr << "  Variable = "
                                      << d_variable_src_name[i]
                                      << " : depth index = " << d << std::endl;
                           tbox::perr << "    result = " << result
                                      << " : correct = " << correct << std::endl;
                           test_failed = true;
                        }
                     }
                  }
               }
            }
            tested_neighbors.spliceBack(neighbor_face_ghost);
         }
      }

      std::shared_ptr<hier::PatchGeometry> pgeom(patch.getPatchGeometry());

      for (int b = 0; b < d_dim.getValue(); ++b) {
         const std::vector<hier::BoundaryBox> bdry =
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

                  correct = (FACE_MBLK_KERNEL_TYPE)bdry[k].getLocationIndex() + 200.0;

               } else {

                  correct /= (FACE_MBLK_KERNEL_TYPE)num_sing_neighbors;

               }

            } else {
               correct = (FACE_MBLK_KERNEL_TYPE)(bdry[k].getLocationIndex() + 100);
            }

            for (tbox::Dimension::dir_t axis = 0; axis < d_dim.getValue(); ++axis) {
               hier::Box patch_face_box =
                  pdat::FaceGeometry::toFaceBox(pbox, axis);

               pdat::FaceIterator ciend(pdat::FaceGeometry::end(fill_box, axis));
               for (pdat::FaceIterator ci(pdat::FaceGeometry::begin(fill_box, axis));
                    ci != ciend; ++ci) {

                  if (!patch_face_box.contains(*ci)) {

                     bool use_index = true;
                     for (tbox::Dimension::dir_t n = 0; n < d_dim.getValue(); ++n) {
                        if (axis == n && bdry[k].getBox().numberCells(n) ==
                            1) {
                           if ((*ci)(0) == patch_face_box.lower() (0) ||
                               (*ci)(0) == patch_face_box.upper() (0)) {
                              use_index = false;
                              break;
                           }
                        }
                     }

                     if (use_index) {
                        for (int d = 0; d < depth; ++d) {
                           FACE_MBLK_KERNEL_TYPE result = (*face_data)(*ci, d);

                           if (!tbox::MathUtilities<FACE_MBLK_KERNEL_TYPE>::equalEps(correct,
                                  result)) {
                              tbox::perr << "Test FAILED: ...."
                                         << " : face index = " << *ci << std::endl;
                              tbox::perr << "  Variable = "
                                         << d_variable_src_name[i]
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
         }
      }
   }

   if (!test_failed) {
      tbox::plog << "FaceMultiblockTest Successful!" << std::endl;
   } else {
      tbox::perr << "Multiblock FaceMultiblockTest FAILED: .\n" << std::endl;
   }

   solution.reset();   // just to be anal...

   tbox::plog << "\nExiting FaceMultiblockTest::verifyResults..." << std::endl;
   tbox::plog << "level_number = " << level_number << std::endl;
   tbox::plog << "Patch box = " << patch.getBox() << std::endl << std::endl;

   return !test_failed;
}

void FaceMultiblockTest::postprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const std::shared_ptr<hier::VariableContext>& context,
   const hier::Box& fine_box,
   const hier::IntVector& ratio) const
{
   pdat::FaceDoubleConstantRefine ref_op;

   hier::BoxContainer fine_box_list(fine_box);
   hier::BoxContainer empty_box_list;

   xfer::BoxGeometryVariableFillPattern fill_pattern;

   for (int i = 0; i < static_cast<int>(d_variables.size()); ++i) {

      int id = hier::VariableDatabase::getDatabase()->
         mapVariableAndContextToIndex(d_variables[i], context);

      std::shared_ptr<hier::PatchDataFactory> fine_pdf(
         fine.getPatchDescriptor()->getPatchDataFactory(id));

      std::shared_ptr<hier::BoxOverlap> fine_overlap(
         fill_pattern.computeFillBoxesOverlap(
            fine_box_list,
            empty_box_list,
            fine.getBox(),
            fine.getPatchData(id)->getGhostBox(),
            *fine_pdf));

      ref_op.refine(fine, coarse, id, id, *fine_overlap, ratio);
   }
}

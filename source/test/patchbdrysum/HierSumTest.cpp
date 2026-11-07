/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   SAMRAI interface class for hierarchy node and edge sum test
 *
 ************************************************************************/

#include "HierSumTest.h"

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/EdgeData.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeIndex.h"
#include "SAMRAI/pdat/NodeIterator.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/algs/PatchBoundaryNodeSum.h"
#include "SAMRAI/algs/PatchBoundaryEdgeSum.h"
#include "SAMRAI/hier/CoarseFineBoundary.h"
#include "SAMRAI/tbox/Collectives.h"
#include "SAMRAI/tbox/NVTXUtilities.h"

extern "C" {
void SAMRAI_F77_FUNC(setedges2d, SETEDGES2D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double *,
   double *,
   double *);

void SAMRAI_F77_FUNC(checkedges2d, CHECKEDGES2D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double&,
   int&,
   const double *,
   const double *);
void SAMRAI_F77_FUNC(setedges3d, SETEDGES3D) (const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const double *,
   double *,
   double *,
   double *);

void SAMRAI_F77_FUNC(checkedges3d, CHECKEDGES3D) (const int&, const int&,
   const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const double&,
   int&,
   const double *,
   const double *,
   const double *);
}

using namespace SAMRAI;
using namespace tbox;
using namespace hier;
using namespace pdat;
using namespace geom;
using namespace algs;

/*************************************************************************
 *
 * Constructor and Destructor.
 *
 ************************************************************************/

HierSumTest::HierSumTest(
   const std::string& object_name,
   const tbox::Dimension& dim,
   std::shared_ptr<Database> input_db
#ifdef HAVE_HDF5
   ,
   std::shared_ptr<appu::VisItDataWriter> viz_writer
#endif
   ):
   d_object_name(object_name),
   d_dim(dim),
   d_depth(1),
   d_ucell_var(new CellVariable<double>(dim, "ucell", d_depth)),
   d_unode_var(new NodeVariable<double>(dim, "unode", d_depth)),
   d_uedge_var(new EdgeVariable<double>(dim, "uedge", d_depth)),
   d_node_ghosts(dim),
   d_edge_ghosts(dim)
{
   /*
    * Initialize object with data read from given input databases.
    */
   getFromInput(input_db);

   /*
    * Set up variables and contexts.
    *
    * Vars:
    *    ucell - cell-centered variable U
    *    unode - node-centered variable U
    *    uedge - edge-centered variable U
    *
    * Contexts:
    *
    *    NOGHOST  - zero ghosts (unode)
    *    ONEGHOST - one ghost (ucell)
    *
    * What the test case does:
    *   1. Set ucell = 1.0 on cells of all levels
    *   2. Set ucell = 0.0 on cells of L < LN that are
    *      covered by refined cells.
    *   3. Set node weight on patch interiors = sum(cell weights)
    *   3. Do a hier sum transaction
    *   4. Correct result - all nodes on all levels = 2^d_dim
    *
    *
    * Below we construct u variable and its contexts.
    */
   VariableDatabase* variable_db = VariableDatabase::getDatabase();

   std::shared_ptr<VariableContext> cxt1(
      variable_db->getContext("CONTEXT1"));
   std::shared_ptr<VariableContext> cxt2(
      variable_db->getContext("CONTEXT2"));

   IntVector one_ghost(dim, 1);

   d_ucell_node_id =
      variable_db->registerVariableAndContext(d_ucell_var,
         cxt1,
         one_ghost);

   d_unode_id = variable_db->registerVariableAndContext(d_unode_var,
         cxt1,
         d_node_ghosts);
   d_ucell_edge_id =
      variable_db->registerVariableAndContext(d_ucell_var,
         cxt2,
         one_ghost);

   d_uedge_id = variable_db->registerVariableAndContext(d_uedge_var,
         cxt1,
         d_edge_ghosts);

#ifdef HAVE_HDF5
   /*
    * Register u values to be written by the viz writer.
    */
   if (viz_writer) {
      viz_writer->registerPlotQuantity("ucell::node", "SCALAR",
         d_ucell_node_id, 0);
      viz_writer->registerPlotQuantity("ucell::edge", "SCALAR",
         d_ucell_edge_id, 0);
      viz_writer->registerPlotQuantity("unode", "SCALAR", d_unode_id);
   }
#endif

}

HierSumTest::~HierSumTest()
{
}

/*************************************************************************
 *
 * Set initial node values based on sum of surrounding cell weights
 *
 ************************************************************************/
int
HierSumTest::setInitialNodeValues(
   const std::shared_ptr<PatchHierarchy> hierarchy)
{
   int fail_count = 0;

   /*
    * Set node weight on patch interiors = sum(cell weights)
    */
   // loop over hierarchy levels
   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln) {
      std::shared_ptr<PatchLevel> level(hierarchy->getPatchLevel(ln));

      // loop over patches on level
      for (PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<Patch>& patch = *ip;

         std::shared_ptr<NodeData<double> > unode(
            SAMRAI_SHARED_PTR_CAST<NodeData<double>, PatchData>(
               patch->getPatchData(d_unode_id)));
         std::shared_ptr<CellData<double> > ucell(
            SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
               patch->getPatchData(d_ucell_node_id)));
         TBOX_ASSERT(unode);
         TBOX_ASSERT(ucell);

         // output initial cell values
         int level_number = level->getLevelNumber();
         tbox::plog << "INITIAL Cell values for NODE - Level: " << level_number
                    << "\tPatch: " << patch->getBox() << std::endl;
         ucell->print(ucell->getGhostBox(), plog);

         // loop over nodes of patch
         Box pbox = patch->getBox();
         double cell_val;
         int d;
         NodeIterator niend(NodeGeometry::end(pbox));
         for (NodeIterator ni(NodeGeometry::begin(pbox)); ni != niend; ++ni) {
            NodeIndex node = *ni;
            for (d = 0; d < ucell->getDepth(); ++d) {

               (*unode)(node, d) = 0.;

               /*
                * Sum contributions from surrounding cells for each node
                * value.
                */
               if (d_dim == tbox::Dimension(2)) {
                  for (int j = 0; j <= 1; ++j) {
                     for (int i = 0; i <= 1; ++i) {
                        CellIndex cell(*ni);
                        cell(0) -= i;
                        cell(1) -= j;
                        cell_val = (*ucell)(cell, d);
                        (*unode)(node, d) += cell_val;
                     }
                  }
               }
               if (d_dim == tbox::Dimension(3)) {
                  for (int k = 0; k <= 1; ++k) {
                     for (int j = 0; j <= 1; ++j) {
                        for (int i = 0; i <= 1; ++i) {
                           CellIndex cell(*ni);
                           cell(0) -= i;
                           cell(1) -= j;
                           cell(2) -= k;
                           cell_val = (*ucell)(cell, d);
                           (*unode)(node, d) += cell_val;
                        }
                     }
                  }
               }
            } // loop over depth
         } // loop over nodes
      } // loop over patches

      /*
       * Any nodes that are *inside* the complement region (inside
       * meaning at least one cell away from the coarse-fine boundary)
       * should not contribute to the nodal sum on finer levels.  Here
       * we verify this is the case by resetting all nodes inside this
       * region to -999.
       *
       * The so-called fine_overlap_shrunk is computed by finding the
       * region of the level that is not overlapped by fine
       * patches (complement), growing it, and removing intersection
       * with the level boxes.  This is effectively a shrunken coarse-fine
       * overlap region, on which the nodes shouldn't participate in any
       * communication.
       */
      BoxContainer fine_overlap_shrunk = level->getBoxes();
      BoxContainer complement(fine_overlap_shrunk);
      if (level->getLevelNumber() != hierarchy->getFinestLevelNumber()) {
         std::shared_ptr<PatchLevel> fine_level(
            hierarchy->getPatchLevel(ln + 1));
         BoxContainer fine_level_boxes = fine_level->getBoxes();
         fine_level_boxes.coarsen(fine_level->getRatioToCoarserLevel());
         complement.removeIntersections(fine_level_boxes);
         complement.grow(IntVector::getOne(d_dim));
         fine_overlap_shrunk.removeIntersections(complement);

         for (PatchLevel::iterator ip(level->begin());
              ip != level->end(); ++ip) {
            const std::shared_ptr<Patch>& patch = *ip;

            std::shared_ptr<NodeData<double> > unode(
               SAMRAI_SHARED_PTR_CAST<NodeData<double>, PatchData>(
                  patch->getPatchData(d_unode_id)));
            TBOX_ASSERT(unode);

            for (BoxContainer::iterator b = fine_overlap_shrunk.begin();
                 b != fine_overlap_shrunk.end(); ++b) {
               Box fine_overlap = *b;
               Box patch_interior = patch->getBox();
               Box data_box = fine_overlap * patch_interior;
               NodeIterator niend(NodeGeometry::end(data_box));
               for (NodeIterator ni(NodeGeometry::begin(data_box));
                    ni != niend; ++ni) {
                  NodeIndex node = *ni;
                  for (int d = 0; d < unode->getDepth(); ++d) {
                     double node_val = (*unode)(node, d);
                     if (tbox::MathUtilities<double>::equalEps(node_val,
                            0.0)) {
                        (*unode)(node, d) = -999.;
                     }
                  }
               }
            } // loop over complement boxes
         } // loop over patches
      } // if a finer level exists

      for (PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<Patch>& patch = *ip;

         std::shared_ptr<NodeData<double> > unode(
            SAMRAI_SHARED_PTR_CAST<NodeData<double>, PatchData>(
               patch->getPatchData(d_unode_id)));
         TBOX_ASSERT(unode);

         // output initial node values
         tbox::plog << "INITIAL Node values - Level: " << level->getLevelNumber()
                    << "\tPatch: " << patch->getBox() << std::endl;
         unode->print(unode->getGhostBox(), plog);

      } // loop over patches
   } // loop over levels

   return fail_count;
}

/*************************************************************************
 *
 * Set initial edge values based on sum of surrounding cell weights
 *
 ************************************************************************/
int
HierSumTest::setInitialEdgeValues(
   const std::shared_ptr<PatchLevel> level)
{
   int fail_count = 0;

   double correct_val = 1.;
   for (int i = 0; i < d_dim.getValue() - 1; ++i) {
      correct_val *= 2.;
   }

   /*
    * Set edge weight on patch interiors = sum(cell weights)
    */
   for (PatchLevel::iterator ip(level->begin());
        ip != level->end(); ++ip) {
      const std::shared_ptr<Patch>& patch = *ip;

      std::shared_ptr<EdgeData<double> > uedge(
         SAMRAI_SHARED_PTR_CAST<EdgeData<double>, PatchData>(
            patch->getPatchData(d_uedge_id)));
      std::shared_ptr<CellData<double> > ucell(
         SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
            patch->getPatchData(d_ucell_edge_id)));
      TBOX_ASSERT(uedge);
      TBOX_ASSERT(ucell);

      // output initial cell values
      int level_number = level->getLevelNumber();
      tbox::plog << "INITIAL Cell values for EDGE - Level: " << level_number
                 << "\tPatch: " << patch->getBox() << std::endl;
      ucell->print(ucell->getGhostBox(), plog);

      const Index ifirst(patch->getBox().lower());
      const Index ilast(patch->getBox().upper());

      IntVector cellg(ucell->getGhostCellWidth());
      IntVector edgeg(uedge->getGhostCellWidth());

      for (int d = 0; d < uedge->getDepth(); ++d) {

         if (d_dim == tbox::Dimension(2)) {
            SAMRAI_F77_FUNC(setedges2d, SETEDGES2D) (ifirst(0), ifirst(1),
               ilast(0), ilast(1),
               cellg(0), cellg(1),
               edgeg(0), edgeg(1),
               ucell->getPointer(d),
               uedge->getPointer(0, d),
               uedge->getPointer(1, d));
         }
         if (d_dim == tbox::Dimension(3)) {
            SAMRAI_F77_FUNC(setedges3d, SETEDGES3D) (ifirst(0), ifirst(1), ifirst(2),
               ilast(0), ilast(1), ilast(2),
               cellg(0), cellg(1), cellg(2),
               edgeg(0), edgeg(1), edgeg(2),
               ucell->getPointer(d),
               uedge->getPointer(0, d),
               uedge->getPointer(1, d),
               uedge->getPointer(2, d));
         }

         /*
          * If you want to check edges BEFORE communication (to make sure
          * the communication is actually doing something) then do the
          * check here.  Be forwarned that it dumps a lot of data because,
          * if things are working right, the data before communication will
          * have many errors, that the communication will fix.
          */

         if (d_check_data_before_communication) {

            int fort_all_correct = 1;

            if (d_dim == tbox::Dimension(2)) {
               SAMRAI_F77_FUNC(checkedges2d, CHECKEDGES2D) (ifirst(0), ifirst(1),
                  ilast(0), ilast(1),
                  edgeg(0), edgeg(1),
                  correct_val,
                  fort_all_correct,
                  uedge->getPointer(0, d),
                  uedge->getPointer(1, d));
            }
            if (d_dim == tbox::Dimension(3)) {
               SAMRAI_F77_FUNC(checkedges3d, CHECKEDGES3D) (ifirst(0), ifirst(1),
                  ifirst(2),
                  ilast(0), ilast(1), ilast(2),
                  edgeg(0), edgeg(1), edgeg(2),
                  correct_val,
                  fort_all_correct,
                  uedge->getPointer(0, d),
                  uedge->getPointer(1, d),
                  uedge->getPointer(2, d));
            }

            if (fort_all_correct == 0) {
               ++fail_count;
               tbox::perr
               << "PatchBdrySum Edge test FAILED:  Errors on Level: "
               << level->getLevelNumber()
               << "\t Patch: " << patch->getBox()
               << "\nAll edges are not correct value." << std::endl;
            } else {
#if (TESTING == 1)
               tbox::plog
#else
               tbox::pout
#endif
               << "All edges on Level: " << level->getLevelNumber()
               << "\t Patch: " << patch->getBox()
               << "\tare correct." << std::endl;
            }
         }

      } // loop over depth

#if (TESTING == 1)
      tbox::plog << "INITIAL Edge values - Level: " << level->getLevelNumber()
                 << "\tPatch: " << patch->getBox() << std::endl;
      uedge->print(uedge->getGhostBox(), plog);
#endif

   } // loop over patches

   tbox::SAMRAI_MPI::getSAMRAIWorld().Barrier();

   return fail_count;

}

/*************************************************************************
 *
 * Setup outer node sum.
 *
 ************************************************************************/
void
HierSumTest::setupOuternodeSum(
   const std::shared_ptr<PatchHierarchy> hierarchy)
{
   d_node_sum_util.reset(new PatchBoundaryNodeSum("Node Sum Util"));

   d_node_sum_util->registerSum(d_unode_id);

   int num_levels = hierarchy->getNumberOfLevels();
   if (num_levels > 1) {

      int coarsest_level_number = 0;
      d_node_sum_util->setupSum(hierarchy,
         coarsest_level_number,
         hierarchy->getFinestLevelNumber());
   } else {
      d_node_sum_util->setupSum(hierarchy->getPatchLevel(0));
   }

}

/*************************************************************************
 *
 * Perform outer node sum.
 *
 ************************************************************************/
void
HierSumTest::doOuternodeSum()
{
   TBOX_ASSERT(d_node_sum_util);

   bool fill_hanging_nodes = true;
   d_node_sum_util->computeSum(fill_hanging_nodes);

}

/*************************************************************************
 *
 * Setup Outeredge sum.
 *
 ************************************************************************/
void
HierSumTest::setupOuteredgeSum(
   const std::shared_ptr<PatchHierarchy> hierarchy,
   const int level_num)
{
   if (level_num >= static_cast<int>(d_edge_sum_util.size())) {
      d_edge_sum_util.resize(level_num + 1);
   }

   d_edge_sum_util[level_num].reset(
      new PatchBoundaryEdgeSum("Level Edge Sum Util"));

   d_edge_sum_util[level_num]->registerSum(d_uedge_id);

   d_edge_sum_util[level_num]->setupSum(hierarchy->getPatchLevel(level_num));
}

/*************************************************************************
 *
 * Setup and perform Outeredge sum.
 *
 ************************************************************************/
void
HierSumTest::doOuteredgeSum(
   const int level_num)
{
   TBOX_ASSERT(level_num < static_cast<int>(d_edge_sum_util.size()));
   TBOX_ASSERT(d_edge_sum_util[level_num]);

   d_edge_sum_util[level_num]->computeSum();

}

/*************************************************************************
 *
 * Check correctness of hierarchy node sum operation
 *
 ************************************************************************/

int HierSumTest::checkNodeResult(
   const std::shared_ptr<PatchHierarchy> hierarchy)
{

   int fail_count = 0;

   /*
    * After the communication the sum on every node of every level should be
    * 2^d_dim.  Check this here...
    */

   double correct_val = 1.;
   for (int i = 0; i < d_dim.getValue(); ++i) {
      correct_val *= 2.;
   }

   // loop over hierarchy levels
   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln) {
      std::shared_ptr<PatchLevel> level(hierarchy->getPatchLevel(ln));

      BoxContainer level_boxes_complement = level->getBoxes();

      // If a finer level exists, remove overlap boxes by computing complement
      if (level->getLevelNumber() != hierarchy->getFinestLevelNumber()) {
         std::shared_ptr<PatchLevel> fine_level(
            hierarchy->getPatchLevel(ln + 1));
         BoxContainer fine_level_boxes = fine_level->getBoxes();
         fine_level_boxes.coarsen(fine_level->getRatioToCoarserLevel());
         level_boxes_complement.removeIntersections(fine_level_boxes);
      }

      // loop over patches on level
      for (PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         bool all_correct = true;
         const std::shared_ptr<Patch>& patch = *ip;

         std::shared_ptr<NodeData<double> > unode(
            SAMRAI_SHARED_PTR_CAST<NodeData<double>, PatchData>(
               patch->getPatchData(d_unode_id)));
         TBOX_ASSERT(unode);

         // loop over Level complement boxlist
         for (BoxContainer::iterator b = level_boxes_complement.begin();
              b != level_boxes_complement.end(); ++b) {
            Box complement = *b;

            // intersect patch box with level box complement
            Box patch_interior = patch->getBox();
            Box data_box = complement * patch_interior;

            /*
             * Iterate over nodes and check correctness of result.
             */
            NodeIterator iend(NodeGeometry::end(data_box));
            for (NodeIterator i(NodeGeometry::begin(data_box));
                 i != iend; ++i) {
               NodeIndex node = *i;  // i,j
               for (int d = 0; d < unode->getDepth(); ++d) {

                  bool node_correct = false;
                  double node_val = (*unode)(node, d);

                  if (tbox::MathUtilities<double>::equalEps(node_val,
                         correct_val)) {
                     node_correct = true;
                  }

                  if (!node_correct) {
                     tbox::perr << "BAD NODE = " << node_val << " at index "
                                << *i
                                << " in L" << ln << " " << patch->getBox()
                                << " depth = " << d << " should be "
                                << correct_val << std::endl;
                     all_correct = false;
                     break;
                  }
               }
            }
         } // loop over complement boxes

         if (!all_correct) {
            ++fail_count;
            tbox::perr << "PatchBdrySum Node test FAILED:  Errors on Level: "
                       << level->getLevelNumber()
                       << "\t Patch: " << patch->getBox()
                       << "\nAll nodes are not correct value." << std::endl;
         } else {
#if (TESTING == 1)
            tbox::plog
#else
            tbox::pout
#endif
            << "All nodes on Level: " << level->getLevelNumber()
            << "\t Patch: " << patch->getBox()
            << "\tare correct." << std::endl;
         }

#if (TESTING == 1)
         std::shared_ptr<CellData<double> > ucell_node(
            SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
               patch->getPatchData(d_ucell_node_id)));
         TBOX_ASSERT(ucell_node);

         tbox::plog << "FINAL Cell values for NODE - Level: "
                    << level->getLevelNumber()
                    << "\tPatch: " << patch->getBox() << std::endl;
         ucell_node->print(ucell_node->getGhostBox(), plog);

         tbox::plog << "FINAL Node values - Level: " << level->getLevelNumber()
                    << "\tPatch " << patch->getBox()
                    << std::endl;
         unode->print(unode->getBox(), plog);
#endif

      } // loop over patches

   } // loop over levels

   return fail_count;

}

/*************************************************************************
 *
 * Check correctness of level edge sum operation
 *
 ************************************************************************/

int HierSumTest::checkEdgeResult(
   const std::shared_ptr<PatchLevel> level)
{

   int fail_count = 0;

   /*
    * After the communication the sum on every edge of every level should be
    * 2^(d_dim-1).  Check this here...
    */

   double correct_val = 1.;
   for (int i = 0; i < d_dim.getValue() - 1; ++i) {
      correct_val *= 2.;
   }

   // loop over patches on level
   for (PatchLevel::iterator ip(level->begin());
        ip != level->end(); ++ip) {
      const std::shared_ptr<Patch>& patch = *ip;

      std::shared_ptr<EdgeData<double> > uedge(
         SAMRAI_SHARED_PTR_CAST<EdgeData<double>, PatchData>(
            patch->getPatchData(d_uedge_id)));
      TBOX_ASSERT(uedge);

      const Index ifirst(patch->getBox().lower());
      const Index ilast(patch->getBox().upper());

      IntVector edgeg(uedge->getGhostCellWidth());

      for (int d = 0; d < uedge->getDepth(); ++d) {

         /*
          * In the fortran, set "fort_all_correct" to 0 if we
          * detect differences between uedge data and "correct_val".
          */
         int fort_all_correct = 1;

         if (d_dim == tbox::Dimension(2)) {
            SAMRAI_F77_FUNC(checkedges2d, CHECKEDGES2D) (ifirst(0), ifirst(1),
               ilast(0), ilast(1),
               edgeg(0), edgeg(1),
               correct_val,
               fort_all_correct,
               uedge->getPointer(0, d),
               uedge->getPointer(1, d));
         }
         if (d_dim == tbox::Dimension(3)) {
            SAMRAI_F77_FUNC(checkedges3d, CHECKEDGES3D) (ifirst(0), ifirst(1),
               ifirst(2),
               ilast(0), ilast(1), ilast(2),
               edgeg(0), edgeg(1), edgeg(2),
               correct_val,
               fort_all_correct,
               uedge->getPointer(0, d),
               uedge->getPointer(1, d),
               uedge->getPointer(2, d));
         }

         if (fort_all_correct == 0) {
            ++fail_count;
            tbox::perr << "PatchBdrySum Edge test FAILED:  Errors on Level: "
                       << level->getLevelNumber()
                       << "\t Patch: " << patch->getBox()
                       << "\nAll edges are not correct value." << std::endl;
         } else {
#if (TESTING == 1)
            tbox::plog
#else
            tbox::pout
#endif
            << "All edges on Level: " << level->getLevelNumber()
            << "\t Patch: " << patch->getBox()
            << "\tare correct." << std::endl;
         }

      } // loop over depth

#if (TESTING == 1)
      std::shared_ptr<CellData<double> > ucell_edge(
         SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
            patch->getPatchData(d_ucell_edge_id)));
      TBOX_ASSERT(ucell_edge);

      tbox::plog << "FINAL Cell values for EDGE - Level: "
                 << level->getLevelNumber()
                 << "\tPatch: " << patch->getBox() << std::endl;
      ucell_edge->print(ucell_edge->getGhostBox(), plog);

      tbox::plog << "FINAL Edge values - Level: "
                 << level->getLevelNumber()
                 << "\tPatch: " << patch->getBox() << std::endl;
      uedge->print(uedge->getGhostBox(), plog);
#endif

   } // loop over patches

   return fail_count;

}

/*************************************************************************
 *
 * Methods inherited from StandardTagAndInitStrategy.
 *
 ************************************************************************/

/*
 * Allocate storage and initialize data on the level.
 */

void HierSumTest::initializeLevelData(
   const std::shared_ptr<PatchHierarchy>& hierarchy,
   const int level_number,
   const double time,
   const bool can_be_refined,
   const bool initial_time,
   const std::shared_ptr<PatchLevel>& old_level,
   const bool allocate_data)
{
   NULL_USE(can_be_refined);
   NULL_USE(initial_time);
   NULL_USE(old_level);

   std::shared_ptr<PatchHierarchy> local_hierarchy(hierarchy);

   /*
    * Set initial data on hierarchy level.
    *   1. Set ucell = 1.0 on cells of all levels
    *   2. For NODE data only, set ucell = 0.0 on cells of L < LN that are
    *      covered by refined cells.
    */
   std::shared_ptr<PatchLevel> level(
      hierarchy->getPatchLevel(level_number));

   /*
    * Allocate storage for cell and node data.
    */
   if (allocate_data) {
      level->allocatePatchData(d_ucell_node_id, time);
      level->allocatePatchData(d_ucell_edge_id, time);
      level->allocatePatchData(d_unode_id, time);
      level->allocatePatchData(d_uedge_id, time);
   }

   /*
    * Set edge/node values to zero initially.
    */
   for (PatchLevel::iterator p0(level->begin());
        p0 != level->end(); ++p0) {
      const std::shared_ptr<Patch>& patch = *p0;

      std::shared_ptr<NodeData<double> > unode(
         SAMRAI_SHARED_PTR_CAST<NodeData<double>, PatchData>(
            patch->getPatchData(d_unode_id)));
      std::shared_ptr<EdgeData<double> > uedge(
         SAMRAI_SHARED_PTR_CAST<EdgeData<double>, PatchData>(
            patch->getPatchData(d_uedge_id)));
      TBOX_ASSERT(unode);
      TBOX_ASSERT(uedge);
      unode->fillAll(0.0);
      uedge->fillAll(0.0);
   }

#if defined(HAVE_RAJA)
   tbox::parallel_synchronize();
#endif
   /*
    * Set cell weights.  We want interior cells set to 1.0 and
    * ghost cells set to 0.0.  (Eventually, we will reset the ghosts
    * based on overlap with neighboring patches but this will be
    * in the next step).
    */
   for (PatchLevel::iterator p0(level->begin());
        p0 != level->end(); ++p0) {
      const std::shared_ptr<Patch>& patch = *p0;

      std::shared_ptr<CellData<double> > ucell_node(
         SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
            patch->getPatchData(d_ucell_node_id)));
      std::shared_ptr<CellData<double> > ucell_edge(
         SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
            patch->getPatchData(d_ucell_edge_id)));
      TBOX_ASSERT(ucell_node);
      TBOX_ASSERT(ucell_edge);

      ucell_node->fillAll(0.0, ucell_node->getGhostBox()); // ghost box
      ucell_node->fillAll(1.0, patch->getBox());          // interior patch box

      ucell_edge->fillAll(0.0, ucell_edge->getGhostBox()); // ghost box
      ucell_edge->fillAll(1.0, patch->getBox());          // interior patch box
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif


      // set cell values at physical boundary
      const std::shared_ptr<CartesianPatchGeometry> patch_geom(
         SAMRAI_SHARED_PTR_CAST<CartesianPatchGeometry, PatchGeometry>(
            patch->getPatchGeometry()));
      TBOX_ASSERT(patch_geom);
      const std::vector<BoundaryBox>& node_bdry =
         patch_geom->getCodimensionBoundaries(d_dim.getValue());
      const std::vector<BoundaryBox>& edge_bdry =
         patch_geom->getCodimensionBoundaries(d_dim.getValue() - 1);
      std::vector<BoundaryBox> empty_vector(0, BoundaryBox(d_dim));
      const std::vector<BoundaryBox>& face_bdry =
         d_dim == tbox::Dimension(3) ?
         patch_geom->getCodimensionBoundaries(1) : empty_vector;
      // node cell values
      setBoundaryConditions(*patch,
         node_bdry,
         edge_bdry,
         face_bdry,
         d_ucell_node_id);

      // edge cell values
      setBoundaryConditions(*patch,
         node_bdry,
         edge_bdry,
         face_bdry,
         d_ucell_edge_id);

   }

   if (level_number > 0) {

      /*
       * For node data, set the cell weights to zero on coarser level
       * where there is overlap with fine level patches.
       */
      std::shared_ptr<PatchLevel> coarser_level(
         hierarchy->getPatchLevel(level_number - 1));
      BoxContainer fine_level_boxes = level->getBoxes();

      IntVector ratio(level->getRatioToCoarserLevel());
      fine_level_boxes.coarsen(ratio);

      for (PatchLevel::iterator p1(coarser_level->begin());
           p1 != coarser_level->end(); ++p1) {
         const std::shared_ptr<Patch>& cpatch = *p1;

         std::shared_ptr<CellData<double> > ucell_node(
            SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
               cpatch->getPatchData(d_ucell_node_id)));
         TBOX_ASSERT(ucell_node);

         Box cpbox = cpatch->getBox();
         for (BoxContainer::iterator fine_level_itr = fine_level_boxes.begin();
              fine_level_itr != fine_level_boxes.end(); ++fine_level_itr) {
            Box setbox = cpbox * *fine_level_itr;
            if (!setbox.empty()) {
               ucell_node->fillAll(0.0, setbox);
            }
         }
#if defined(HAVE_RAJA)
         tbox::parallel_synchronize();
#endif


         // zero out cells on the boundary that lie at coarse-fine
         // interface.
         zeroOutPhysicalBoundaryCellsAtCoarseFineBoundary(*cpatch,
            d_ucell_node_id);

      } // loop over coarser level patches

      /*
       * For edge data, set the ghosts of the cell data equal to 1.0 at the
       * coarse-fine boundaries.
       */
      IntVector max_ghosts(d_dim, 1);
      CoarseFineBoundary cfbdry(*hierarchy,
                                level_number,
                                max_ghosts);

      level = hierarchy->getPatchLevel(level_number);
      for (PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {
         const std::shared_ptr<Patch>& patch = *p;
         Box pbox = patch->getBox();
         const GlobalId global_id = patch->getGlobalId();

         const std::vector<BoundaryBox>& node_bdry =
            cfbdry.getNodeBoundaries(global_id);
         const std::vector<BoundaryBox>& edge_bdry =
            cfbdry.getEdgeBoundaries(global_id);
         std::vector<BoundaryBox> empty_vector(0, BoundaryBox(d_dim));
         const std::vector<BoundaryBox>& face_bdry =
            d_dim == tbox::Dimension(3) ?
            cfbdry.getFaceBoundaries(global_id) : empty_vector;

         setBoundaryConditions(*patch,
            node_bdry,
            edge_bdry,
            face_bdry,
            d_ucell_edge_id);

      } // loop over level patches

   } // if level_number > 0

}

/*
 * Perform operations necessary when grid changes for dynamic grid
 * calculations (to be added later...)
 */
void
HierSumTest::resetHierarchyConfiguration(
   const std::shared_ptr<PatchHierarchy>& hierarchy,
   const int coarsest_level,
   const int finest_level)
{
   NULL_USE(hierarchy);
   NULL_USE(coarsest_level);
   NULL_USE(finest_level);
}

/*
 * Routine to do cell tagging.  Add later...
 */
void
HierSumTest::applyGradientDetector(
   const std::shared_ptr<PatchHierarchy>& hierarchy,
   const int level_number,
   const double time,
   const int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(hierarchy);
   NULL_USE(level_number);
   NULL_USE(time);
   NULL_USE(tag_index);
   NULL_USE(initial_time);
   NULL_USE(uses_richardson_extrapolation_too);
}

/*
 * Set boundary conditions by shifting patch appropriately and
 * finding intersection with boundary fill box.
 */
void
HierSumTest::setBoundaryConditions(
   Patch& patch,
   const std::vector<BoundaryBox>& node_bdry,
   const std::vector<BoundaryBox>& edge_bdry,
   const std::vector<BoundaryBox>& face_bdry,
   const int cell_data_id)
{
   const int num_node_bdry_boxes = static_cast<int>(node_bdry.size());
   const int num_edge_bdry_boxes = static_cast<int>(edge_bdry.size());
   const int num_face_bdry_boxes = static_cast<int>(face_bdry.size());

   const std::shared_ptr<CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<CartesianPatchGeometry, PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);

   /*
    * std::shared_ptr to data in ghost regions.
    */
   std::shared_ptr<CellData<double> > ucell(
      SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
         patch.getPatchData(cell_data_id)));

   TBOX_ASSERT(ucell);

   IntVector ghost_cells(ucell->getGhostCellWidth());
   const Box pbox(patch.getBox());

   int i, d;
   IntVector shift(d_dim, 0);
   Box shifted_pbox = pbox;

   if (d_dim == tbox::Dimension(3)) {
      /*
       * Set cell weights to 1.0 on FACES of patch.
       */
      for (i = 0; i < num_face_bdry_boxes; ++i) {
         Box fill_box = patch_geom->getBoundaryFillBox(face_bdry[i],
               pbox,
               ghost_cells);
         /*
          * location index:
          *    0,1 - X lower,upper
          *    2,3 - Y lower,upper
          *    4,5 - Z lower,upper
          */
         int loc_indx = face_bdry[i].getLocationIndex();
         shift = IntVector(d_dim, 0);
         shifted_pbox = pbox;
         if (loc_indx == 0) {
            shift(0) = -1;
         } else if (loc_indx == 1) {
            shift(0) = 1;
         } else if (loc_indx == 2) {
            shift(1) = -1;
         } else if (loc_indx == 3) {
            shift(1) = 1;
         } else if (loc_indx == 4) {
            shift(2) = -1;
         } else if (loc_indx == 5) {
            shift(2) = 1;
         }
         shifted_pbox.shift(shift);
         fill_box = fill_box * shifted_pbox;

         CellIterator ciend(CellGeometry::end(fill_box));
         for (CellIterator ci(CellGeometry::begin(fill_box));
              ci != ciend; ++ci) {
            CellIndex cell = *ci;
            for (d = 0; d < ucell->getDepth(); ++d) {
               (*ucell)(cell, d) = 1.0;
            }
         }
      }
   }

   /*
    * Set cell weights to 1.0 on EDGES of patch.
    */
   for (i = 0; i < num_edge_bdry_boxes; ++i) {
      Box fill_box = patch_geom->getBoundaryFillBox(edge_bdry[i],
            pbox,
            ghost_cells);

      int loc_indx = edge_bdry[i].getLocationIndex();
      shift = IntVector(d_dim, 0);
      shifted_pbox = pbox;

      if (d_dim == tbox::Dimension(3)) {
         /*
          * location index:
          * 3D:
          *    0,1,2,3 - XloYlo, XhiYLO, XloYhi, XhiYhi
          *    4,5,6,7 - XloZlo, XhiZlo, XloZhi, XhiZhi
          *    8,9,10,11 - YloZlo, YhiZlo, YloZhi, YhiZhi
          */
         if (loc_indx == 0) {
            shift(0) = -1;
            shift(1) = -1;
         } else if (loc_indx == 1) {
            shift(0) = 1;
            shift(1) = -1;
         } else if (loc_indx == 2) {
            shift(0) = -1;
            shift(1) = 1;
         } else if (loc_indx == 3) {
            shift(0) = 1;
            shift(1) = 1;
         } else if (loc_indx == 4) {
            shift(0) = -1;
            shift(2) = -1;
         } else if (loc_indx == 5) {
            shift(0) = 1;
            shift(2) = -1;
         } else if (loc_indx == 6) {
            shift(0) = -1;
            shift(2) = 1;
         } else if (loc_indx == 7) {
            shift(0) = 1;
            shift(2) = 1;
         } else if (loc_indx == 8) {
            shift(1) = -1;
            shift(2) = -1;
         } else if (loc_indx == 9) {
            shift(1) = 1;
            shift(2) = -1;
         } else if (loc_indx == 10) {
            shift(1) = -1;
            shift(2) = 1;
         } else if (loc_indx == 11) {
            shift(1) = 1;
            shift(2) = 1;
         }
      }
      if (d_dim == tbox::Dimension(2)) {

         /*
          * location index:
          * 2D:
          *    0,1 - Xlo, Xhi
          *    2,3 - Ylo, Yhi
          */
         if (loc_indx == 0) {
            shift(0) = -1;
         } else if (loc_indx == 1) {
            shift(0) = 1;
         } else if (loc_indx == 2) {
            shift(1) = -1;
         } else if (loc_indx == 3) {
            shift(1) = 1;
         }
      }
      shifted_pbox.shift(shift);
      fill_box = fill_box * shifted_pbox;

      CellIterator ciend(CellGeometry::end(fill_box));
      for (CellIterator ci(CellGeometry::begin(fill_box));
           ci != ciend; ++ci) {
         CellIndex cell = *ci;
         for (d = 0; d < ucell->getDepth(); ++d) {
            (*ucell)(cell, d) = 1.0;
         }
      }
   }

   /*
    * Set cell weights to 1.0 on NODES of patch.
    */

   for (i = 0; i < num_node_bdry_boxes; ++i) {
      Box fill_box = patch_geom->getBoundaryFillBox(node_bdry[i],
            pbox,
            ghost_cells);
      CellIterator ciend(CellGeometry::end(fill_box));
      for (CellIterator ci(CellGeometry::begin(fill_box));
           ci != ciend; ++ci) {
         CellIndex cell = *ci;  //i,j

         for (d = 0; d < ucell->getDepth(); ++d) {
            (*ucell)(cell, d) = 1.0;
         }
      }
   }

}

/*
 * Zero out the cells on the physical boundary that lie at the
 * coarse-fine boundary.  This is needed for cases in which the
 * fine patch intersects the physical boundary.
 */
void HierSumTest::zeroOutPhysicalBoundaryCellsAtCoarseFineBoundary(
   Patch& cpatch,
   const int cell_data_id)
{

   const std::shared_ptr<CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<CartesianPatchGeometry, PatchGeometry>(
         cpatch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);

   /*
    * Get node and edge boundary boxes.
    */

   const std::vector<BoundaryBox>& edge_bdry =
      patch_geom->getCodimensionBoundaries(d_dim.getValue() - 1);
   const int num_edge_bdry_boxes = static_cast<int>(edge_bdry.size());

   std::vector<BoundaryBox> empty_vector(0, BoundaryBox(d_dim));
   const std::vector<BoundaryBox>& face_bdry =
      d_dim == tbox::Dimension(3) ?
      patch_geom->getCodimensionBoundaries(1) : empty_vector;
   const int num_face_bdry_boxes = static_cast<int>(face_bdry.size());

   /*
    * std::shared_ptr to data in ghost regions.
    */
   std::shared_ptr<CellData<double> > ucell(
      SAMRAI_SHARED_PTR_CAST<CellData<double>, PatchData>(
         cpatch.getPatchData(cell_data_id)));

   TBOX_ASSERT(ucell);

   IntVector ghost_cells(ucell->getGhostCellWidth());

   int d;
   Box cpbox = cpatch.getBox();
   double interior_value;

   if (d_dim == tbox::Dimension(3)) {
      /*
       * Zero out values on the FACE boundary when they border a value of
       * zero on the interior.
       */
      for (int i_box = 0; i_box < num_face_bdry_boxes; ++i_box) {
         Box fill_box = patch_geom->getBoundaryFillBox(face_bdry[i_box],
               cpbox,
               ghost_cells);
         /*
          * location index:
          *    0,1 - X lower,upper
          *    2,3 - Y lower,upper
          *    4,5 - Z lower,upper
          */
         int loc_indx = face_bdry[i_box].getLocationIndex();
         CellIterator ciend(CellGeometry::end(fill_box));
         for (CellIterator ci(CellGeometry::begin(fill_box));
              ci != ciend; ++ci) {
            CellIndex boundary_cell = *ci;
            CellIndex interior = boundary_cell;
            if (loc_indx == 0) {
               interior(0) += 1;
            } else if (loc_indx == 1) {
               interior(0) -= 1;
            } else if (loc_indx == 2) {
               interior(1) += 1;
            } else if (loc_indx == 3) {
               interior(1) -= 1;
            } else if (loc_indx == 4) {
               interior(2) += 1;
            } else if (loc_indx == 5) {
               interior(2) -= 1;
            }

            for (d = 0; d < ucell->getDepth(); ++d) {
               interior_value = (*ucell)(interior, d);

               if (tbox::MathUtilities<double>::equalEps(interior_value,
                      0.0)) {
                  (*ucell)(boundary_cell, d) = 0.0;
               }
            }
         }
      }
   }
   /*
    * Zero out values on the EDGE boundary when they border a value of
    * zero on the interior.
    */
   for (int i = 0; i < num_edge_bdry_boxes; ++i) {
      Box fill_box = patch_geom->getBoundaryFillBox(edge_bdry[i],
            cpbox,
            ghost_cells);

      /*
       * location index:
       *    0,1,2,3 - XloYlo, XhiYLo, XloYhi, XhiYhi
       *    4,5,6,7 - XloZlo, XhiZlo, XloZhi, XhiZhi
       *    8,9,10,11 - YloZlo, YhiZlo, YloZhi, YhiZhi
       * 2D:
       *    0,1 - Xlo, Xhi
       *    2,3 - Ylo, Yhi
       */
      int loc_indx = edge_bdry[i].getLocationIndex();
      CellIterator ciend(CellGeometry::end(fill_box));
      for (CellIterator ci(CellGeometry::begin(fill_box));
           ci != ciend; ++ci) {
         CellIndex boundary_cell = *ci;
         CellIndex interior = boundary_cell;
         if (d_dim == tbox::Dimension(3)) {
            if (loc_indx == 0) {
               interior(0) += 1;
               interior(1) += 1;
            } else if (loc_indx == 1) {
               interior(0) -= 1;
               interior(1) += 1;
            } else if (loc_indx == 2) {
               interior(0) += 1;
               interior(1) -= 1;
            } else if (loc_indx == 3) {
               interior(0) -= 1;
               interior(1) -= 1;
            } else if (loc_indx == 4) {
               interior(0) += 1;
               interior(2) += 1;
            } else if (loc_indx == 5) {
               interior(0) -= 1;
               interior(2) += 1;
            } else if (loc_indx == 6) {
               interior(0) += 1;
               interior(2) -= 1;
            } else if (loc_indx == 7) {
               interior(0) -= 1;
               interior(2) -= 1;
            } else if (loc_indx == 8) {
               interior(1) += 1;
               interior(2) += 1;
            } else if (loc_indx == 9) {
               interior(1) -= 1;
               interior(2) += 1;
            } else if (loc_indx == 10) {
               interior(1) += 1;
               interior(2) -= 1;
            } else if (loc_indx == 11) {
               interior(1) -= 1;
               interior(2) -= 1;
            }
         }
         if (d_dim == tbox::Dimension(2)) {
            if (loc_indx == 0) {
               interior(0) += 1;
            } else if (loc_indx == 1) {
               interior(0) -= 1;
            } else if (loc_indx == 2) {
               interior(1) += 1;
            } else if (loc_indx == 3) {
               interior(1) -= 1;
            }
         }
         for (d = 0; d < ucell->getDepth(); ++d) {
            interior_value = (*ucell)(interior, d);

            if (tbox::MathUtilities<double>::equalEps(interior_value, 0.0)) {
               (*ucell)(boundary_cell, d) = 0.0;
            }
         }
      }
   }

   /*
    * Since we never use cell values at corners (nodes), no need to set
    * them to anything.
    */
}

/*
 *************************************************************************
 *
 * Get data from input database.
 *
 *************************************************************************
 */
void
HierSumTest::getFromInput(
   std::shared_ptr<Database> input_db)
{
   /*
    * Set number of ghosts for node and edge data.
    */
   std::vector<int> tmp_array;
   if (input_db->keyExists("node_ghosts")) {
      tmp_array = input_db->getIntegerVector("node_ghosts");
      if (static_cast<int>(tmp_array.size()) != d_dim.getValue()) {
         TBOX_ERROR("HierSumTest::getFromInput()"
            << "invalid 'node_ghosts' entry - must be integer"
            << "array of size d_dim" << std::endl);
      }
   } else {
      tmp_array.resize(d_dim.getValue());
      for (int i = 0; i < d_dim.getValue(); ++i) {
         tmp_array[i] = 0;
      }
   }

   for (int i = 0; i < d_dim.getValue(); ++i) {
      d_node_ghosts(i) = tmp_array[i];
   }

   if (input_db->keyExists("edge_ghosts")) {
      tmp_array = input_db->getIntegerVector("edge_ghosts");
      if (static_cast<int>(tmp_array.size()) != d_dim.getValue()) {
         TBOX_ERROR("HierSumTest::getFromInput()"
            << "invalid 'edge_ghosts' entry - must be integer"
            << "array of size d_dim" << std::endl);
      }
   } else {
      tmp_array.resize(d_dim.getValue());
      for (int i = 0; i < d_dim.getValue(); ++i) {
         tmp_array[i] = 0;
      }
   }

   for (int i = 0; i < d_dim.getValue(); ++i) {
      d_edge_ghosts(i) = tmp_array[i];
   }

   if (input_db->keyExists("var_depth")) {
      d_depth = input_db->getInteger("var_depth");
   }

   /*
    * See if we want to check data before communication
    */
   d_check_data_before_communication = false;
   if (input_db->keyExists("check_data_before_communication")) {
      d_check_data_before_communication =
         input_db->getBool("check_data_before_communication");
   }

}

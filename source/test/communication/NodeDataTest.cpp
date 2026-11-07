/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   AMR communication tests for node-centered patch data
 *
 ************************************************************************/

#include "NodeDataTest.h"

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#include "SAMRAI/pdat/NodeIndex.h"
#include "SAMRAI/pdat/NodeIterator.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "CommTester.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"

#include <vector>

namespace SAMRAI {

NodeDataTest::NodeDataTest(
   const std::string& object_name,
   const tbox::Dimension& dim,
   std::shared_ptr<tbox::Database> main_input_db,
   bool do_refine,
   bool do_coarsen,
   const std::string& refine_option):
   PatchDataTestStrategy(dim),
   d_dim(dim)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(main_input_db);
   TBOX_ASSERT(!refine_option.empty());

   d_object_name = object_name;

   d_do_refine = do_refine;
   d_do_coarsen = false;
   if (!do_refine) {
      d_do_coarsen = do_coarsen;
   }

   d_refine_option = refine_option;

   d_Acoef = 0.0;
   d_Bcoef = 0.0;
   d_Ccoef = 0.0;
   d_Dcoef = 0.0;

   d_finest_level_number = main_input_db->
      getDatabase("PatchHierarchy")->
      getInteger("max_levels") - 1;

   d_cart_grid_geometry.reset(
      new geom::CartesianGridGeometry(
         dim,
         "CartesianGridGeometry",
         main_input_db->getDatabase("CartesianGridGeometry")));

   setGridGeometry(d_cart_grid_geometry);

   readTestInput(main_input_db->getDatabase("NodePatchDataTest"));

}

NodeDataTest::~NodeDataTest()
{
}

void NodeDataTest::readTestInput(
   std::shared_ptr<tbox::Database> db)
{
   TBOX_ASSERT(db);

   /*
    * Read coeeficients of linear profile to test interpolation.
    */
   if (db->keyExists("Acoef")) {
      d_Acoef = db->getDouble("Acoef");
   } else {
      TBOX_ERROR(d_object_name << " input error: No `Acoeff' found." << std::endl);
   }
   if (db->keyExists("Dcoef")) {
      d_Dcoef = db->getDouble("Dcoef");
   } else {
      TBOX_ERROR(d_object_name << " input error: No `Dcoef' found." << std::endl);
   }
   if (d_dim > tbox::Dimension(1)) {
      if (db->keyExists("Bcoef")) {
         d_Bcoef = db->getDouble("Bcoef");
      } else {
         TBOX_ERROR(d_object_name << " input error: No `Bcoef' found." << std::endl);
      }
   }
   if (d_dim > tbox::Dimension(2)) {
      if (db->keyExists("Ccoef")) {
         d_Ccoef = db->getDouble("Ccoef");
      } else {
         TBOX_ERROR(d_object_name << " input error: No `Ccoef' found." << std::endl);
      }
   }

   /*
    * Base class reads variable parameters and boxes to refine.
    */

   readVariableInput(db->getDatabase("VariableData"));
}

void NodeDataTest::registerVariables(
   CommTester* commtest)
{
   TBOX_ASSERT(commtest != 0);

   int nvars = static_cast<int>(d_variable_src_name.size());

   d_variables.resize(nvars);

   for (int i = 0; i < nvars; ++i) {
      d_variables[i].reset(
         new pdat::NodeVariable<NODE_KERNEL_TYPE>(d_dim,
            d_variable_src_name[i],
            d_variable_depth[i]));

      if (d_do_refine) {
         commtest->registerVariable(d_variables[i],
            d_variables[i],
            d_variable_src_ghosts[i],
            d_variable_dst_ghosts[i],
            d_cart_grid_geometry,
            d_variable_refine_op[i]);
      } else if (d_do_coarsen) {
         commtest->registerVariable(d_variables[i],
            d_variables[i],
            d_variable_src_ghosts[i],
            d_variable_dst_ghosts[i],
            d_cart_grid_geometry,
            d_variable_coarsen_op[i]);
      }

   }

}

void NodeDataTest::setLinearData(
   std::shared_ptr<pdat::NodeData<NODE_KERNEL_TYPE> > data,
   const hier::Box& box,
   const hier::Patch& patch) const
{
   TBOX_ASSERT(data);

   std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(pgeom);
   const pdat::NodeIndex loweri(
      patch.getBox().lower(), (pdat::NodeIndex::Corner)0);
   const double* dx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   NODE_KERNEL_TYPE x, y, z;

   const int depth = data->getDepth();

   const hier::Box sbox = data->getGhostBox() * box;

   pdat::NodeIterator ciend(pdat::NodeGeometry::end(sbox));
   for (pdat::NodeIterator ci(pdat::NodeGeometry::begin(sbox));
        ci != ciend; ++ci) {

      /*
       * Compute spatial location of node center and
       * set data to linear profile.
       */

      x = lowerx[0] + dx[0] * ((*ci)(0) - loweri(0));
      y = z = 0.;
      if (d_dim > tbox::Dimension(1)) {
         y = lowerx[1] + dx[1] * ((*ci)(1) - loweri(1));
      }
      if (d_dim > tbox::Dimension(2)) {
         z = lowerx[2] + dx[2] * ((*ci)(2) - loweri(2));
      }

      for (int d = 0; d < depth; ++d) {
         (*data)(*ci, d) = static_cast<NODE_KERNEL_TYPE>(d_Dcoef + d_Acoef * x + d_Bcoef * y + d_Ccoef * z);
      }

   }

}

void NodeDataTest::setPeriodicData(
   std::shared_ptr<pdat::NodeData<NODE_KERNEL_TYPE> > data,
   const hier::Box& box,
   const hier::Patch& patch) const
{
   TBOX_ASSERT(data);

   NULL_USE(patch);

   const double* xlo = d_cart_grid_geometry->getXLower();
   const double* xup = d_cart_grid_geometry->getXUpper();
   std::vector<double> domain_len(d_dim.getValue());
   for (int d = 0; d < d_dim.getValue(); ++d) {
      domain_len[d] = xup[d] - xlo[d];
   }

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);
   const double* dx = patch_geom->getDx();

   const int depth = data->getDepth();

   const hier::Box sbox = data->getGhostBox() * box;

   pdat::NodeIterator niend(pdat::NodeGeometry::end(sbox));
   for (pdat::NodeIterator ni(pdat::NodeGeometry::begin(sbox));
        ni != niend; ++ni) {

      NODE_KERNEL_TYPE val = 1.0;
      for (int d = 0; d < d_dim.getValue(); ++d) {
         double tmpf = dx[d] * (*ni)(d) / domain_len[d];
         tmpf = sin(2 * M_PI * tmpf);
         val *= static_cast<NODE_KERNEL_TYPE>(tmpf);
      }
      val = val + 20.0; // Shift function range to [1,3] to avoid bad floating point compares.
      for (int d = 0; d < depth; ++d) {
         (*data)(*ni, d) = val;
      }

   }

}

void NodeDataTest::initializeDataOnPatch(
   const hier::Patch& patch,
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number,
   char src_or_dst)
{
   NULL_USE(src_or_dst);
   NULL_USE(level_number);
   NULL_USE(hierarchy);

   const hier::IntVector periodic_shift(
      d_cart_grid_geometry->getPeriodicShift(hier::IntVector(d_dim, 1)));
   bool is_periodic = periodic_shift.max() > 0;

   if (d_do_refine) {

      for (int i = 0; i < static_cast<int>(d_variables.size()); ++i) {

         std::shared_ptr<pdat::NodeData<NODE_KERNEL_TYPE> > node_data(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<NODE_KERNEL_TYPE>, hier::PatchData>(
               patch.getPatchData(d_variables[i], getDataContext())));
         TBOX_ASSERT(node_data);

         hier::Box dbox = node_data->getBox();

         if (is_periodic) {
            setPeriodicData(node_data, dbox, patch);
         } else {
            setLinearData(node_data, dbox, patch);
         }

      }

   } else if (d_do_coarsen) {

      for (int i = 0; i < static_cast<int>(d_variables.size()); ++i) {

         std::shared_ptr<pdat::NodeData<NODE_KERNEL_TYPE> > node_data(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<NODE_KERNEL_TYPE>, hier::PatchData>(
               patch.getPatchData(d_variables[i], getDataContext())));
         TBOX_ASSERT(node_data);

         hier::Box dbox = node_data->getGhostBox();

         if (is_periodic) {
            setPeriodicData(node_data, dbox, patch);
         } else {
            setLinearData(node_data, dbox, patch);
         }

      }

   }

}

void NodeDataTest::checkPatchInteriorData(
   const std::shared_ptr<pdat::NodeData<NODE_KERNEL_TYPE> >& data,
   const hier::Box& interior,
   const hier::Patch& patch) const
{
   TBOX_ASSERT(data);

   const bool is_periodic =
      d_cart_grid_geometry->getPeriodicShift(hier::IntVector(d_dim,
            1)).max() > 0;

   const int depth = data->getDepth();

   std::shared_ptr<pdat::NodeData<NODE_KERNEL_TYPE> > correct_data(
      new pdat::NodeData<NODE_KERNEL_TYPE>(
         data->getBox(),
         depth,
         data->getGhostCellWidth()));
   if (is_periodic) {
      setPeriodicData(correct_data, correct_data->getGhostBox(), patch);
   } else {
      setLinearData(correct_data, correct_data->getGhostBox(), patch);
   }

   pdat::NodeIterator niend(pdat::NodeGeometry::end(interior));
   for (pdat::NodeIterator ni(pdat::NodeGeometry::begin(interior));
        ni != niend; ++ni) {
      for (int d = 0; d < depth; ++d) {
         if (!(tbox::MathUtilities<NODE_KERNEL_TYPE>::equalEps((*data)(*ni, d),
                  (*correct_data)(*ni, d)))) {
            tbox::perr << "FAILED: -- patch interior not properly filled"
                       << std::endl;
         }
      }
   }

}

void NodeDataTest::setPhysicalBoundaryConditions(
   const hier::Patch& patch,
   const double time,
   const hier::IntVector& gcw) const
{
   NULL_USE(time);

   const hier::IntVector periodic_shift =
      d_cart_grid_geometry->getPeriodicShift(hier::IntVector(d_dim, 1));
   bool is_periodic = periodic_shift.max() > 0;

   std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(pgeom);

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

      std::shared_ptr<pdat::NodeData<NODE_KERNEL_TYPE> > node_data(
         SAMRAI_SHARED_PTR_CAST<pdat::NodeData<NODE_KERNEL_TYPE>, hier::PatchData>(
            patch.getPatchData(d_variables[i], getDataContext())));
      TBOX_ASSERT(node_data);

      hier::Box patch_interior = node_data->getBox();
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      checkPatchInteriorData(node_data, patch_interior, patch);

      /*
       * Set node boundary data.
       */
      for (int ni = 0; ni < num_node_bdry_boxes; ++ni) {

         hier::Box fill_box = pgeom->getBoundaryFillBox(node_bdry[ni],
               patch.getBox(),
               gcw);

         if (is_periodic) {
            setPeriodicData(node_data, fill_box, patch);
         } else {
            setLinearData(node_data, fill_box, patch);
         }
      }

      if (d_dim > tbox::Dimension(1)) {
         /*
          * Set edge boundary data.
          */
         for (int ei = 0; ei < num_edge_bdry_boxes; ++ei) {

            hier::Box fill_box = pgeom->getBoundaryFillBox(edge_bdry[ei],
                  patch.getBox(),
                  gcw);

            if (is_periodic) {
               setPeriodicData(node_data, fill_box, patch);
            } else {
               setLinearData(node_data, fill_box, patch);
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
                  gcw);

            if (is_periodic) {
               setPeriodicData(node_data, fill_box, patch);
            } else {
               setLinearData(node_data, fill_box, patch);
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
bool NodeDataTest::verifyResults(
   const hier::Patch& patch,
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number)
{
   NULL_USE(hierarchy);

   bool test_failed = false;

   const hier::IntVector periodic_shift(
      d_cart_grid_geometry->getPeriodicShift(hier::IntVector(d_dim, 1)));
   bool is_periodic = periodic_shift.max() > 0;

   if (d_do_refine || d_do_coarsen) {

      tbox::plog << "\nEntering NodeDataTest::verifyResults..." << std::endl;
      tbox::plog << "level_number = " << level_number << std::endl;
      tbox::plog << "Patch box = " << patch.getBox() << std::endl;

      hier::IntVector tgcw(periodic_shift.getDim(), 0);
      for (int i = 0; i < static_cast<int>(d_variables.size()); ++i) {
         tgcw.max(patch.getPatchData(d_variables[i], getDataContext())->
            getGhostCellWidth());
      }
      hier::Box pbox = patch.getBox();

      std::shared_ptr<pdat::NodeData<NODE_KERNEL_TYPE> > solution(
         new pdat::NodeData<NODE_KERNEL_TYPE>(pbox, 1, tgcw));

      hier::Box gbox(pbox);
      gbox.grow(tgcw);

      if (d_do_refine) {
         if (is_periodic) {
            setPeriodicData(solution, gbox, patch);
         } else {
            setLinearData(solution, gbox, patch);
         }
      } else {
         if (is_periodic) {
            setPeriodicData(solution, gbox, patch);
         } else {
            setLinearData(solution, gbox, patch);
         }
      }

      for (int i = 0; i < static_cast<int>(d_variables.size()); ++i) {

         std::shared_ptr<pdat::NodeData<NODE_KERNEL_TYPE> > node_data(
            SAMRAI_SHARED_PTR_CAST<pdat::NodeData<NODE_KERNEL_TYPE>, hier::PatchData>(
               patch.getPatchData(d_variables[i], getDataContext())));
         TBOX_ASSERT(node_data);
         int depth = node_data->getDepth();
         hier::Box dbox = node_data->getGhostBox();

         pdat::NodeIterator ciend(pdat::NodeGeometry::end(dbox));
         for (pdat::NodeIterator ci(pdat::NodeGeometry::begin(dbox));
              ci != ciend; ++ci) {
            NODE_KERNEL_TYPE correct = (*solution)(*ci);
            for (int d = 0; d < depth; ++d) {
               NODE_KERNEL_TYPE result = (*node_data)(*ci, d);
               if (!tbox::MathUtilities<NODE_KERNEL_TYPE>::equalEps(correct, result)) {
                  tbox::perr << "Test FAILED: ...."
                             << " : node index = " << *ci
                             << " of L" << level_number
                             << " P" << patch.getLocalId()
                             << " " << patch.getBox() << std::endl;
                  tbox::perr << "    hier::Variable = "
                             << d_variable_src_name[i]
                             << " : depth index = " << d << std::endl;
                  tbox::perr << "    result = " << result
                             << " : correct = " << correct << std::endl;
                  test_failed = true;
               }
            }
         }

      }
      if (!test_failed) {
         tbox::plog << "Node test Successful!" << std::endl;
      }

      solution.reset();   // just to be anal...

      tbox::plog << "\nExiting NodeDataTest::verifyResults..." << std::endl;
      tbox::plog << "level_number = " << level_number << std::endl;
      tbox::plog << "Patch box = " << patch.getBox() << std::endl << std::endl;

   }

   return !test_failed;

}

#ifdef SAMRAI_HAVE_CONDUIT
void NodeDataTest::addFields(
   conduit::Node& node,
   int domain_id,
   const std::shared_ptr<hier::Patch>& patch)
{

   std::shared_ptr<hier::VariableContext> source =
      hier::VariableDatabase::getDatabase()->getContext("SOURCE");

   std::shared_ptr<pdat::NodeData<NODE_KERNEL_TYPE> > node_data(
      SAMRAI_SHARED_PTR_CAST<pdat::NodeData<NODE_KERNEL_TYPE>, hier::PatchData>(
         patch->getPatchData(d_variables[0], source)));

   size_t data_size = node_data->getGhostBox().size();

   std::string mesh_name =
      "domain_" + tbox::Utilities::intToString(domain_id, 6);

   for (int d = 0; d < node_data->getDepth(); ++d) {
      std::string data_name = "node_data_" + tbox::Utilities::intToString(d);
      node_data->putBlueprintField(node[mesh_name], data_name, "mesh", d);
   }

}
#endif


}

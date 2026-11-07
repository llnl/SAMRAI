/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   AMR communication tests for cell-centered patch data
 *
 ************************************************************************/

#include "CellDataTest.h"

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "CommTester.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/Database.h"

#include <vector>
#include <iostream>


namespace SAMRAI {

CellDataTest::CellDataTest(
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
      new geom::CartesianGridGeometry(dim,
         "CartesianGridGeometry",
         main_input_db->getDatabase("CartesianGridGeometry")));

   setGridGeometry(d_cart_grid_geometry);

   readTestInput(main_input_db->getDatabase("CellPatchDataTest"));

}

CellDataTest::~CellDataTest()
{
}

void CellDataTest::readTestInput(
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

void CellDataTest::registerVariables(
   CommTester* commtest)
{
   TBOX_ASSERT(commtest != 0);

   int nvars = static_cast<int>(d_variable_src_name.size());

   d_variables.resize(nvars);

   for (int i = 0; i < nvars; ++i) {
      d_variables[i].reset(
         new pdat::CellVariable<CELL_KERNEL_TYPE>(d_dim,
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

void CellDataTest::setLinearData(
   std::shared_ptr<pdat::CellData<CELL_KERNEL_TYPE> > data,
   const hier::Box& box,
   const hier::Patch& patch) const
{
   TBOX_ASSERT(data);

   std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(pgeom);
   const pdat::CellIndex loweri(patch.getBox().lower());
   const pdat::CellIndex upperi(patch.getBox().upper());
   const double* pdx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x, y, z;

   const int depth = data->getDepth();

   const hier::Box sbox = data->getGhostBox() * box;

   pdat::CellIterator ciend(pdat::CellGeometry::end(sbox));
   for (pdat::CellIterator ci(pdat::CellGeometry::begin(sbox));
        ci != ciend; ++ci) {

      /*
       * Compute spatial location of cell center and
       * set data to linear profile.
       */

      x = lowerx[0] + pdx[0] * ((*ci)(0) - loweri(0) + 0.5);
      y = z = 0.;
      if (d_dim > tbox::Dimension(1)) {
         y = lowerx[1] + pdx[1] * ((*ci)(1) - loweri(1) + 0.5);
      }
      if (d_dim > tbox::Dimension(2)) {
         z = lowerx[2] + pdx[2] * ((*ci)(2) - loweri(2) + 0.5);
      }

      for (int d = 0; d < depth; ++d) {
         (*data)(*ci, d) = d_Dcoef + d_Acoef * x + d_Bcoef * y + d_Ccoef * z;
      }

   }

}

void CellDataTest::setConservativeData(
   std::shared_ptr<pdat::CellData<CELL_KERNEL_TYPE> > data,
   const hier::Box& box,
   const hier::Patch& patch,
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= hierarchy->getFinestLevelNumber()));

   std::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(level_number));

   const hier::BoxContainer& domain =
      level->getPhysicalDomain(hier::BlockId::zero());
   size_t ncells = 0;
   for (hier::BoxContainer::const_iterator i = domain.begin();
        i != domain.end(); ++i) {
      ncells += i->size();
   }

   const int depth = data->getDepth();

   const hier::Box sbox = data->getGhostBox() * box;

   if (level_number == 0) {

      /*
       * Set cell value on level zero to u(i,j,k) = (i + j + k)/ncells.
       */

      pdat::CellIterator fiend(pdat::CellGeometry::end(sbox));
      for (pdat::CellIterator fi(pdat::CellGeometry::begin(sbox));
           fi != fiend; ++fi) {
         double value = 0.0;
         for (int d = 0; d < d_dim.getValue(); ++d) {
            value += (double)((*fi)(d));
         }
         value /= static_cast<double>(ncells);
         for (int dep = 0; dep < depth; ++dep) {
            (*data)(*fi, dep) = value;
         }
      }

   } else {

      /*
       * Set cell value on level > 0 to
       *    u(i,j,k) = u_c + ci*del_i + cj*del_j + ck*del_k
       * where u_c is underlying coarse value, (ci,cj,ck) is
       * the underlying coarse cell index, and (del_i,del_j,del_k)
       * is the vector between the coarse and fine cell centers.
       */

      const hier::IntVector& ratio = level->getRatioToLevelZero();

      std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
         SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
      TBOX_ASSERT(pgeom);
      const double* dx = pgeom->getDx();

      size_t coarse_ncells = ncells;
      std::vector<std::vector<double> > delta(d_dim.getValue());
      for (int d = 0; d < d_dim.getValue(); ++d) {
         delta[d].resize(ratio(d), 0.0);
         coarse_ncells /= ratio(d);
         double coarse_dx = dx[d] * ratio(d);
         for (int i = 0; i < ratio(d); ++i) {
            /*
             * delta[d][i] is the physical distance from i-th fine
             * cell centroid in d-direction to coarse cell centroid.
             * The distance is the d-th component of the displacement
             * vector.
             */
            delta[d][i] = (i + 0.5) * dx[d] - coarse_dx * 0.5;
         }
      }

      pdat::CellIterator fiend(pdat::CellGeometry::end(sbox));
      for (pdat::CellIterator fi(pdat::CellGeometry::begin(sbox));
           fi != fiend; ++fi) {

         const hier::IntVector ci(hier::Index::coarsen(*fi, ratio));
         hier::IntVector del(ci.getDim());  // Index vector from ci to fi.
         double value = 0.0;
         for (int d = 0; d < d_dim.getValue(); ++d) {
            del(d) = (int)delta[d][(*fi)(d) - ci(d) * ratio(d)];
            value += (double)(ci(d));
         }
         value /= static_cast<double>(coarse_ncells);

         for (int d = 0; d < d_dim.getValue(); ++d) {
            value += ci(d) * del(d);
         }

         for (int dep = 0; dep < depth; ++dep) {
            (*data)(*fi, dep) = value;
         }

      }

   }

}

void CellDataTest::setPeriodicData(
   std::shared_ptr<pdat::CellData<CELL_KERNEL_TYPE> > data,
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

   pdat::CellIterator ciend(pdat::CellGeometry::end(sbox));
   for (pdat::CellIterator ci(pdat::CellGeometry::begin(sbox));
        ci != ciend; ++ci) {

      double val = 1.0;
      for (int d = 0; d < d_dim.getValue(); ++d) {
         double tmpf = dx[d] * ((*ci)(d) + 0.5) / domain_len[d];
         tmpf = sin(2 * M_PI * tmpf);
         val *= tmpf;
      }
      for (int d = 0; d < depth; ++d) {
         (*data)(*ci, d) = val;
      }

   }

}

void CellDataTest::initializeDataOnPatch(
   const hier::Patch& patch,
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number,
   char src_or_dst)
{
   NULL_USE(src_or_dst);

   const hier::IntVector periodic_shift(d_cart_grid_geometry->getPeriodicShift(
                                           hier::IntVector(d_dim, 1)));
   bool is_periodic = periodic_shift.max() > 0;

   if (d_do_refine) {

      for (int i = 0; i < static_cast<int>(d_variables.size()); ++i) {

         std::shared_ptr<pdat::CellData<CELL_KERNEL_TYPE> > cell_data(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<CELL_KERNEL_TYPE>, hier::PatchData>(
               patch.getPatchData(d_variables[i], getDataContext())));
         TBOX_ASSERT(cell_data);

         hier::Box dbox = cell_data->getBox();

         if (is_periodic) {
            setPeriodicData(cell_data, dbox, patch);
         } else {
            setLinearData(cell_data, dbox, patch);
         }

      }

   } else if (d_do_coarsen) {

      for (int i = 0; i < static_cast<int>(d_variables.size()); ++i) {

         std::shared_ptr<pdat::CellData<CELL_KERNEL_TYPE> > cell_data(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<CELL_KERNEL_TYPE>, hier::PatchData>(
               patch.getPatchData(d_variables[i], getDataContext())));
         TBOX_ASSERT(cell_data);

         hier::Box dbox = cell_data->getGhostBox();

         if (is_periodic) {
            setPeriodicData(cell_data, dbox, patch);
         } else {
            setConservativeData(cell_data, dbox, patch, hierarchy, level_number);
         }

      }

   }

}

void CellDataTest::checkPatchInteriorData(
   const std::shared_ptr<pdat::CellData<CELL_KERNEL_TYPE> >& data,
   const hier::Box& interior,
   const hier::Patch& patch) const
{
   TBOX_ASSERT(data);

   const bool is_periodic =
      d_cart_grid_geometry->getPeriodicShift(hier::IntVector(d_dim,
            1)).max() > 0;

   const int depth = data->getDepth();

   std::shared_ptr<pdat::CellData<CELL_KERNEL_TYPE> > correct_data(
      new pdat::CellData<CELL_KERNEL_TYPE>(
         data->getBox(),
         depth,
         data->getGhostCellWidth()));
   if (is_periodic) {
      setPeriodicData(correct_data, correct_data->getGhostBox(), patch);
   } else {
      setLinearData(correct_data, correct_data->getGhostBox(), patch);
   }

   pdat::CellIterator ciend(pdat::CellGeometry::end(interior));
   for (pdat::CellIterator ci(pdat::CellGeometry::begin(interior));
        ci != ciend; ++ci) {
      for (int d = 0; d < depth; ++d) {
         if (!(tbox::MathUtilities<CELL_KERNEL_TYPE>::equalEps((*data)(*ci, d),
                  (*correct_data)(*ci, d)))) {
            tbox::perr << "FAILED: -- patch interior not properly filled: "
                       << (*data)(*ci,d)
                       << " vs "
                       << (*correct_data)(*ci,d)
                       << std::endl;
         }
      }
   }

}

void CellDataTest::setPhysicalBoundaryConditions(
   const hier::Patch& patch,
   const double time,
   const hier::IntVector& gcw_to_fill) const
{
   NULL_USE(time);

   const hier::IntVector periodic_shift(
      d_cart_grid_geometry->getPeriodicShift(hier::IntVector(d_dim, 1)));
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

      std::shared_ptr<pdat::CellData<CELL_KERNEL_TYPE> > cell_data(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<CELL_KERNEL_TYPE>, hier::PatchData>(
            patch.getPatchData(d_variables[i], getDataContext())));
      TBOX_ASSERT(cell_data);

      hier::Box patch_interior = cell_data->getBox();

      checkPatchInteriorData(cell_data, patch_interior, patch);

      /*
       * Set node boundary data.
       */
      for (int ni = 0; ni < num_node_bdry_boxes; ++ni) {

         hier::Box fill_box = pgeom->getBoundaryFillBox(node_bdry[ni],
               patch.getBox(),
               gcw_to_fill);

         if (is_periodic) {
            setPeriodicData(cell_data, fill_box, patch);
         } else {
            setLinearData(cell_data, fill_box, patch);
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

            if (is_periodic) {
               setPeriodicData(cell_data, fill_box, patch);
            } else {
               setLinearData(cell_data, fill_box, patch);
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

            if (is_periodic) {
               setPeriodicData(cell_data, fill_box, patch);
            } else {
               setLinearData(cell_data, fill_box, patch);
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
bool CellDataTest::verifyResults(
   const hier::Patch& patch,
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number)
{

   bool test_failed = false;

   const hier::IntVector periodic_shift(
      d_cart_grid_geometry->getPeriodicShift(hier::IntVector(d_dim, 1)));
   bool is_periodic = periodic_shift.max() > 0;

   if (d_do_refine || d_do_coarsen) {

      tbox::plog << "\nEntering CellDataTest::verifyResults..." << std::endl;
      tbox::plog << "level_number = " << level_number << std::endl;
      tbox::plog << "Patch box = " << patch.getBox() << std::endl;

      hier::IntVector tgcw(periodic_shift.getDim(), 0);
      for (int i = 0; i < static_cast<int>(d_variables.size()); ++i) {
         tgcw.max(patch.getPatchData(d_variables[i], getDataContext())->
            getGhostCellWidth());
      }
      hier::Box pbox = patch.getBox();

      std::shared_ptr<pdat::CellData<CELL_KERNEL_TYPE> > solution(
         new pdat::CellData<CELL_KERNEL_TYPE>(pbox, 1, tgcw));

      hier::Box tbox(pbox);
      tbox.grow(tgcw);

      if (d_do_refine) {
         if (is_periodic) {
            setPeriodicData(solution, tbox, patch);
         } else {
            setLinearData(solution, tbox, patch);
         }
      } else {
         if (is_periodic) {
            setPeriodicData(solution, tbox, patch);
         } else {
            setConservativeData(solution, tbox, patch, hierarchy, level_number);
         }
      }

      for (int i = 0; i < static_cast<int>(d_variables.size()); ++i) {

         std::shared_ptr<pdat::CellData<CELL_KERNEL_TYPE> > cell_data(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<CELL_KERNEL_TYPE>, hier::PatchData>(
               patch.getPatchData(d_variables[i], getDataContext())));
         TBOX_ASSERT(cell_data);
         int depth = cell_data->getDepth();
         hier::Box dbox = cell_data->getGhostBox();

         pdat::CellIterator ciend(pdat::CellGeometry::end(dbox));
         for (pdat::CellIterator ci(pdat::CellGeometry::begin(dbox));
              ci != ciend; ++ci) {
            CELL_KERNEL_TYPE correct = (*solution)(*ci);
            for (int d = 0; d < depth; ++d) {
               CELL_KERNEL_TYPE result = (*cell_data)(*ci, d);
               if (!tbox::MathUtilities<CELL_KERNEL_TYPE>::equalEps(correct, result)) {
                  tbox::perr << "Test FAILED: ...."
                             << " : cell index = " << *ci
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
            if (test_failed) break;
         }
      }
      if (!test_failed) {
         tbox::plog << "CellDataTest Successful!" << std::endl;
      }

      solution.reset();   // just to be anal...

      tbox::plog << "\nExiting CellDataTest::verifyResults..." << std::endl;
      tbox::plog << "level_number = " << level_number << std::endl;
      tbox::plog << "Patch box = " << patch.getBox() << std::endl << std::endl;

   }

   return !test_failed;

}

void CellDataTest::setDataIds(std::list<int>& data_ids)
{
   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
   for (int i = 0; i < static_cast<int>(d_variables.size()); ++i) {
      int data_id = variable_db->mapVariableAndContextToIndex(
         d_variables[i],
         getDataContext());
      data_ids.push_back(data_id);
   }
}

bool CellDataTest::verifyCompositeBoundaryData(
   const hier::Patch& patch,
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   int data_id,
   int level_number,
   const std::vector<std::shared_ptr<hier::PatchData> >& bdry_data)
{
   bool test_failed = false;

   const hier::IntVector periodic_shift(
      d_cart_grid_geometry->getPeriodicShift(hier::IntVector(d_dim, 1)));
   bool is_periodic = periodic_shift.max() > 0;

   if (d_do_refine && !is_periodic) {
      for (std::vector<std::shared_ptr<hier::PatchData> >::const_iterator
           itr = bdry_data.begin(); itr != bdry_data.end(); ++itr) {
         std::shared_ptr<pdat::CellData<CELL_KERNEL_TYPE> > cell_bdry_data(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<CELL_KERNEL_TYPE>, hier::PatchData>(*itr));
         TBOX_ASSERT(cell_bdry_data);

         hier::Patch solution_patch(
            cell_bdry_data->getBox(), patch.getPatchDescriptor());

         d_cart_grid_geometry->setGeometryDataOnPatch(
            solution_patch,
            hierarchy->getPatchLevel(level_number+1)->getRatioToLevelZero(),
            hier::PatchGeometry::TwoDimBool(patch.getDim(), false));

         std::shared_ptr<pdat::CellData<CELL_KERNEL_TYPE> > solution(
            new pdat::CellData<CELL_KERNEL_TYPE>(
               solution_patch.getBox(), 1,
               hier::IntVector::getZero(patch.getDim())));

         setLinearData(solution, solution_patch.getBox(), solution_patch);

         int depth = cell_bdry_data->getDepth();
         hier::Box dbox(cell_bdry_data->getBox());

         pdat::CellIterator ciend(pdat::CellGeometry::end(dbox));
         for (pdat::CellIterator ci(pdat::CellGeometry::begin(dbox));
              ci != ciend; ++ci) {
            CELL_KERNEL_TYPE correct = (*solution)(*ci);
            for (int d = 0; d < depth; ++d) {
               CELL_KERNEL_TYPE result = (*cell_bdry_data)(*ci, d);
               if (!tbox::MathUtilities<CELL_KERNEL_TYPE>::equalEps(correct, result)) {
                  tbox::perr << "Test FAILED: ...."
                             << " : cell index = " << *ci
                             << " on composite boundary stencil"
                             << " " << solution->getBox() << std::endl;
                  tbox::perr << "    patch data id " << data_id
                             << " : depth index = " << d << std::endl;
                  tbox::perr << "    result = " << result
                             << " : correct = " << correct << std::endl;
                  test_failed = true;
               }
            }
         }
      }
   }

   return !test_failed;

}

#ifdef SAMRAI_HAVE_CONDUIT
void CellDataTest::addFields(conduit::Node& node, int domain_id, const std::shared_ptr<hier::Patch>& patch)
{

   std::shared_ptr<hier::VariableContext> source =
      hier::VariableDatabase::getDatabase()->getContext("SOURCE");

   std::shared_ptr<pdat::CellData<CELL_KERNEL_TYPE> > cell_data(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<CELL_KERNEL_TYPE>, hier::PatchData>(
         patch->getPatchData(d_variables[0], source)));

   size_t data_size = cell_data->getGhostBox().size();

   std::string mesh_name =
      "domain_" + tbox::Utilities::intToString(domain_id, 6);

   for (int d = 0; d < cell_data->getDepth(); ++d) {
      std::string data_name = "cell_data_" + tbox::Utilities::intToString(d); 
      cell_data->putBlueprintField(node[mesh_name], data_name, "mesh", d);
   }
}
#endif

}

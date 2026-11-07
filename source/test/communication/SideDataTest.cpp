/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   AMR communication tests for side-centered patch data
 *
 ************************************************************************/

#include "SideDataTest.h"

#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "CommTester.h"
#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideIndex.h"
#include "SAMRAI/pdat/SideIterator.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/VariableDatabase.h"

#include <vector>

namespace SAMRAI {


SideDataTest::SideDataTest(
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

   readTestInput(main_input_db->getDatabase("SidePatchDataTest"));

}

SideDataTest::~SideDataTest()
{
}

void SideDataTest::readTestInput(
   std::shared_ptr<tbox::Database> db)
{
   TBOX_ASSERT(db);

   /*
    * Base class reads variable parameters and boxes to refine.
    */

   readVariableInput(db->getDatabase("VariableData"));

   std::shared_ptr<tbox::Database> var_data(
      db->getDatabase("VariableData"));
   std::vector<std::string> var_keys = var_data->getAllKeys();
   int nkeys = static_cast<int>(var_keys.size());

   d_test_direction.resize(nkeys, hier::IntVector::getZero(d_dim));
   d_use_fine_value_at_interface.resize(nkeys);

   for (int i = 0; i < nkeys; ++i) {
      std::shared_ptr<tbox::Database> var_db(
         var_data->getDatabase(var_keys[i]));

      int test_direction = var_db->getIntegerWithDefault("test_direction", -1);
      for (int j = 0; j < d_dim.getValue(); ++j) {
         if (test_direction == -1 || test_direction == j) {
            d_test_direction[i][j] = 1;
         }
      }

      if (var_db->keyExists("use_fine_value_at_interface")) {
         d_use_fine_value_at_interface[i] =
            var_db->getBool("use_fine_value_at_interface");
      } else {
         d_use_fine_value_at_interface[i] = true;
      }

   }

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

}

void SideDataTest::registerVariables(
   CommTester* commtest)
{
   TBOX_ASSERT(commtest != 0);

   int nvars = static_cast<int>(d_variable_src_name.size());

   d_variables.resize(nvars);

   for (int i = 0; i < nvars; ++i) {
      d_variables[i].reset(
         new pdat::SideVariable<double>(
            d_dim,
            d_variable_src_name[i],
            d_test_direction[i],
            d_variable_depth[i],
            d_use_fine_value_at_interface[i]));

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

void SideDataTest::setConservativeData(
   std::shared_ptr<pdat::SideData<double> > data,
   const hier::Box& box,
   const hier::Patch& patch,
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= hierarchy->getFinestLevelNumber()));

   int i, j;
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

   const hier::IntVector& directions(data->getDirectionVector());

   if (level_number == 0) {

      /*
       * Set side value on level zero as follows:
       *
       *    u0(i,j,k) = (j + k)/ncells
       *    u1(i,j,k) = (i + k)/ncells
       *    u2(i,j,k) = (i + j)/ncells
       */

      for (int axis = 0; axis < d_dim.getValue(); ++axis) {
         if (directions(axis)) {
            pdat::CellIterator ciend(pdat::CellGeometry::end(sbox));
            for (pdat::CellIterator ci(pdat::CellGeometry::begin(sbox));
                 ci != ciend; ++ci) {
               double value = 0.0;
               for (i = 0; i < d_dim.getValue(); ++i) {
                  if (i != axis) {
                     value += (double)((*ci)(i));
                  }
               }
               value /= static_cast<double>(ncells);
               for (int side = pdat::SideIndex::Lower;
                    side <= pdat::SideIndex::Upper; ++side) {
                  pdat::SideIndex si(*ci, axis, side);
                  for (int d = 0; d < depth; ++d) {
                     (*data)(si, d) = value;
                  }
               }
            }
         }
      }

   } else {

      /*
       * Set side value on level > 0 to
       *    u(i,j,k) = u_c + ci*del_i + cj*del_j + ck*del_k
       * where u_c is value on the underlying coarse side, (ci,cj,ck) is
       * the underlying coarse side index, and (del_i,del_j,del_k)
       * is the vector between the coarse and fine cell side centers.
       */

      hier::IntVector ratio(level->getRatioToLevelZero());
      const int max_ratio = ratio.max();

      std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
         SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
      TBOX_ASSERT(pgeom);
      const double* dx = pgeom->getDx();

      size_t coarse_ncells = ncells;
      double* delta = new double[max_ratio * d_dim.getValue()];
      for (j = 0; j < d_dim.getValue(); ++j) {
         coarse_ncells /= ratio(j);
         double coarse_dx = dx[j] * ratio(j);
         for (i = 0; i < ratio(j); ++i) {
            delta[j * max_ratio + i] = (i + 0.5) * dx[j] - coarse_dx * 0.5;
         }
      }

      for (tbox::Dimension::dir_t axis = 0; axis < d_dim.getValue(); ++axis) {
         if (directions(axis)) {
            hier::IntVector ci(ratio.getDim());
            hier::IntVector del(ratio.getDim());
            pdat::CellIterator fiend(pdat::CellGeometry::end(sbox));
            for (pdat::CellIterator fi(pdat::CellGeometry::begin(sbox));
                 fi != fiend; ++fi) {
               double value = 0.0;
               for (i = 0; i < d_dim.getValue(); ++i) {
                  if (i != axis) {
                     int findx = (*fi)(i);
                     ci(i) = ((findx < 0) ? (findx + 1) / ratio(i) - 1
                              : findx / ratio(i));
                     del(i) =
                        (int)delta[i * max_ratio + findx - ci(i) * ratio(i)];
                     value += (double)(ci(i));
                  }
               }
               value /= static_cast<double>(coarse_ncells);

               for (j = 0; j < d_dim.getValue(); ++j) {
                  if (j != axis) {
                     value += ci(j) * del(j);
                  }
               }

               for (int side = pdat::SideIndex::Lower;
                    side <= pdat::SideIndex::Upper; ++side) {
                  pdat::SideIndex si(*fi, axis, side);
                  for (int d = 0; d < depth; ++d) {
                     (*data)(si, d) = value;
                  }
               }
            }
         }
      }
      delete[] delta;

   }

}

void SideDataTest::initializeDataOnPatch(
   const hier::Patch& patch,
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number,
   char src_or_dst)
{
   NULL_USE(src_or_dst);
   NULL_USE(level_number);
   NULL_USE(hierarchy);

   const hier::IntVector periodic_shift(
      d_cart_grid_geometry->getPeriodicShift(hier::IntVector(tbox::Dimension(
               d_dim), 1)));
   const bool is_periodic = periodic_shift.max() > 0;

   if (d_do_refine) {

      for (int i = 0; i < static_cast<int>(d_variables.size()); ++i) {

         std::shared_ptr<pdat::SideData<double> > side_data(
            SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
               patch.getPatchData(d_variables[i], getDataContext())));
         TBOX_ASSERT(side_data);

         hier::Box dbox = side_data->getBox();

         if (is_periodic) {
            setPeriodicData(side_data, dbox, patch);
         } else {
            setLinearData(side_data, dbox, patch);
         }
      }

   } else if (d_do_coarsen) {

      for (int i = 0; i < static_cast<int>(d_variables.size()); ++i) {

         std::shared_ptr<pdat::SideData<double> > side_data(
            SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
               patch.getPatchData(d_variables[i], getDataContext())));
         TBOX_ASSERT(side_data);

         hier::Box dbox = side_data->getGhostBox();

         setConservativeData(side_data, dbox,
            patch, hierarchy, level_number);

      }

   }

}

void SideDataTest::checkPatchInteriorData(
   const std::shared_ptr<pdat::SideData<double> >& data,
   const hier::Box& interior,
   const hier::Patch& patch) const
{
   TBOX_ASSERT(data);

   const bool is_periodic =
      d_cart_grid_geometry->getPeriodicShift(hier::IntVector(tbox::Dimension(
               d_dim), 1)).max() > 0;

   const int depth = data->getDepth();

   std::shared_ptr<pdat::SideData<double> > correct_data(
      new pdat::SideData<double>(
         data->getBox(),
         depth,
         data->getGhostCellWidth(),
         data->getDirectionVector()));

   if (d_do_refine) {
      if (is_periodic) {
         setPeriodicData(correct_data, correct_data->getGhostBox(), patch);
      } else {
         setLinearData(correct_data, correct_data->getGhostBox(), patch);
      }
   }

   const hier::IntVector& directions(data->getDirectionVector());

   for (tbox::Dimension::dir_t axis = 0; axis < d_dim.getValue(); ++axis) {
      if (directions(axis)) {
         const pdat::SideIndex loweri(interior.lower(), axis, 0);
         pdat::SideIterator siend(pdat::SideGeometry::end(interior, axis));
         for (pdat::SideIterator si(pdat::SideGeometry::begin(interior, axis));
              si != siend; ++si) {
            for (int d = 0; d < depth; ++d) {
               if (!(tbox::MathUtilities<double>::equalEps((*data)(*si, d),
                        (*correct_data)(*si, d)))) {
                  tbox::perr << "FAILED: -- patch interior not properly filled"
                             << " : side_data index = "
                             << si->getAxis() << '/' << *si
                             << std::endl;
               }
            }
         }
      }
   }
}

void SideDataTest::setPhysicalBoundaryConditions(
   const hier::Patch& patch,
   const double time,
   const hier::IntVector& gcw) const
{
   NULL_USE(time);

   const hier::IntVector periodic_shift(
      d_cart_grid_geometry->getPeriodicShift(hier::IntVector(gcw.getDim(), 1)));

   const bool is_periodic = periodic_shift.max() > 0;

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

      std::shared_ptr<pdat::SideData<double> > side_data(
         SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
            patch.getPatchData(d_variables[i], getDataContext())));
      TBOX_ASSERT(side_data);

      hier::Box patch_interior = side_data->getBox();
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      checkPatchInteriorData(side_data, patch_interior, patch);

      /*
       * Set node boundary data.
       */
      for (int ni = 0; ni < num_node_bdry_boxes; ++ni) {

         hier::Box fill_box = pgeom->getBoundaryFillBox(node_bdry[ni],
               patch.getBox(),
               gcw);

         if (is_periodic) {
            setPeriodicData(side_data, fill_box, patch);
         } else {
            setLinearData(side_data, fill_box, patch);
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
               setPeriodicData(side_data, fill_box, patch);
            } else {
               setLinearData(side_data, fill_box, patch);
            }
         }
      }

      if (d_dim == tbox::Dimension(3)) {
         /*
          * Set face boundary data.
          */
         for (int fi = 0; fi < num_face_bdry_boxes; ++fi) {
            hier::Box fbox(face_bdry[fi].getBox());

            hier::Box fill_box = pgeom->getBoundaryFillBox(face_bdry[fi],
                  patch.getBox(),
                  gcw);

            if (is_periodic) {
               setPeriodicData(side_data, fill_box, patch);
            } else {
               setLinearData(side_data, fill_box, patch);
            }
         }
      }

   }

}

void SideDataTest::setLinearData(
   std::shared_ptr<pdat::SideData<double> > data,
   const hier::Box& box,
   const hier::Patch& patch) const
{
   TBOX_ASSERT(data);

   std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(pgeom);
   const double* dx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x, y, z;

   const int depth = data->getDepth();

   const hier::Box sbox = data->getGhostBox() * box;

   hier::IntVector directions(data->getDirectionVector());

   for (tbox::Dimension::dir_t axis = 0; axis < d_dim.getValue(); ++axis) {
      if (directions(axis)) {
         const pdat::SideIndex loweri(patch.getBox().lower(), axis, 0);
         pdat::SideIterator eiend(pdat::SideGeometry::end(sbox, axis));
         for (pdat::SideIterator ei(pdat::SideGeometry::begin(sbox, axis));
              ei != eiend; ++ei) {

            /*
             * Compute spatial location of cell center and
             * set data to linear profile.
             */

            if (axis == 0) {
               x = lowerx[0] + dx[0] * ((*ei)(0) - loweri(0));
            } else {
               x = lowerx[0] + dx[0] * ((*ei)(0) - loweri(0) + 0.5);
            }
            y = z = 0.;
            if (d_dim > tbox::Dimension(1)) {
               if (axis == 1) {
                  y = lowerx[1] + dx[1] * ((*ei)(1) - loweri(1));
               } else {
                  y = lowerx[1] + dx[1] * ((*ei)(1) - loweri(1) + 0.5);
               }
            }
            if (d_dim > tbox::Dimension(2)) {
               if (axis == 2) {
                  z = lowerx[2] + dx[2] * ((*ei)(2) - loweri(2));
               } else {
                  z = lowerx[2] + dx[2] * ((*ei)(2) - loweri(2) + 0.5);
               }
            }

            for (int d = 0; d < depth; ++d) {
               (*data)(*ei,
                       d) = d_Dcoef + d_Acoef * x + d_Bcoef * y + d_Ccoef * z;
            }

         }
      }
   }
}

void SideDataTest::setPeriodicData(
   std::shared_ptr<pdat::SideData<double> > data,
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

   hier::IntVector directions(data->getDirectionVector());

   for (tbox::Dimension::dir_t axis = 0; axis < d_dim.getValue(); ++axis) {
      if (directions(axis)) {
         const pdat::SideIndex loweri(patch.getBox().lower(), axis, 0);
         pdat::SideIterator siend(pdat::SideGeometry::end(sbox, axis));
         for (pdat::SideIterator si(pdat::SideGeometry::begin(sbox, axis));
              si != siend; ++si) {

            double val = 1.0;
            for (int d = 0; d < d_dim.getValue(); ++d) {
               double tmpf = d == axis ? (*si)(d) : 0.5 + (*si)(d);
               tmpf = tmpf * dx[d] / domain_len[d];
               tmpf = sin(2 * M_PI * tmpf);
               val *= tmpf;
            }
            val = val + 2.0; // Shift function range to [1,3] to avoid bad floating point compares.
            for (int d = 0; d < depth; ++d) {
               (*data)(*si, d) = val;
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

bool SideDataTest::verifyResults(
   const hier::Patch& patch,
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number)
{
   NULL_USE(hierarchy);

   bool test_failed = false;

   const hier::IntVector periodic_shift(
      d_cart_grid_geometry->getPeriodicShift(hier::IntVector(tbox::Dimension(
               d_dim), 1)));
   bool is_periodic = periodic_shift.max() > 0;

   if (d_do_refine || d_do_coarsen) {

      tbox::plog << "\nEntering SideDataTest::verifyResults..." << std::endl;
      tbox::plog << "level_number = " << level_number << std::endl;
      tbox::plog << "Patch box = " << patch.getBox() << std::endl;

      hier::IntVector tgcw(periodic_shift.getDim(), 0);
      for (int i = 0; i < static_cast<int>(d_variables.size()); ++i) {
         tgcw.max(patch.getPatchData(d_variables[i], getDataContext())->
            getGhostCellWidth());
      }
      const hier::Box& pbox = patch.getBox();

      std::shared_ptr<pdat::SideData<double> > solution(
         new pdat::SideData<double>(pbox, 1, tgcw));

      hier::Box tbox(pbox);
      tbox.grow(tgcw);

      if (d_do_refine) {
         if (is_periodic) {
            setPeriodicData(solution, tbox, patch);
         } else {
            setLinearData(solution, tbox, patch);
         }
      }
      if (d_do_coarsen) {
         setConservativeData(solution, tbox,
            patch, hierarchy, level_number);
      }

      for (int i = 0; i < static_cast<int>(d_variables.size()); ++i) {

         std::shared_ptr<pdat::SideData<double> > side_data(
            SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
               patch.getPatchData(d_variables[i], getDataContext())));
         TBOX_ASSERT(side_data);
         int depth = side_data->getDepth();
         hier::Box dbox = side_data->getGhostBox();

         hier::IntVector directions(side_data->getDirectionVector());

         for (tbox::Dimension::dir_t id = 0; id < d_dim.getValue(); ++id) {
            if (directions(id)) {
               pdat::SideIterator siend(pdat::SideGeometry::end(dbox, id));
               for (pdat::SideIterator si(pdat::SideGeometry::begin(dbox, id));
                    si != siend; ++si) {
                  double correct = (*solution)(*si);
                  for (int d = 0; d < depth; ++d) {
                     double result = (*side_data)(*si, d);
                     if (!tbox::MathUtilities<double>::equalEps(correct,
                            result)) {
                        test_failed = true;
                        tbox::perr << "Test FAILED: ...."
                                   << " : side_data index = "
                                   << si->getAxis() << '/' << *si << std::endl;
                        tbox::perr << "    hier::Variable = "
                                   << d_variable_src_name[i]
                                   << " : depth index = " << d << std::endl;
                        tbox::perr << "    result = " << result
                                   << " : correct = " << correct << std::endl;
                     }
                  }
               }
            }
         }

      }

      solution.reset();   // just to be anal...

      tbox::plog << "\nExiting SideDataTest::verifyResults..." << std::endl;
      tbox::plog << "level_number = " << level_number << std::endl;
      tbox::plog << "Patch box = " << patch.getBox() << std::endl << std::endl;

   }
   return !test_failed;
}

}

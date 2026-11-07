/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   AMR communication tests for outerface-centered patch data
 *
 ************************************************************************/

#include "OuterfaceDataTest.h"

#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "CommTester.h"
#include "SAMRAI/pdat/FaceGeometry.h"
#include "SAMRAI/pdat/FaceIndex.h"
#include "SAMRAI/pdat/FaceIterator.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/pdat/OuterfaceGeometry.h"
#include "SAMRAI/pdat/OuterfaceVariable.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/VariableDatabase.h"

namespace SAMRAI {


OuterfaceDataTest::OuterfaceDataTest(
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

   d_use_fine_value_at_interface.resize(0);

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

   readTestInput(main_input_db->getDatabase("OuterfacePatchDataTest"));

}

OuterfaceDataTest::~OuterfaceDataTest()
{
}

void OuterfaceDataTest::readTestInput(
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

   d_use_fine_value_at_interface.resize(nkeys);

   for (int i = 0; i < nkeys; ++i) {
      std::shared_ptr<tbox::Database> var_db(
         var_data->getDatabase(var_keys[i]));

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

void OuterfaceDataTest::registerVariables(
   CommTester* commtest)
{
   TBOX_ASSERT(commtest != 0);

   int nvars = static_cast<int>(d_variable_src_name.size());

   d_variables_src.resize(nvars);
   d_variables_dst.resize(nvars);

   for (int i = 0; i < nvars; ++i) {
      d_variables_src[i].reset(
         new pdat::OuterfaceVariable<OUTERFACE_KERNEL_TYPE>(d_dim,
            d_variable_src_name[i],
            d_variable_depth[i]));

      if (i % 2 == 0) {
         d_variables_dst[i].reset(
            new pdat::FaceVariable<OUTERFACE_KERNEL_TYPE>(
               d_dim,
               d_variable_dst_name[i],
               d_variable_depth[i],
               d_use_fine_value_at_interface[i]));
      } else {
         d_variables_dst[i].reset(
            new pdat::OuterfaceVariable<OUTERFACE_KERNEL_TYPE>(
               d_dim,
               d_variable_dst_name[i],
               d_variable_depth[i]));
      }

      if (d_do_refine) {
         commtest->registerVariable(d_variables_src[i],
            d_variables_dst[i],
            d_variable_src_ghosts[i],
            d_variable_dst_ghosts[i],
            d_cart_grid_geometry,
            d_variable_refine_op[i]);
      } else if (d_do_coarsen) {
         commtest->registerVariable(d_variables_src[i],
            d_variables_dst[i],
            d_variable_src_ghosts[i],
            d_variable_dst_ghosts[i],
            d_cart_grid_geometry,
            d_variable_coarsen_op[i]);
      }

   }

}

void OuterfaceDataTest::initializeDataOnPatch(
   const hier::Patch& patch,
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number,
   char src_or_dst)
{
   NULL_USE(hierarchy);
   NULL_USE(level_number);
   hier::VariableDatabase* variable_db =
      hier::VariableDatabase::getDatabase();
   variable_db->printClassData();
   std::vector<std::shared_ptr<hier::Variable> >& variables(
      src_or_dst == 's' ? d_variables_src : d_variables_dst);

   if (d_do_refine) {

      for (int i = 0; i < static_cast<int>(variables.size()); ++i) {

         std::shared_ptr<hier::PatchData> data(
            patch.getPatchData(variables[i], getDataContext()));

         TBOX_ASSERT(data);

         std::shared_ptr<pdat::OuterfaceData<OUTERFACE_KERNEL_TYPE> > oface_data(
            std::dynamic_pointer_cast<pdat::OuterfaceData<OUTERFACE_KERNEL_TYPE>,
                                        hier::PatchData>(data));
         std::shared_ptr<pdat::FaceData<OUTERFACE_KERNEL_TYPE> > face_data(
            std::dynamic_pointer_cast<pdat::FaceData<OUTERFACE_KERNEL_TYPE>,
                                        hier::PatchData>(data));

         hier::Box dbox = data->getBox();

         if (face_data) {
            setLinearData(face_data, dbox, patch);
         }
         if (oface_data) {
            setLinearData(oface_data, dbox, patch);
         }

      }

   } else if (d_do_coarsen) {

      for (int i = 0; i < static_cast<int>(variables.size()); ++i) {

         std::shared_ptr<hier::PatchData> data(
            patch.getPatchData(variables[i], getDataContext()));

         TBOX_ASSERT(data);

         std::shared_ptr<pdat::OuterfaceData<OUTERFACE_KERNEL_TYPE> > oface_data(
            std::dynamic_pointer_cast<pdat::OuterfaceData<OUTERFACE_KERNEL_TYPE>,
                                        hier::PatchData>(data));
         std::shared_ptr<pdat::FaceData<OUTERFACE_KERNEL_TYPE> > face_data(
            std::dynamic_pointer_cast<pdat::FaceData<OUTERFACE_KERNEL_TYPE>,
                                        hier::PatchData>(data));

         hier::Box dbox = data->getGhostBox();

         if (face_data) {
            setLinearData(face_data, dbox, patch);
         }
         if (oface_data) {
            setLinearData(oface_data, dbox, patch);
         }

      }

   }

}

void OuterfaceDataTest::checkPatchInteriorData(
   const std::shared_ptr<pdat::OuterfaceData<OUTERFACE_KERNEL_TYPE> >& data,
   const hier::Box& interior,
   const std::shared_ptr<geom::CartesianPatchGeometry>& pgeom) const
{
   TBOX_ASSERT(data);

   const double* dx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x = 0., y = 0., z = 0.;

   const int depth = data->getDepth();

   for (tbox::Dimension::dir_t axis = 0; axis < d_dim.getValue(); ++axis) {
      const pdat::FaceIndex loweri(interior.lower(), axis, 0);
      pdat::FaceIterator fiend(pdat::FaceGeometry::end(interior, axis));
      for (pdat::FaceIterator fi(pdat::FaceGeometry::begin(interior, axis));
           fi != fiend; ++fi) {

         /*
          * Compute spatial location of face and
          * set data to linear profile.
          */

         if (axis == 0) {
            x = lowerx[0] + dx[0] * ((*fi)(0) - loweri(0));
            if (d_dim > tbox::Dimension(1)) {
               y = lowerx[1] + dx[1] * ((*fi)(1) - loweri(1) + 0.5);
            }
            if (d_dim > tbox::Dimension(2)) {
               z = lowerx[2] + dx[2] * ((*fi)(2) - loweri(2) + 0.5);
            }
         } else if (axis == 1) {
            x = lowerx[0] + dx[0]
               * ((*fi)(d_dim.getValue() - 1) - loweri(d_dim.getValue() - 1) + 0.5);
            if (d_dim > tbox::Dimension(1)) {
               y = lowerx[1] + dx[1] * ((*fi)(0) - loweri(0));
            }
            if (d_dim > tbox::Dimension(2)) {
               z = lowerx[2] + dx[2] * ((*fi)(1) - loweri(1) + 0.5);
            }
         } else if (axis == 2) {
            x = lowerx[0] + dx[0] * ((*fi)(1) - loweri(1) + 0.5);
            if (d_dim > tbox::Dimension(1)) {
               y = lowerx[1] + dx[1] * ((*fi)(2) - loweri(2) + 0.5);
            }
            if (d_dim > tbox::Dimension(2)) {
               z = lowerx[2] + dx[2] * ((*fi)(0) - loweri(0));
            }
         }

         OUTERFACE_KERNEL_TYPE value;
         for (int d = 0; d < depth; ++d) {
            value = d_Dcoef + d_Acoef * x + d_Bcoef * y + d_Ccoef * z;
            if (!(tbox::MathUtilities<OUTERFACE_KERNEL_TYPE>::equalEps((*data)(*fi,
                                                                d), value))) {
               tbox::perr << "FAILED: -- patch interior not properly filled"
                          << std::endl;
            }
         }
      }
   }
}

void OuterfaceDataTest::setPhysicalBoundaryConditions(
   const hier::Patch& patch,
   const double time,
   const hier::IntVector& gcw) const
{
   NULL_USE(patch);
   NULL_USE(time);
   NULL_USE(gcw);
}

void OuterfaceDataTest::setLinearData(
   std::shared_ptr<pdat::FaceData<OUTERFACE_KERNEL_TYPE> > data,
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
   double x = 0., y = 0., z = 0.;

   const int depth = data->getDepth();

   const hier::Box sbox = data->getGhostBox() * box;

   for (tbox::Dimension::dir_t axis = 0; axis < d_dim.getValue(); ++axis) {
      const pdat::FaceIndex loweri(patch.getBox().lower(), axis, 0);
      pdat::FaceIterator fiend(pdat::FaceGeometry::end(sbox, axis));
      for (pdat::FaceIterator fi(pdat::FaceGeometry::begin(sbox, axis));
           fi != fiend; ++fi) {

         /*
          * Compute spatial location of cell center and
          * set data to linear profile.
          */

         if (axis == 0) {
            x = lowerx[0] + dx[0] * ((*fi)(0) - loweri(0));
            if (d_dim > tbox::Dimension(1)) {
               y = lowerx[1] + dx[1] * ((*fi)(1) - loweri(1) + 0.5);
            }
            if (d_dim > tbox::Dimension(2)) {
               z = lowerx[2] + dx[2] * ((*fi)(2) - loweri(2) + 0.5);
            }
         } else if (axis == 1) {
            x = lowerx[0] + dx[0]
               * ((*fi)(d_dim.getValue() - 1) - loweri(d_dim.getValue() - 1) + 0.5);
            if (d_dim > tbox::Dimension(1)) {
               y = lowerx[1] + dx[1] * ((*fi)(0) - loweri(0));
            }
            if (d_dim > tbox::Dimension(2)) {
               z = lowerx[2] + dx[2] * ((*fi)(1) - loweri(1) + 0.5);
            }
         } else if (axis == 2) {
            x = lowerx[0] + dx[0] * ((*fi)(1) - loweri(1) + 0.5);
            if (d_dim > tbox::Dimension(1)) {
               y = lowerx[1] + dx[1] * ((*fi)(2) - loweri(2) + 0.5);
            }
            if (d_dim > tbox::Dimension(2)) {
               z = lowerx[2] + dx[2] * ((*fi)(0) - loweri(0));
            }
         }

         for (int d = 0; d < depth; ++d) {
            (*data)(*fi,
                    d) = d_Dcoef + d_Acoef * x + d_Bcoef * y + d_Ccoef * z;
         }

      }
   }

}

void OuterfaceDataTest::setLinearData(
   std::shared_ptr<pdat::OuterfaceData<OUTERFACE_KERNEL_TYPE> > data,
   const hier::Box& box,
   const hier::Patch& patch) const
{
   NULL_USE(box);

   TBOX_ASSERT(data);

   std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(pgeom);
   const double* dx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x = 0., y = 0., z = 0.;

   const int depth = data->getDepth();

   for (int axis = 0; axis < d_dim.getValue(); ++axis) {
      for (int f = 0; f < 2; ++f) {
         const hier::Box databox = data->getArrayData(axis, f).getBox();

         const pdat::FaceIndex loweri(patch.getBox().lower(), axis, 0);
         hier::Box::iterator biend(databox.end());
         for (hier::Box::iterator bi(databox.begin()); bi != biend; ++bi) {

            /*
             * Compute spatial location of cell center and
             * set data to linear profile.
             */

            if (axis == 0) {
               x = lowerx[0] + dx[0] * ((*bi)(0) - loweri(0));
               if (d_dim > tbox::Dimension(1)) {
                  y = lowerx[1] + dx[1] * ((*bi)(1) - loweri(1) + 0.5);
               }
               if (d_dim > tbox::Dimension(2)) {
                  z = lowerx[2] + dx[2] * ((*bi)(2) - loweri(2) + 0.5);
               }
            } else if (axis == 1) {
               x = lowerx[0] + dx[0]
                  * ((*bi)(d_dim.getValue() - 1) - loweri(d_dim.getValue() - 1) + 0.5);
               if (d_dim > tbox::Dimension(1)) {
                  y = lowerx[1] + dx[1] * ((*bi)(0) - loweri(0));
               }
               if (d_dim > tbox::Dimension(2)) {
                  z = lowerx[2] + dx[2] * ((*bi)(1) - loweri(1) + 0.5);
               }
            } else if (axis == 2) {
               x = lowerx[0] + dx[0] * ((*bi)(1) - loweri(1) + 0.5);
               if (d_dim > tbox::Dimension(1)) {
                  y = lowerx[1] + dx[1] * ((*bi)(2) - loweri(2) + 0.5);
               }
               if (d_dim > tbox::Dimension(2)) {
                  z = lowerx[2] + dx[2] * ((*bi)(0) - loweri(0));
               }
            }
            double value = d_Dcoef + d_Acoef * x + d_Bcoef * y + d_Ccoef * z;

            for (int d = 0; d < depth; ++d) {
               data->getArrayData(axis, f) (*bi, d) = value;
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

bool OuterfaceDataTest::verifyResults(
   const hier::Patch& patch,
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number)
{
   NULL_USE(hierarchy);
   bool test_failed = false;
   if (d_do_refine || d_do_coarsen) {

      tbox::plog << "\nEntering OuterfaceDataTest::verifyResults..." << std::endl;
      tbox::plog << "level_number = " << level_number << std::endl;
      tbox::plog << "Patch box = " << patch.getBox() << std::endl;

      hier::IntVector tgcw(d_dim, 0);
      for (int i = 0; i < static_cast<int>(d_variables_dst.size()); ++i) {
         tgcw.max(patch.getPatchData(d_variables_dst[i], getDataContext())->
            getGhostCellWidth());
      }
      hier::Box pbox = patch.getBox();

      std::shared_ptr<pdat::FaceData<OUTERFACE_KERNEL_TYPE> > solution(
         new pdat::FaceData<OUTERFACE_KERNEL_TYPE>(pbox, 1, tgcw));

      hier::Box tbox(pbox);
      tbox.grow(tgcw);

      if (d_do_refine) {
         setLinearData(solution, tbox, patch);
      } else {
         setLinearData(solution, tbox, patch); //, hierarchy, level_number);
      }

      for (int i = 0; i < static_cast<int>(d_variables_dst.size()); ++i) {

         if (i % 2 == 0) {
            std::shared_ptr<pdat::FaceData<OUTERFACE_KERNEL_TYPE> > face_data(
               SAMRAI_SHARED_PTR_CAST<pdat::FaceData<OUTERFACE_KERNEL_TYPE>, hier::PatchData>(
                  patch.getPatchData(d_variables_dst[i], getDataContext())));
            TBOX_ASSERT(face_data);
            int depth = face_data->getDepth();
            hier::Box dbox = face_data->getGhostBox();

            for (tbox::Dimension::dir_t id = 0; id < d_dim.getValue(); ++id) {
               pdat::FaceIterator fiend(pdat::FaceGeometry::end(dbox, id));
               for (pdat::FaceIterator fi(pdat::FaceGeometry::begin(dbox, id));
                    fi != fiend; ++fi) {
                  OUTERFACE_KERNEL_TYPE correct = (*solution)(*fi);
                  for (int d = 0; d < depth; ++d) {
                     OUTERFACE_KERNEL_TYPE result = (*face_data)(*fi, d);
                     if (!tbox::MathUtilities<OUTERFACE_KERNEL_TYPE>::equalEps(correct,
                            result)) {
                        tbox::perr << "Test FAILED: ...."
                                   << " : face_data index = " << *fi << std::endl;
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
         } else {
            std::shared_ptr<pdat::OuterfaceData<OUTERFACE_KERNEL_TYPE> > oface_data(
               SAMRAI_SHARED_PTR_CAST<pdat::OuterfaceData<OUTERFACE_KERNEL_TYPE>, hier::PatchData>(
                  patch.getPatchData(d_variables_dst[i], getDataContext())));
            TBOX_ASSERT(oface_data);
            int depth = oface_data->getDepth();
            hier::Box dbox = oface_data->getGhostBox();

            for (tbox::Dimension::dir_t id = 0; id < d_dim.getValue(); ++id) {
               hier::Box dbox_lo(dbox);
               dbox_lo.setUpper(id, dbox_lo.lower(id));
               hier::BoxIterator loend(dbox_lo.end());
               for (hier::BoxIterator si(dbox_lo.begin()); si != loend; ++si) {
                  pdat::FaceIndex fndx(*si, id, 0);
                  OUTERFACE_KERNEL_TYPE correct = (*solution)(fndx);
                  for (int d = 0; d < depth; ++d) {
                     OUTERFACE_KERNEL_TYPE result = (*oface_data)(fndx, 0, d);
                     if (!tbox::MathUtilities<OUTERFACE_KERNEL_TYPE>::equalEps(correct,
                            result)) {
                        tbox::perr << "Test FAILED: ...."
                                   << " : oface_data index = " << fndx << std::endl;
                        tbox::perr << "    hier::Variable = "
                                   << d_variable_src_name[i]
                                   << " : depth index = " << d << std::endl;
                        tbox::perr << "    result = " << result
                                   << " : correct = " << correct << std::endl;
                        test_failed = true;
                     }
                  }
               }

               hier::Box dbox_hi(dbox);
               dbox_hi.setLower(id, dbox_hi.upper(id));
               hier::BoxIterator hiend(dbox_hi.end());
               for (hier::BoxIterator si(dbox_hi.begin()); si != hiend; ++si) {
                  pdat::FaceIndex fndx(*si, id, 1);
                  OUTERFACE_KERNEL_TYPE correct = (*solution)(fndx);
                  for (int d = 0; d < depth; ++d) {
                     OUTERFACE_KERNEL_TYPE result = (*oface_data)(fndx, 1, d);
                     if (!tbox::MathUtilities<OUTERFACE_KERNEL_TYPE>::equalEps(correct,
                            result)) {
                        tbox::perr << "Test FAILED: ...."
                                   << " : oface_data index = " << fndx << std::endl;
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
         }
      }
      if (!test_failed) {
         tbox::plog << "Outerface test Successful!" << std::endl;
      }

      solution.reset();   // just to be anal...

      tbox::plog << "\nExiting OuterfaceDataTest::verifyResults..." << std::endl;
      tbox::plog << "level_number = " << level_number << std::endl;
      tbox::plog << "Patch box = " << patch.getBox() << std::endl << std::endl;

   }

   return !test_failed;
}

}

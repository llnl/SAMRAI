/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Main program to test cell patch data operations.
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include <typeinfo>

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <memory>

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"

#include "SAMRAI/tbox/SAMRAIManager.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellDataFactory.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include <string>

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/VariableContext.h"

#include "SAMRAI/math/PatchCellDataOpsReal.h"


#include <cmath>

using namespace SAMRAI;

/* Helper function declarations */
static bool
doubleDataSameAsValue(
   int desc_id,
   double value,
   std::shared_ptr<hier::Patch> patch);

int main(
   int argc,
   char* argv[]) {

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   if (argc < 2) {
      TBOX_ERROR("Usage: " << argv[0] << " [dimension]");
   }

   const unsigned short d = static_cast<unsigned short>(atoi(argv[1]));
   TBOX_ASSERT(d > 0);
   TBOX_ASSERT(d <= SAMRAI::MAX_DIM_VAL);
   const tbox::Dimension dim(d);

   int num_failures = 0;

   const std::string log_fn = std::string("cell_patchtest.")
      + tbox::Utilities::intToString(dim.getValue(), 1) + "d.log";
   tbox::PIO::logAllNodes(log_fn);

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      /* Make a dummy mesh domain with one patch */
      double lo[SAMRAI::MAX_DIM_VAL];
      double hi[SAMRAI::MAX_DIM_VAL];
      for (int i = 0; i < dim.getValue(); ++i) {
         lo[i] = 0.0;
         if (i == 1) {
            hi[i] = 0.5;
         } else {
            hi[i] = 1.0;
         }
      }
      hier::Index indxlo(dim, 0);
      hier::Index indxhi(dim, 9);
      indxhi(1) = 4;
      hier::Box patch_box(indxlo, indxhi, hier::BlockId(0));
      hier::BoxContainer grid_domain(patch_box);
      hier::IntVector ratio(dim, 1);

      geom::CartesianGridGeometry geometry("CartesianGeometry",
                                           lo, hi, grid_domain);
      hier::ComponentSelector patch_components;

      hier::Box patch_node(patch_box, hier::LocalId::getZero(), mpi.getRank());
      std::shared_ptr<hier::Patch> tpatch(
         new hier::Patch(
            patch_node,
            hier::VariableDatabase::getDatabase()->getPatchDescriptor()));

      /* Make a variety of data on the patch. */

      /* Make three contexts for patch */
      std::shared_ptr<hier::VariableContext> ghost_width_1_context(
         hier::VariableDatabase::getDatabase()->getContext("ghost_width_1"));
      std::shared_ptr<hier::VariableContext> ghost_width_2_context(
         hier::VariableDatabase::getDatabase()->getContext("ghost_width_2"));
      std::shared_ptr<hier::VariableContext> ghost_width_3_context(
         hier::VariableDatabase::getDatabase()->getContext("ghost_width_3"));

      /* Make ghost cell IntVectors which are used when variables
       * and contexts are registered
       */
      hier::IntVector nghosts_1(dim, 1);
      hier::IntVector nghosts_2(dim, 2);
      hier::IntVector nghosts_3(dim, 3);

      /* Make cell-centered double variable for patch */
      std::shared_ptr<pdat::CellVariable<double> > cell_double_variable(
         new pdat::CellVariable<double>(
            dim,
            "cell_double_variable",
            1));

      int cdvindx[3];

      /*
       * *Register cell-centered double variable and 3 contexts with
       * hier::VariableDatabase.
       */
      cdvindx[0] =
         hier::VariableDatabase::getDatabase()->registerVariableAndContext(
            cell_double_variable, ghost_width_1_context, nghosts_1);
      patch_components.setFlag(cdvindx[0]);

      cdvindx[1] =
         hier::VariableDatabase::getDatabase()->registerVariableAndContext(
            cell_double_variable, ghost_width_2_context, nghosts_2);
      patch_components.setFlag(cdvindx[1]);

      cdvindx[2] =
         hier::VariableDatabase::getDatabase()->registerVariableAndContext(
            cell_double_variable, ghost_width_3_context, nghosts_3);
      patch_components.setFlag(cdvindx[2]);

      /* Make control volume for cell-centered patch variables */
      std::shared_ptr<hier::VariableContext> ghost_width_0_context(
         hier::VariableDatabase::getDatabase()->getContext("ghost_width_0"));
      hier::IntVector nghosts_0(dim, 0);
      std::shared_ptr<pdat::CellVariable<double> > cwgt(
         new pdat::CellVariable<double>(dim, "cwgt", 1));
      int cwgt_id =
         hier::VariableDatabase::getDatabase()->registerVariableAndContext(
            cwgt, ghost_width_0_context, nghosts_0);
      patch_components.setFlag(cwgt_id);

      int ccvindx[3];

      ccvindx[0] =
         hier::VariableDatabase::getDatabase()->registerVariableAndContext(
            cell_double_variable, ghost_width_1_context, nghosts_1);
      patch_components.setFlag(ccvindx[0]);

      ccvindx[1] =
         hier::VariableDatabase::getDatabase()->registerVariableAndContext(
            cell_double_variable, ghost_width_2_context, nghosts_2);
      patch_components.setFlag(ccvindx[1]);

      // Make two cell-centered int variables for the patch
      std::shared_ptr<pdat::CellVariable<int> > cell_int_variable(
         new pdat::CellVariable<int>(
            dim,
            "cell_int_variable",
            1));

      int civindx[3];

      civindx[0] =
         hier::VariableDatabase::getDatabase()->registerVariableAndContext(
            cell_int_variable, ghost_width_1_context, nghosts_1);
      patch_components.setFlag(civindx[0]);

      civindx[1] =
         hier::VariableDatabase::getDatabase()->registerVariableAndContext(
            cell_int_variable, ghost_width_2_context, nghosts_2);
      patch_components.setFlag(civindx[1]);

      // Test #1: Check the state of hier::PatchDescriptor
      int desc_id;
      std::string var_ctxt_name;
      bool descriptor_test_passed = true;
      bool name_error_indx[6];
      bool factory_error_indx[6];
      for (desc_id = 0; desc_id < 6; ++desc_id) {
         name_error_indx[desc_id] = false;
         factory_error_indx[desc_id] = false;
      }
      //make strings to be used for comparison during tests

      std::string cell_double_variable1("cell_double_variable##ghost_width_1");
      std::string cell_double_variable2("cell_double_variable##ghost_width_2");
      std::string cell_double_variable3("cell_double_variable##ghost_width_3");
      std::string cwgt_variable0("cwgt##ghost_width_0");
      std::string cell_int_variable1("cell_int_variable##ghost_width_1");
      std::string cell_int_variable2("cell_int_variable##ghost_width_2");

      for (desc_id = 0; desc_id < 6; ++desc_id) {
         var_ctxt_name =
            hier::VariableDatabase::getDatabase()->getPatchDescriptor()->
            mapIndexToName(desc_id);

         // test flag is used to overcome a compiler bug in GCC which
         // complained if the the comparison was done inside the if;
         // must be too complicated for it.
         bool test;

         switch (desc_id) {
            case 0:
               {

               if (var_ctxt_name != cell_double_variable1) {
                  descriptor_test_passed = false;
                  name_error_indx[desc_id] = true;
               }

               auto& pdf = *(hier::VariableDatabase::getDatabase()
                                            ->getPatchDescriptor()
                                            ->getPatchDataFactory(desc_id));

               test = ( typeid(pdf) == typeid(pdat::CellDataFactory<double>) );

               if (!test) {
                  descriptor_test_passed = false;
                  factory_error_indx[desc_id] = true;
               }

               } // END case 0
               break;

            case 1:
               {


               if (var_ctxt_name != cell_double_variable2) {
                  descriptor_test_passed = false;
                  name_error_indx[desc_id] = true;
               }

               auto& pdf = *(hier::VariableDatabase::getDatabase()
                                            ->getPatchDescriptor()
                                            ->getPatchDataFactory(desc_id));

               test = (typeid(pdf) == typeid(pdat::CellDataFactory<double>) );

               if (!test) {
                  descriptor_test_passed = false;
                  factory_error_indx[desc_id] = true;
               }

               } // END case 1
               break;

            case 2:
               {

               if (var_ctxt_name != cell_double_variable3) {
                  descriptor_test_passed = false;
                  name_error_indx[desc_id] = true;
               }

               auto& pdf = *(hier::VariableDatabase::getDatabase()
                                            ->getPatchDescriptor()
                                            ->getPatchDataFactory(desc_id));

               test = ( typeid(pdf) == typeid(pdat::CellDataFactory<double>) );

               if (!test) {
                  descriptor_test_passed = false;
                  factory_error_indx[desc_id] = true;
               }

               } // END case 2
               break;

            case 3:
               {

               if (var_ctxt_name != cwgt_variable0) {
                  descriptor_test_passed = false;
                  name_error_indx[desc_id] = true;
               }

               auto& pdf = *(hier::VariableDatabase::getDatabase()
                                            ->getPatchDescriptor()
                                            ->getPatchDataFactory(desc_id));

               test = ( typeid(pdf) == typeid(pdat::CellDataFactory<double>) );

               if (!test) {
                  descriptor_test_passed = false;
                  factory_error_indx[desc_id] = true;
               }

               } // END case 3
               break;

            case 4:
               {

               if (var_ctxt_name != cell_int_variable1) {
                  descriptor_test_passed = false;
                  name_error_indx[desc_id] = true;
               }

               auto& pdf = *(hier::VariableDatabase::getDatabase()
                                            ->getPatchDescriptor()
                                            ->getPatchDataFactory(desc_id));

               test = ( typeid(pdf) == typeid(pdat::CellDataFactory<int>) );

               if (!test) {
                  descriptor_test_passed = false;
                  factory_error_indx[desc_id] = true;
               }

               } // END case 4
               break;

            case 5:
               {
               if (var_ctxt_name != cell_int_variable2) {
                  descriptor_test_passed = false;
                  name_error_indx[desc_id] = true;
               }

               auto& pdf = *(hier::VariableDatabase::getDatabase()
                                            ->getPatchDescriptor()
                                            ->getPatchDataFactory(desc_id));

               test = ( typeid(pdf) == typeid(pdat::CellDataFactory<int>) );

               if (!test) {
                  descriptor_test_passed = false;
                  factory_error_indx[desc_id] = true;
               }

               } // END case 5
               break;
         }
      }

      if (!descriptor_test_passed) {
         ++num_failures;
         tbox::perr << "FAILED: - Test #1: State of PatchDescriptor"
                    << std::endl;

         for (desc_id = 0; desc_id < 6; ++desc_id) {
            if (name_error_indx[desc_id] == true) {
               tbox::plog << "Name for index = " << desc_id << " incorrect"
                          << std::endl;
            }
            if (factory_error_indx[desc_id] == true) {
               tbox::plog << "Factory for index = " << desc_id
                          << " incorrect" << std::endl;
            }
         }
      }

      // Test #2: Check state of hier::Patch before allocating storage
      if (!tpatch->getBox().isSpatiallyEqual(patch_box)) {
         ++num_failures;
         tbox::perr << "FAILED: - Test #2a: hier::Patch box incorrectly set\n"
                    << "Expected: d_box = " << patch_box << "\n"
                    << "Set to: d_box = " << tpatch->getBox() << std::endl;
      }
      if (tpatch->getLocalId() != patch_node.getLocalId()) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #2b: hier::Patch number incorrectly set\n"
         << "Expected: d_patch_number = "
         << patch_node.getLocalId() << "\n"
         << "Set to: d_patch_number = "
         << tpatch->getLocalId() << std::endl;
      }
/*   HOW TO GET NUMBER OF COMPONENTS ON PATCH
 *   FOR NOW JUST CHECK THAT getPatchData(0) returns NULL
 *
 *   ++num_failures;
 *   tbox::perr << "FAILED: - Test #2c: Number of components allocated incorrect\n"
 *   << "Expected: number of components = 0\n"
 *   << "Got: number of components = 0\n"
 */
      for (desc_id = 0; desc_id < 6; ++desc_id) {
         if (tpatch->checkAllocated(desc_id)) {
            ++num_failures;
            tbox::perr << "FAILED: - Test #2c." << desc_id
                       << ": Descriptor slot " << desc_id
                       << " should not be allocated but is!" << std::endl;
         }
      }

      // Allocate all data on patch
      tpatch->allocatePatchData(patch_components);

      // Test #3: Check state of hier::Patch after allocating storage
      if (!tpatch->getBox().isSpatiallyEqual(patch_box)) {
         ++num_failures;
         tbox::perr << "FAILED: - Test #3a: hier::Patch box incorrectly set\n"
                    << "Expected: d_box = " << patch_box << "\n"
                    << "Set to: d_box = " << tpatch->getBox() << std::endl;
      }
      if (tpatch->getLocalId() != patch_node.getLocalId()) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #3b: hier::Patch number incorrectly set\n"
         << "Expected: d_patch_number = "
         << patch_node.getLocalId() << "\n"
         << "Set to: d_patch_number = "
         << tpatch->getLocalId() << std::endl;
      }
/* SAME ISSUE AS ABOVE FOR NUMBER OF COMPONENTS */
      for (desc_id = 0; desc_id < 6; ++desc_id) {

         if (!tpatch->checkAllocated(desc_id)) {
            ++num_failures;
            tbox::perr << "FAILED: - Test #3c.0: Descriptor index " << desc_id
                       << " should be allocated but isn't!" << std::endl;
         } else {

            auto& p = *tpatch->getPatchData(desc_id);
            std::string patch_data_name = typeid(p).name();

            hier::IntVector ghost_width(tpatch->getPatchData(
                                           desc_id)->getGhostCellWidth());

            switch (desc_id) {
               case 0:
                  {
                  auto& p = *tpatch->getPatchData(desc_id);
                  if (typeid(p) != typeid(pdat::CellData<double>)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #3c.0.a: hier::Patch Data name incorrect\n"
                     << "Expected: pdat::CellData<double >\n"
                     << "Actual: " << patch_data_name << std::endl;
                  }
                  if (ghost_width != hier::IntVector(dim, 1)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #3c.0.b: Ghost width incorrect\n"
                     << "Expected: (1,1)\n"
                     << "Actual: " << ghost_width << std::endl;
                  }

                  } // END case 0
                  break;

               case 1:
                  {
                  auto& p = *tpatch->getPatchData(desc_id);
                  if (typeid(p) != typeid(pdat::CellData<double>)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #3c.1.a: hier::Patch Data name incorrect\n"
                     << "Expected: pdat::CellData<double >\n"
                     << "Actual: " << patch_data_name << std::endl;
                  }
                  if (ghost_width != hier::IntVector(dim, 2)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #3c.1.b: Ghost width incorrect\n"
                     << "Expected: (2,2)\n"
                     << "Actual: " << ghost_width << std::endl;
                  }

                  } // END case 1
                  break;

               case 2:
                  {
                  auto& p = *tpatch->getPatchData(desc_id);
                  if (typeid(p) != typeid(pdat::CellData<double>)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #3c.2.a: hier::Patch Data name incorrect\n"
                     << "Expected: pdat::CellData<double >\n"
                     << "Actual: " << patch_data_name << std::endl;
                  }
                  if (ghost_width != hier::IntVector(dim, 3)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #3c.2.b: Ghost width incorrect\n"
                     << "Expected: (3,3)\n"
                     << "Actual: " << ghost_width << std::endl;
                  }

                  } // END case 2
                  break;

               case 3:
                  {
                  auto& p = *tpatch->getPatchData(desc_id);
                  if (typeid(p) != typeid(pdat::CellData<double>)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #3c.3.a: hier::Patch Data name incorrect\n"
                     << "Expected: pdat::CellData<double >\n"
                     << "Actual: " << patch_data_name << std::endl;
                  }
                  if (ghost_width != hier::IntVector(dim, 0)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #3c.3.b: Ghost width incorrect\n"
                     << "Expected: (0,0)\n"
                     << "Actual: " << ghost_width << std::endl;
                  }

                  } // END case 3
                  break;

               case 4:
                  {
                  auto& p = *tpatch->getPatchData(desc_id);
                  if (typeid(p) != typeid(pdat::CellData<int>)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #3c.4.a: hier::Patch Data name incorrect\n"
                     << "Expected: pdat::CellData<int >\n"
                     << "Actual: " << patch_data_name << std::endl;
                  }
                  if (ghost_width != hier::IntVector(dim, 1)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #3c.4.b: Ghost width incorrect\n"
                     << "Expected: (1,1)\n"
                     << "Actual: " << ghost_width << std::endl;
                  }

                  } // END case 4
                  break;

               case 5:
                  {
                  auto& p = *tpatch->getPatchData(desc_id);
                  if (typeid(p) != typeid(pdat::CellData<int>)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #3c.5.a: hier::Patch Data name incorrect\n"
                     << "Expected: pdat::CellData<int >\n"
                     << "Actual: " << patch_data_name << std::endl;
                  }
                  if (ghost_width != hier::IntVector(dim, 2)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #3c.5.b: Ghost width incorrect\n"
                     << "Expected: (2,2)\n"
                     << "Actual: " << ghost_width << std::endl;
                  }

                  } // END case 5
                  break;

            }
         }
      }

      // Initialize control volume data for cell-centered data
      const double* dx = geometry.getDx();
      double cell_vol = dx[0];
      for (int i = 1; i < dim.getValue(); ++i) {
         cell_vol *= dx[i];
      }

      std::shared_ptr<pdat::CellData<double> > weight(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            tpatch->getPatchData(cwgt_id)));
      TBOX_ASSERT(weight);
      weight->fillAll(cell_vol);

      // Simple tests of cell data operations

      math::PatchCellDataOpsReal<double> cdops_double;

      // Get pointers to patch data objects
      std::shared_ptr<pdat::CellData<double> > cddata0(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            tpatch->getPatchData(cdvindx[0])));
      std::shared_ptr<pdat::CellData<double> > cddata1(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            tpatch->getPatchData(cdvindx[1])));
      std::shared_ptr<pdat::CellData<double> > cddata2(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            tpatch->getPatchData(cdvindx[2])));

      TBOX_ASSERT(cddata0);
      TBOX_ASSERT(cddata1);
      TBOX_ASSERT(cddata2);

      // Test #4a: math::PatchCellDataOpsReal::setToScalar()
      // Expected: cddata0 = 0.0
      cdops_double.setToScalar(cddata0, 0.0, cddata0->getGhostBox());
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      double val0 = 0.0;
      if (!doubleDataSameAsValue(cdvindx[0], val0, tpatch)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #4a: math::PatchCellDataOpsReal::setToScalar()\n"
         << "Expected: cddata0 = " << val0 << std::endl;
         cdops_double.printData(cddata0, cddata0->getGhostBox(), tbox::plog);
      }

      // Test #4b: math::PatchCellDataOpsReal::setToScalar()
      // Expected: cddata1 = 1.0
      cdops_double.setToScalar(cddata1, 1.0, cddata1->getGhostBox());
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      double val1 = 1.0;
      if (!doubleDataSameAsValue(cdvindx[1], val1, tpatch)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #4b: math::PatchCellDataOpsReal::setToScalar()\n"
         << "Expected: cddata1 = " << val1 << std::endl;
         cdops_double.printData(cddata1, cddata1->getGhostBox(), tbox::plog);
      }

      // Test #4c: math::PatchCellDataOpsReal::setToScalar()
      // Expected: cddata2 = 2.0
      cdops_double.setToScalar(cddata2, 2.0, cddata2->getGhostBox());
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      double val2 = 2.0;
      if (!doubleDataSameAsValue(cdvindx[2], val2, tpatch)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #4b: math::PatchCellDataOpsReal::setToScalar()\n"
         << "Expected: cddata2 = " << val2 << std::endl;
         cdops_double.printData(cddata2, cddata2->getGhostBox(), tbox::plog);
      }

      // Test #5: math::PatchCellDataOpsReal::add()
      // Expected: cddata0 =  cddata1 + cddata2
      cdops_double.add(cddata0, cddata1, cddata2, tpatch->getBox());
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      double val_add = 3.0;
      if (!doubleDataSameAsValue(cdvindx[0], val_add, tpatch)) {
         ++num_failures;
         tbox::perr << "FAILED: - Test #5: math::PatchCellDataOpsReal::add()\n"
                    << "Expected: cddata0 = " << val_add << std::endl;
         cdops_double.printData(cddata0, cddata0->getGhostBox(), tbox::plog);
      }

      // Test #6: math::PatchCellDataOpsReal::subtract() on [(3,1),(5,2)]
      // Expected: cddata0 = cddata0 - cddata2
      hier::Index indx0(dim, 1);
      hier::Index indx1(dim, 2);
      indx0(0) = 3;
      indx1(0) = 5;
      cdops_double.subtract(cddata0, cddata0, cddata2,
         hier::Box(indx0, indx1, hier::BlockId(0)));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      bool subtract_inbox_test_passed = true;
      hier::Box inbox(indx0, indx1, hier::BlockId(0));
      double val_inbox = 1.0;
      double val_not_inbox = 3.0;
      std::shared_ptr<pdat::CellData<double> > cvdata(
         SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
            tpatch->getPatchData(cdvindx[0])));

      TBOX_ASSERT(cvdata);

      pdat::CellIterator cend(pdat::CellGeometry::end(cvdata->getBox()));
      for (pdat::CellIterator c(pdat::CellGeometry::begin(cvdata->getBox()));
           c != cend && subtract_inbox_test_passed; ++c) {
         pdat::CellIndex cell_index = *c;

         double value;
         if (inbox.contains(cell_index)) {
            value = val_inbox;
         } else {
            value = val_not_inbox;
         }

         if (!tbox::MathUtilities<double>::equalEps((*cvdata)(cell_index),
                value)) {
            subtract_inbox_test_passed = false;
         }
      }

      if (!subtract_inbox_test_passed) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #6: math::PatchCellDataOpsReal::subtract() on [(3,1),(5,2)]\n"
         << "Expected: cddata0 = 1.0 in [(3,1),(5,2)]\n"
         << "          cddata0 = 3.0 outside box\n" << std::endl;
         cdops_double.printData(cddata0, tpatch->getBox(), tbox::plog);
      }

      // Test #7: math::PatchCellDataOpsReal::scale()
      // Expected: cddata0 = 0.4 * cddata2
      cdops_double.scale(cddata0, 0.4, cddata2, tpatch->getBox());
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      double val_scale = 0.8;
      if (!doubleDataSameAsValue(cdvindx[0], val_scale, tpatch)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #7: math::PatchCellDataOpsReal::scale()\n"
         << "Expected: cddata0 = " << val_scale << std::endl;
         cdops_double.printData(cddata0, cddata0->getGhostBox(), tbox::plog);
      }

      // Test #8: math::PatchCellDataOpsReal::multiply()
      // Expected: cddata0 = cddata0 * cddata2
      cdops_double.multiply(cddata0, cddata0, cddata2, tpatch->getBox());
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      double val_mult = 1.6;
      if (!doubleDataSameAsValue(cdvindx[0], val_mult, tpatch)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #8: math::PatchCellDataOpsReal::multiply()\n"
         << "Expected: cddata0 = " << val_mult << std::endl;
         cdops_double.printData(cddata0, cddata0->getGhostBox(), tbox::plog);
      }

      // Test #9: math::PatchCellDataOpsReal::divide() in box [(3,1),(5,2)]
      // Expected: cddata0 = cddata0/cddata2
      cdops_double.divide(cddata0, cddata0, cddata2,
         hier::Box(indx0, indx1, hier::BlockId(0)));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      bool divide_inbox_test_passed = true;
      val_inbox = 0.8;
      val_not_inbox = 1.6;
      cvdata = SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>,
                          hier::PatchData>(tpatch->getPatchData(cdvindx[0]));
      TBOX_ASSERT(cvdata);

      pdat::CellIterator ccend(pdat::CellGeometry::end(cvdata->getBox()));
      for (pdat::CellIterator cc(pdat::CellGeometry::begin(cvdata->getBox()));
           cc != ccend && divide_inbox_test_passed; ++cc) {
         pdat::CellIndex cell_index = *cc;

         double value;
         if (inbox.contains(cell_index)) {
            value = val_inbox;
         } else {
            value = val_not_inbox;
         }

         if (!tbox::MathUtilities<double>::equalEps((*cvdata)(cell_index),
                value)) {
            divide_inbox_test_passed = false;
         }
      }

      if (!divide_inbox_test_passed) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #9: math::PatchCellDataOpsReal::divide() on [(3,1),(5,2)]\n"
         << "Expected: cddata0 = 1.0 in [(3,1),(5,2)]\n"
         << "          cddata0 = 3.0 outside box\n" << std::endl;
         cdops_double.printData(cddata0, tpatch->getBox(), tbox::plog);
      }

      // Test #10: math::PatchCellDataOpsReal::reciprocal()
      // Expected: cddata0 = 1/cddata2
      cdops_double.reciprocal(cddata0, cddata2, tpatch->getBox());
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      double val_rec = 0.5;
      if (!doubleDataSameAsValue(cdvindx[0], val_rec, tpatch)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #10: math::PatchCellDataOpsReal::reciprocal()\n"
         << "Expected: cddata0 = " << val_rec << std::endl;
         cdops_double.printData(cddata0, cddata0->getGhostBox(), tbox::plog);
      }

      // Reset cddata1 and cddata2
      cdops_double.setToScalar(cddata1, 1.0, cddata1->getGhostBox());
      cdops_double.setToScalar(cddata2, 2.0, cddata2->getGhostBox());

#if defined(HAVE_RAJA)
      tbox::parallel_synchronize(); 
#endif
      // Test #11: math::PatchCellDataOpsReal::linearSum()
      // Expected: cddata0 = 10*cddata1 + 20*cddata2
      cdops_double.linearSum(cddata0,
         10.0,
         cddata1,
         20.0,
         cddata2,
         tpatch->getBox());
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      double val_linearSum = 50.0;
      if (!doubleDataSameAsValue(cdvindx[0], val_linearSum, tpatch)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #11: math::PatchCellDataOpsReal::linearSum()\n"
         << "Expected: cddata0 = " << val_linearSum << std::endl;
         cdops_double.printData(cddata0, cddata0->getGhostBox(), tbox::plog);
      }

      cdops_double.setRandomValues(cddata0, 1.0, 0.001, hier::Box(indx0, indx1, hier::BlockId(0)));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      tbox::plog << "\ncddata0 = random " << std::endl;
      cdops_double.printData(cddata0, hier::Box(indx0, indx1, hier::BlockId(0)), tbox::plog);

      cdops_double.setRandomValues(cddata0, 1.0, 0.001, hier::Box(indx0, indx1, hier::BlockId(0)));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      tbox::plog << "\ncddata0 = random " << std::endl;
      cdops_double.printData(cddata0, hier::Box(indx0, indx1, hier::BlockId(0)), tbox::plog);

      // Reset cddata0 to 0.0
      cdops_double.setToScalar(cddata0, 0.0, cddata0->getGhostBox());
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif

      // Test #12: math::PatchCellDataOpsReal::linearSum() on box = [(3,1),(5,2)]
      // Expected: cddata0 = 10*cddata1 + 20*cddata2 on [(3,1),(5,2)]
      cdops_double.linearSum(cddata0, 10.0, cddata1, 20.0, cddata2,
         hier::Box(indx0, indx1, hier::BlockId(0)));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      bool restricted_linSum_test_passed = true;
      val_inbox = 50.0;
      val_not_inbox = 0.0;
      cvdata = SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>,
                          hier::PatchData>(tpatch->getPatchData(cdvindx[0]));
      TBOX_ASSERT(cvdata);

      pdat::CellIterator cciend(pdat::CellGeometry::end(cvdata->getBox()));
      for (pdat::CellIterator cci(pdat::CellGeometry::begin(cvdata->getBox()));
           cci != cciend && restricted_linSum_test_passed; ++cci) {
         pdat::CellIndex cell_index = *cci;

         double value;
         if (inbox.contains(cell_index)) {
            value = val_inbox;
         } else {
            value = val_not_inbox;
         }

         if (!tbox::MathUtilities<double>::equalEps((*cvdata)(cell_index),
                value)) {
            restricted_linSum_test_passed = false;
         }
      }

      if (!restricted_linSum_test_passed) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #12: math::PatchCellDataOpsReal::linearSum()\n"
         << "Expected: cddata0 = " << val_linearSum
         << " on box = [(3,1),(5,2)]" << std::endl;
         cdops_double.printData(cddata0, tpatch->getBox(), tbox::plog);
      }

// set individual data points and check min/max routines
      hier::Index newindx0(indx0);
      hier::Index newindx1(indx0);
      newindx0(1) = 2;
      cdops_double.setToScalar(cddata1, 0.0003, hier::Box(indx0, newindx0, hier::BlockId(0)));
      newindx0(0) = 1;
      cdops_double.setToScalar(cddata1, 12345.0, hier::Box(newindx0, newindx0, hier::BlockId(0)));
      newindx0(0) = 5;
      newindx0(1) = 3;
      newindx1(0) = 5;
      newindx1(1) = 4;
      cdops_double.setToScalar(cddata1, 21.0, hier::Box(newindx0, newindx1, hier::BlockId(0)));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif

      // Test #13: math::PatchCellDataOpsReal::setToScalar() on box
      // Expected: cddata1 = 0.0003 in [(3,1),(3,2)]
      //           cddata1 = 12345.0 in [(1,2),(1,2)]
      //           cddata1 = 21.0 in [(5,3),(5,4)]
      //           cddata1 = 1.0 everywhere else
      bool setToScalar_onBox_test_passed = true;
      newindx0 = indx0;
      newindx0(1) = 2;
      hier::Box box1(indx0, newindx0, hier::BlockId(0));
      newindx0(0) = 1;
      hier::Box box2(newindx0, newindx0, hier::BlockId(0));
      newindx0(0) = 5;
      newindx0(1) = 3;
      newindx1(0) = 5;
      newindx1(1) = 4;
      hier::Box box3(newindx0, newindx1, hier::BlockId(0));
      double val_inbox1 = 0.0003;
      double val_inbox2 = 12345.0;
      double val_inbox3 = 21.0;
      val_not_inbox = 1.0;

      cvdata = SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>,
                          hier::PatchData>(tpatch->getPatchData(cdvindx[1]));
      TBOX_ASSERT(cvdata);
      pdat::CellIterator ciend(pdat::CellGeometry::end(cvdata->getBox()));
      for (pdat::CellIterator ci(pdat::CellGeometry::begin(cvdata->getBox()));
           ci != ciend && setToScalar_onBox_test_passed; ++ci) {
         pdat::CellIndex cell_index = *ci;

         double value;
         if (box1.contains(cell_index)) {
            value = val_inbox1;
         } else {
            if (box2.contains(cell_index)) {
               value = val_inbox2;
            } else {
               if (box3.contains(cell_index)) {
                  value = val_inbox3;
               } else {
                  value = val_not_inbox;
               }
            }
         }

         if (!tbox::MathUtilities<double>::equalEps((*cvdata)(cell_index),
                value)) {
            setToScalar_onBox_test_passed = false;
         }
      }

      if (!setToScalar_onBox_test_passed) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #13: math::PatchCellDataOpsReal::setToScalar() on box\n"
         << "Expected: cddata1 = 0.0003 in [(3,1),(3,2)]\n"
         << "          cddata1 = 12345.0 in [(1,2),(1,2)]\n"
         << "          cddata1 = 21.0 in [(5,3),(5,4)]\n"
         << "          cddata1 = 1.0 everywhere else\n" << std::endl;
         cdops_double.printData(cddata1, tpatch->getBox(), tbox::plog);
      }

      // Test #14: math::PatchCellDataOpsReal::max() on box [(3,1),(7,4)]
      // Expected: lmax = 21.0
      hier::Index indx2(dim, 4);
      indx2(0) = 7;
      double lmax = cdops_double.max(cddata1, hier::Box(indx0, indx2, hier::BlockId(0)));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(lmax, 21.0)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #14: math::PatchCellDataOpsReal::max() on box [(3,1),(7,4)]\n"
         << "Expected value = 21.0, Computed value = "
         << lmax << std::endl;
      }

      // Test #15: math::PatchCellDataOpsReal::max() in box [(0,0),(9,4)]
      // Expected: lmax = 12345.0
      lmax = cdops_double.max(cddata1, tpatch->getBox());
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(lmax, 12345.0)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #15: math::PatchCellDataOpsReal::max() in box [(0,0),(9,4)]\n"
         << "Expected value = 12345.0, Computed value = "
         << lmax << std::endl;
      }

// check axpy, axmy routines

      // set cddata0, cddata1, cddata2
      cdops_double.setToScalar(cddata0, 0.0, cddata0->getGhostBox());
      cdops_double.setToScalar(cddata1, 1.0, cddata1->getGhostBox());
      cdops_double.setToScalar(cddata2, 2.0, cddata2->getGhostBox());

#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      // Test #16: math::PatchCellDataOpsReal::axpy()
      // Expected: cddata0 = 0.5 * 1.0 + 2.0 = 2.5
      cdops_double.axpy(cddata0, 0.5, cddata1, cddata2, tpatch->getBox());
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      double val_axpy = 2.5;
      if (!doubleDataSameAsValue(cdvindx[0], val_axpy, tpatch)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #16: math::PatchCellDataOpsReal::axpy()\n"
         << "Expected: cddata0 = " << val_axpy << std::endl;
         cdops_double.printData(cddata0, cddata0->getBox(), tbox::plog);
      }

      // Test #17: math::PatchCellDataOpsReal::axmy()
      // Expected: cddata0 = 1.5 * 2.0 - 1.0 = 2.0
      cdops_double.axmy(cddata0, 1.5, cddata2, cddata1, tpatch->getBox());
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      double val_axmy = 2.0;
      if (!doubleDataSameAsValue(cdvindx[0], val_axmy, tpatch)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #17: math::PatchCellDataOpsReal::axmy()\n"
         << "Expected: cddata0 = " << val_axmy << std::endl;
         cdops_double.printData(cddata0, cddata0->getBox(), tbox::plog);
      }

// Test the norm ops stuff

      // Test #18a: math::PatchCellDataOpsReal::sumControlVolumes() for cddata1
      // Expected: lsum = 0.5
      double lsum = cdops_double.sumControlVolumes(cddata1,
            weight,
            tpatch->getBox());
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(lsum, 0.5)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #18a: math::PatchCellDataOpsReal::sumControlVolumes() for cddata1\n"
         << "Expected value = 0.5, Computed value = "
         << lsum << std::endl;
      }

      // Test #18b: math::PatchCellDataOpsReal::sumControlVolumes() for cddata2
      // Expected: lsum = 0.5
      lsum = cdops_double.sumControlVolumes(cddata2, weight, tpatch->getBox());
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(lsum, 0.5)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #18b: math::PatchCellDataOpsReal::sumControlVolumes() for cddata2\n"
         << "Expected value = 0.5, Computed value = "
         << lsum << std::endl;
      }

      // Test #19a: math::PatchCellDataOpsReal::L1norm() for cddata1
      // Expected: l1norm = 0.5
      double l1norm = cdops_double.L1Norm(cddata1, tpatch->getBox(), weight);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(l1norm, 0.5)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #19a: math::PatchCellDataOpsReal::L1norm() for cddata1\n"
         << "Expected value = 0.5, Computed value = "
         << l1norm << std::endl;
      }

      // Test #19b: math::PatchCellDataOpsReal::L1norm() for cddata2
      // Expected: l1norm = 1.0
      l1norm = cdops_double.L1Norm(cddata2, tpatch->getBox(), weight);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(l1norm, 1.0)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #19b: math::PatchCellDataOpsReal::L1norm() for cddata2\n"
         << "Expected value = 1.0, Computed value = "
         << l1norm << std::endl;
      }

      // Test #20: math::PatchCellDataOpsReal::L2norm() for cddata2
      // Expected: l2norm = sqrt(2) = 1.4142135623731
      double l2norm = cdops_double.L2Norm(cddata2, tpatch->getBox(), weight);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(l2norm, 1.4142135623731)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #20: math::PatchCellDataOpsReal::L2norm() for cddata2\n"
         << "Expected value = sqrt(2) = 1.4142135623731, Computed value = "
         << std::setprecision(12) << l2norm << std::endl;
      }

      // Reset cddata1 to 0.5
      cdops_double.setToScalar(cddata1, 0.5, tpatch->getBox());
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif

      // Test #21: math::PatchCellDataOpsReal::weightedL2norm() for cddata2
      // Expected: wl2norm = sqrt(0.5) = 0.70710678118655
      double wl2norm = cdops_double.weightedL2Norm(cddata2,
            cddata1,
            tpatch->getBox(),
            weight);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(wl2norm, 0.70710678118655)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #21: math::PatchCellDataOpsReal::weightedL2norm() for cddata2\n"
         << "Expected value = sqrt(0.5) = 0.70710678118655, Computed value = "
         << wl2norm << std::endl;
      }

      // Test #22: math::PatchCellDataOpsReal::RMSNorm() for cddata2
      // Expected: rmsnorm= L2-Norm/sqrt(control volume) = 2.0
      double rmsnorm = cdops_double.RMSNorm(cddata2, tpatch->getBox(), weight);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(rmsnorm, 2.0)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #22: math::PatchCellDataOpsReal::RMSNorm() for cddata2\n"
         << "Expected value = L2-Norm/sqrt(control volume) = 2.0, "
         << "Computed value = " << rmsnorm << std::endl;
      }

      // Test #23: math::PatchCellDataOpsReal::weightedRMSNorm() for cddata2
      // Expected: wrmsnorm= Weighted L2-Norm/sqrt(control volume) = 1.0
      double wrmsnorm = cdops_double.weightedRMSNorm(cddata2,
            cddata1,
            tpatch->getBox(),
            weight);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(wrmsnorm, 1.0)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #23: math::PatchCellDataOpsReal::weightedRMSNorm() for cddata2\n"
         << "Expected value = Weighted L2-Norm/sqrt(control volume) = 1.0, "
         << "Computed value = " << wrmsnorm << std::endl;
      }

      // Test #24: math::PatchCellDataOpsReal::maxNorm() for cddata2
      // Expected: maxnorm = 2.0
      double maxnorm = cdops_double.maxNorm(cddata2, tpatch->getBox(), weight);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(maxnorm, 2.0)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #24: math::PatchCellDataOpsReal::maxNorm() for cddata2\n"
         << "Expected value = 2.0, Computed value = "
         << maxnorm << std::endl;
      }

      // Reset cddata1 and cddata2
      cdops_double.setToScalar(cddata1, 5.0, cddata1->getGhostBox());
      cdops_double.setToScalar(cddata2, 3.0, cddata2->getGhostBox());

#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      // Test #25: math::PatchCellDataOpsReal::dotp() - (cddata1) * (cddata2)
      // Expected: dotp = 7.5
      double dotp = cdops_double.dot(cddata1, cddata2, tpatch->getBox(), weight);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(dotp, 7.5)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #25: math::PatchCellDataOpsReal::dotp() - (cddata1) * (cddata2)\n"
         << "Expected value = 7.5, Computed value = "
         << dotp << std::endl;
      }

      // Test #26: Check state of hier::Patch before deallocating storage
      if (!tpatch->getBox().isSpatiallyEqual(patch_box)) {
         ++num_failures;
         tbox::perr << "FAILED: - Test #26a: hier::Patch box incorrectly set\n"
                    << "Expected: d_box = " << patch_box << "\n"
                    << "Set to: d_box = " << tpatch->getBox() << std::endl;
      }
      if (tpatch->getLocalId() != patch_node.getLocalId()) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #26b: hier::Patch number incorrectly set\n"
         << "Expected: d_patch_number = "
         << patch_node.getLocalId() << "\n"
         << "Set to: d_patch_number = "
         << tpatch->getLocalId() << std::endl;
      }
/* SAME ISSUE AS ABOVE FOR NUMBER OF COMPONENTS */
      for (desc_id = 0; desc_id < 6; ++desc_id) {

         if (!tpatch->checkAllocated(desc_id)) {
            ++num_failures;
            tbox::perr << "FAILED: - Test #26c.0: Descriptor index " << desc_id
                       << " should be allocated but isn't!" << std::endl;
         } else {

            auto& p = *tpatch->getPatchData(desc_id);
            std::string patch_data_name = typeid(p).name();

            hier::IntVector ghost_width(tpatch->getPatchData(
                                           desc_id)->getGhostCellWidth());

            switch (desc_id) {
               case 0:
                  {
                  auto& p = *tpatch->getPatchData(desc_id);
                  if (typeid(p) != typeid(pdat::CellData<double>)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #26c.0.a: hier::Patch Data name incorrect\n"
                     << "Expected: pdat::CellData<double >\n"
                     << "Actual: " << patch_data_name << std::endl;
                  }
                  if (ghost_width != hier::IntVector(dim, 1)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #26c.0.b: Ghost width incorrect\n"
                     << "Expected: (1,1)\n"
                     << "Actual: " << ghost_width << std::endl;
                  }

                  } // END case 0
                  break;

               case 1:
                  {
                  auto& p = *tpatch->getPatchData(desc_id);
                  if (typeid(p) != typeid(pdat::CellData<double>)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #26c.1.a: hier::Patch Data name incorrect\n"
                     << "Expected: pdat::CellData<double >\n"
                     << "Actual: " << patch_data_name << std::endl;
                  }
                  if (ghost_width != hier::IntVector(dim, 2)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #26c.1.b: Ghost width incorrect\n"
                     << "Expected: (2,2)\n"
                     << "Actual: " << ghost_width << std::endl;
                  }

                  } // END case 1
                  break;

               case 2:
                  {
                  auto& p = *tpatch->getPatchData(desc_id);
                  if (typeid(p) != typeid(pdat::CellData<double>)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #26c.2.a: hier::Patch Data name incorrect\n"
                     << "Expected: pdat::CellData<double >\n"
                     << "Actual: " << patch_data_name << std::endl;
                  }
                  if (ghost_width != hier::IntVector(dim, 3)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #26c.2.b: Ghost width incorrect\n"
                     << "Expected: (3,3)\n"
                     << "Actual: " << ghost_width << std::endl;
                  }

                  } // END case 2
                  break;

               case 3:
                  {
                  auto& p = *tpatch->getPatchData(desc_id);
                  if (typeid(p) != typeid(pdat::CellData<double>)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #26c.3.a: hier::Patch Data name incorrect\n"
                     << "Expected: pdat::CellData<double >\n"
                     << "Actual: " << patch_data_name << std::endl;
                  }
                  if (ghost_width != hier::IntVector(dim, 0)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #26c.3.b: Ghost width incorrect\n"
                     << "Expected: (0,0)\n"
                     << "Actual: " << ghost_width << std::endl;
                  }

                  } // END case 3
                  break;

               case 4:
                  {
                  auto& p = *tpatch->getPatchData(desc_id);
                  if (typeid(p) != typeid(pdat::CellData<int>)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #26c.4.a: hier::Patch Data name incorrect\n"
                     << "Expected: pdat::CellData<int >\n"
                     << "Actual: " << patch_data_name << std::endl;
                  }
                  if (ghost_width != hier::IntVector(dim, 1)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #26c.4.b: Ghost width incorrect\n"
                     << "Expected: (1,1)\n"
                     << "Actual: " << ghost_width << std::endl;
                  }

                  } // END case 4
                  break;

               case 5:
                  {
                  auto& p = *tpatch->getPatchData(desc_id);
                  if (typeid(p) != typeid(pdat::CellData<int>)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #26c.5.a: hier::Patch Data name incorrect\n"
                     << "Expected: pdat::CellData<int >\n"
                     << "Actual: " << patch_data_name << std::endl;
                  }
                  if (ghost_width != hier::IntVector(dim, 2)) {
                     ++num_failures;
                     tbox::perr
                     << "FAILED: - Test #26c.5.b: Ghost width incorrect\n"
                     << "Expected: (2,2)\n"
                     << "Actual: " << ghost_width << std::endl;
                  }

                  } // END case 5
                  break;

            }
         }
      }

      // Deallocate all data on patch
      tpatch->deallocatePatchData(patch_components);

      // Test #27: Check state of hier::Patch after deallocating storage
      if (!tpatch->getBox().isSpatiallyEqual(patch_box)) {
         ++num_failures;
         tbox::perr << "FAILED: - Test #27a: hier::Patch box incorrectly set\n"
                    << "Expected: d_box = " << patch_box << "\n"
                    << "Set to: d_box = " << tpatch->getBox() << std::endl;
      }
      if (tpatch->getLocalId() != patch_node.getLocalId()) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #27b: hier::Patch number incorrectly set\n"
         << "Expected: d_patch_number = "
         << patch_node.getLocalId() << "\n"
         << "Set to: d_patch_number = "
         << tpatch->getLocalId() << std::endl;
      }
/* SAME ISSUE AS ABOVE FOR NUMBER OF COMPONENTS */
      for (desc_id = 0; desc_id < 6; ++desc_id) {

         if (tpatch->checkAllocated(desc_id)) {
            ++num_failures;
            tbox::perr << "FAILED: - Test #27c: Descriptor index " << desc_id
                       << " should be deallocated but isn't!" << std::endl;
         }
      }

      cwgt.reset();
      cell_double_variable.reset();
      cell_int_variable.reset();

      if (num_failures == 0) {
         tbox::pout << "\nPASSED:  cell patchtest" << std::endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return num_failures;
}

/*
 * Returns true if all the data in the patch is equal to the specified
 * value.  Returns false otherwise.
 */
static bool
doubleDataSameAsValue(
   int desc_id,
   double value,
   std::shared_ptr<hier::Patch> patch)
{
   bool test_passed = true;

   std::shared_ptr<pdat::CellData<double> > cvdata(
      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
         patch->getPatchData(desc_id)));

   TBOX_ASSERT(cvdata);

   pdat::CellIterator cend(pdat::CellGeometry::end(cvdata->getBox()));
   for (pdat::CellIterator c(pdat::CellGeometry::begin(cvdata->getBox()));
        c != cend && test_passed; ++c) {
      pdat::CellIndex cell_index = *c;
      if (!tbox::MathUtilities<double>::equalEps((*cvdata)(cell_index),
             value)) {
         test_passed = false;
      }
   }

   return test_passed;
}

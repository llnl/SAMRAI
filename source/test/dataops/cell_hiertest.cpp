/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Main program to test cell-centered patch data ops
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

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
#include "SAMRAI/math/HierarchyDataOpsReal.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/VariableContext.h"


using namespace SAMRAI;

/* Helper function declarations */
static bool
doubleDataSameAsValue(
   int desc_id,
   double value,
   std::shared_ptr<hier::PatchHierarchy> hierarchy);

#define NVARS 4

int main(
   int argc,
   char* argv[]) {
   int num_failures = 0;

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

   const std::string log_fn = std::string("cell_hiertest.")
      + tbox::Utilities::intToString(dim.getValue(), 1) + "d.log";
   tbox::PIO::logAllNodes(log_fn);

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      int ln, iv;

      // Make a dummy hierarchy domain
      double lo[SAMRAI::MAX_DIM_VAL];
      double hi[SAMRAI::MAX_DIM_VAL];

      hier::Index clo0(dim);
      hier::Index chi0(dim);
      hier::Index clo1(dim);
      hier::Index chi1(dim);
      hier::Index flo0(dim);
      hier::Index fhi0(dim);
      hier::Index flo1(dim);
      hier::Index fhi1(dim);

      for (int i = 0; i < dim.getValue(); ++i) {
         lo[i] = 0.0;
         clo0(i) = 0;
         flo0(i) = 4;
         fhi0(i) = 7;
         if (i == 1) {
            hi[i] = 0.5;
            chi0(i) = 2;
            clo1(i) = 3;
            chi1(i) = 4;
         } else {
            hi[i] = 1.0;
            chi0(i) = 9;
            clo1(i) = 0;
            chi1(i) = 9;
         }
         if (i == 0) {
            flo1(i) = 8;
            fhi1(i) = 13;
         } else {
            flo1(i) = flo0(i);
            fhi1(i) = fhi0(i);
         }
      }

      hier::Box coarse0(clo0, chi0, hier::BlockId(0));
      hier::Box coarse1(clo1, chi1, hier::BlockId(0));
      hier::Box fine0(flo0, fhi0, hier::BlockId(0));
      hier::Box fine1(flo1, fhi1, hier::BlockId(0));
      hier::IntVector ratio(dim, 2);

      hier::BoxContainer coarse_domain;
      hier::BoxContainer fine_boxes;
      coarse_domain.pushBack(coarse0);
      coarse_domain.pushBack(coarse1);
      fine_boxes.pushBack(fine0);
      fine_boxes.pushBack(fine1);

      std::shared_ptr<geom::CartesianGridGeometry> geometry(
         new geom::CartesianGridGeometry(
            "CartesianGeometry",
            lo,
            hi,
            coarse_domain));

      std::shared_ptr<hier::PatchHierarchy> hierarchy(
         new hier::PatchHierarchy("PatchHierarchy", geometry));

      hierarchy->setMaxNumberOfLevels(2);
      hierarchy->setRatioToCoarserLevel(ratio, 1);

      const int nproc = mpi.getSize();

      const int n_coarse_boxes = coarse_domain.size();
      const int n_fine_boxes = fine_boxes.size();

      std::shared_ptr<hier::BoxLevel> layer0(
         std::make_shared<hier::BoxLevel>(
            hier::IntVector(dim, 1), geometry));
      std::shared_ptr<hier::BoxLevel> layer1(
         std::make_shared<hier::BoxLevel>(ratio, geometry));

      hier::BoxContainer::iterator coarse_itr = coarse_domain.begin();
      for (int ib = 0; ib < n_coarse_boxes; ++ib, ++coarse_itr) {
         if (nproc > 1) {
            if (ib == layer0->getMPI().getRank()) {
               layer0->addBox(hier::Box(*coarse_itr, hier::LocalId(ib),
                     layer0->getMPI().getRank()));
            }
         } else {
            layer0->addBox(hier::Box(*coarse_itr, hier::LocalId(ib), 0));
         }
      }

      hier::BoxContainer::iterator fine_itr = fine_boxes.begin();
      for (int ib = 0; ib < n_fine_boxes; ++ib, ++fine_itr) {
         if (nproc > 1) {
            if (ib == layer1->getMPI().getRank()) {
               layer1->addBox(hier::Box(*fine_itr, hier::LocalId(ib),
                     layer1->getMPI().getRank()));
            }
         } else {
            layer1->addBox(hier::Box(*fine_itr, hier::LocalId(ib), 0));
         }
      }

      hierarchy->makeNewPatchLevel(0, layer0);
      hierarchy->makeNewPatchLevel(1, layer1);

      // Create instance of hier::Variable database
      hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
      std::shared_ptr<hier::VariableContext> dummy(
         variable_db->getContext("dummy"));
      const hier::IntVector no_ghosts(dim, 0);

      // Make some dummy variables and data on the hierarchy
      std::shared_ptr<pdat::CellVariable<double> > cvar[NVARS];
      int cvindx[NVARS];
      cvar[0].reset(new pdat::CellVariable<double>(dim, "cvar0", 1));
      cvindx[0] = variable_db->registerVariableAndContext(
            cvar[0], dummy, no_ghosts);
      cvar[1].reset(new pdat::CellVariable<double>(dim, "cvar1", 1));
      cvindx[1] = variable_db->registerVariableAndContext(
            cvar[1], dummy, no_ghosts);
      cvar[2].reset(new pdat::CellVariable<double>(dim, "cvar2", 1));
      cvindx[2] = variable_db->registerVariableAndContext(
            cvar[2], dummy, no_ghosts);
      cvar[3].reset(new pdat::CellVariable<double>(dim, "cvar3", 1));
      cvindx[3] = variable_db->registerVariableAndContext(
            cvar[3], dummy, no_ghosts);

      std::shared_ptr<pdat::CellVariable<double> > cwgt(
         new pdat::CellVariable<double>(dim, "cwgt", 1));
      int cwgt_id = variable_db->registerVariableAndContext(
            cwgt, dummy, no_ghosts);

      // allocate data on hierarchy
      for (ln = 0; ln < 2; ++ln) {
         hierarchy->getPatchLevel(ln)->allocatePatchData(cwgt_id);
         for (iv = 0; iv < NVARS; ++iv) {
            hierarchy->getPatchLevel(ln)->allocatePatchData(cvindx[iv]);
         }
      }

      std::shared_ptr<math::HierarchyDataOpsReal<double> > cell_ops(
         new math::HierarchyCellDataOpsReal<double>(
            hierarchy,
            0,
            1));
      TBOX_ASSERT(cell_ops);

      std::shared_ptr<math::HierarchyDataOpsReal<double> > cwgt_ops(
         new math::HierarchyCellDataOpsReal<double>(
            hierarchy,
            0,
            1));

      std::shared_ptr<hier::Patch> patch;

      // Initialize control volume data for cell-centered components
      hier::Box coarse_fine = fine0 + fine1;
      coarse_fine.coarsen(ratio);
      for (ln = 0; ln < 2; ++ln) {
         std::shared_ptr<hier::PatchLevel> level(
            hierarchy->getPatchLevel(ln));
         for (hier::PatchLevel::iterator ip(level->begin());
              ip != level->end(); ++ip) {
            patch = *ip;
            std::shared_ptr<geom::CartesianPatchGeometry> pgeom(
               SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                  patch->getPatchGeometry()));
            TBOX_ASSERT(pgeom);
            const double* dx = pgeom->getDx();
            double cell_vol = dx[0];
            for (int i = 1; i < dim.getValue(); ++i) {
               cell_vol *= dx[i];
            }
            std::shared_ptr<pdat::CellData<double> > cvdata(
               SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                  patch->getPatchData(cwgt_id)));
            TBOX_ASSERT(cvdata);
            cvdata->fillAll(cell_vol);
            if (ln == 0) cvdata->fillAll(0.0, (coarse_fine * patch->getBox()));
         }
      }

#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      // Test #1: Print out control volume data and compute its integral

      // Test #1a: Check control volume data set properly
      // Expected: cwgt = 0.01 on coarse (except where finer patch exists) and
      // 0.0025 on fine level
      bool vol_test_passed = true;
      for (ln = 0; ln < 2; ++ln) {

         std::shared_ptr<hier::PatchLevel> level(
            hierarchy->getPatchLevel(ln));
         for (hier::PatchLevel::iterator ip(level->begin());
              ip != level->end(); ++ip) {
            patch = *ip;
            std::shared_ptr<pdat::CellData<double> > cvdata(
               SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                  patch->getPatchData(cwgt_id)));

            TBOX_ASSERT(cvdata);

            pdat::CellIterator cend(pdat::CellGeometry::end(cvdata->getBox()));
            for (pdat::CellIterator c(pdat::CellGeometry::begin(cvdata->getBox()));
                 c != cend && vol_test_passed; ++c) {
               pdat::CellIndex cell_index = *c;

               if (ln == 0) {
                  if ((coarse_fine * patch->getBox()).contains(cell_index)) {
                     if (!tbox::MathUtilities<double>::equalEps((*cvdata)(
                               cell_index), 0.0)) {
                        vol_test_passed = false;
                     }
                  } else {
                     if (!tbox::MathUtilities<double>::equalEps((*cvdata)(
                               cell_index), (dim == tbox::Dimension(2)) ? 0.01 : 0.001)) {
                        vol_test_passed = false;
                     }
                  }
               }

               if (ln == 1) {
                  if (!tbox::MathUtilities<double>::equalEps((*cvdata)(
                            cell_index), (dim == tbox::Dimension(2)) ? 0.0025 : 0.000125)) {
                     vol_test_passed = false;
                  }
               }
            }
         }
      }
      if (!vol_test_passed) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #1a: Check control volume data set properly"
         << std::endl;
         cwgt_ops->printData(cwgt_id, tbox::plog);
      }

      // Test #1b: math::HierarchyCellDataOpsReal::sumControlVolumes()
      // Expected: norm = 0.5
      double norm = cell_ops->sumControlVolumes(cvindx[0], cwgt_id);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(norm, 0.5)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #1b: math::HierarchyCellDataOpsReal::sumControlVolumes()\n"
         << "Expected value = 0.5 , Computed value = "
         << norm << std::endl;
      }

      // Test #2: math::HierarchyCellDataOpsReal::numberOfEntries()
      // Expected: num_data_points = 90 for 2D, 660 for 3D
      size_t num_data_points = cell_ops->numberOfEntries(cvindx[0]);
      if (num_data_points != ((dim == tbox::Dimension(2)) ? 90 : 660)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #2: math::HierarchyCellDataOpsReal::numberOfEntries()\n"
         << "Expected value = " << ((dim == tbox::Dimension(2)) ? 90 : 660)
         << " , Computed value = " << num_data_points << std::endl;
      }

      // Test #3a: math::HierarchyCellDataOpsReal::setToScalar()
      // Expected: v0 = 2.0
      double val0 = 2.0;
      cell_ops->setToScalar(cvindx[0], val0);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!doubleDataSameAsValue(cvindx[0], val0, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #3a: math::HierarchyCellDataOpsReal::setToScalar()\n"
         << "Expected: v0 = " << val0 << std::endl;
         cell_ops->printData(cvindx[0], tbox::plog);
      }

      // Test #3b: math::HierarchyCellDataOpsReal::setToScalar()
      // Expected: v1 = (4.0)
      cell_ops->setToScalar(cvindx[1], 4.0);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      double val1 = 4.0;
      if (!doubleDataSameAsValue(cvindx[1], val1, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #3b: math::HierarchyCellDataOpsReal::setToScalar()\n"
         << "Expected: v1 = " << val1 << std::endl;
         cell_ops->printData(cvindx[1], tbox::plog);
      }

      // Test #4: math::HierarchyCellDataOpsReal::copyData()
      // Expected: v2 = v1 = (4.0)
      cell_ops->copyData(cvindx[2], cvindx[1]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!doubleDataSameAsValue(cvindx[2], val1, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #4: math::HierarchyCellDataOpsReal::copyData()\n"
         << "Expected: v2 = " << val1 << std::endl;
         cell_ops->printData(cvindx[2], tbox::plog);
      }

      // Test #5: math::HierarchyCellDataOpsReal::swapData()
      // Expected: v0 = (4.0), v1 = (2.0)
      cell_ops->swapData(cvindx[0], cvindx[1]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!doubleDataSameAsValue(cvindx[0], val1, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #5a: math::HierarchyCellDataOpsReal::swapData()\n"
         << "Expected: v0 = " << val1 << std::endl;
         cell_ops->printData(cvindx[0], tbox::plog);
      }
      if (!doubleDataSameAsValue(cvindx[1], val0, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #5b: math::HierarchyCellDataOpsReal::swapData()\n"
         << "Expected: v1 = " << val0 << std::endl;
         cell_ops->printData(cvindx[1], tbox::plog);
      }

      // Test #6: math::HierarchyCellDataOpsReal::scale()
      // Expected: v2 = 0.25 * v2 = (1.0)
      cell_ops->scale(cvindx[2], 0.25, cvindx[2]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      double val_scale = 1.0;
      if (!doubleDataSameAsValue(cvindx[2], val_scale, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #6: math::HierarchyCellDataOpsReal::scale()\n"
         << "Expected: v2 = " << val_scale << std::endl;
         cell_ops->printData(cvindx[2], tbox::plog);
      }

      // Test #7: math::HierarchyCellDataOpsReal::add()
      // Expected: v3 = v0 + v1 = (6.0)
      cell_ops->add(cvindx[3], cvindx[0], cvindx[1]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      double val_add = 6.0;
      if (!doubleDataSameAsValue(cvindx[3], val_add, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #7: math::HierarchyCellDataOpsReal::add()\n"
         << "Expected: v3 = " << val_add << std::endl;
         cell_ops->printData(cvindx[3], tbox::plog);
      }

      // Reset v0: v0 = (0.0)
      cell_ops->setToScalar(cvindx[0], 0.0);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif

      // Test #8: math::HierarchyCellDataOpsReal::subtract()
      // Expected: v1 = v3 - v0 = (6.0)
      cell_ops->subtract(cvindx[1], cvindx[3], cvindx[0]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      double val_sub = 6.0;
      if (!doubleDataSameAsValue(cvindx[1], val_sub, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #8: math::HierarchyCellDataOpsReal::subtract()\n"
         << "Expected: v1 = " << val_sub << std::endl;
         cell_ops->printData(cvindx[1], tbox::plog);
      }

      // Test #9a: math::HierarchyCellDataOpsReal::addScalar()
      // Expected: v1 = v1 + (0.0) = (6.0)
      cell_ops->addScalar(cvindx[1], cvindx[1], 0.0);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      double val_addScalar = 6.0;
      if (!doubleDataSameAsValue(cvindx[1], val_addScalar, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #9a: math::HierarchyCellDataOpsReal::addScalar()\n"
         << "Expected: v1 = " << val_addScalar << std::endl;
         cell_ops->printData(cvindx[1], tbox::plog);
      }

      // Test #9b: math::HierarchyCellDataOpsReal::addScalar()
      // Expected: v2 = v2 + (0.0) = (1.0)
      cell_ops->addScalar(cvindx[2], cvindx[2], 0.0);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      val_addScalar = 1.0;
      if (!doubleDataSameAsValue(cvindx[2], val_addScalar, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #9b: math::HierarchyCellDataOpsReal::addScalar()\n"
         << "Expected: v2 = " << val_addScalar << std::endl;
         cell_ops->printData(cvindx[2], tbox::plog);
      }

      // Test #9c: math::HierarchyCellDataOpsReal::addScalar()
      // Expected: v2 = v2 + (3.0) = (4.0)
      cell_ops->addScalar(cvindx[2], cvindx[2], 3.0);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      val_addScalar = 4.0;
      if (!doubleDataSameAsValue(cvindx[2], val_addScalar, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #9c: math::HierarchyCellDataOpsReal::addScalar()\n"
         << "Expected: v2 = " << val_addScalar << std::endl;
         cell_ops->printData(cvindx[2], tbox::plog);
      }

      // Reset v3:  v3 = (0.5)
      cell_ops->setToScalar(cvindx[3], 0.5);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif

      // Test #10: math::HierarchyCellDataOpsReal::multiply()
      // Expected: v1 = v3 * v1 = (3.0)
      cell_ops->multiply(cvindx[1], cvindx[3], cvindx[1]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      double val_mult = 3.0;
      if (!doubleDataSameAsValue(cvindx[1], val_mult, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #10: math::HierarchyCellDataOpsReal::multiply()\n"
         << "Expected: v1 = " << val_mult << std::endl;
         cell_ops->printData(cvindx[1], tbox::plog);
      }

      // Test #11: math::HierarchyCellDataOpsReal::divide()
      // Expected: v0 = v2 / v1 = 1.3333333333
      cell_ops->divide(cvindx[0], cvindx[2], cvindx[1]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      double val_div = 1.33333333333;
      if (!doubleDataSameAsValue(cvindx[0], val_div, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #11: math::HierarchyCellDataOpsReal::divide()\n"
         << "Expected: v0 = " << val_div << std::endl;
         cell_ops->printData(cvindx[0], tbox::plog);
      }

      // Test #12: math::HierarchyCellDataOpsReal::reciprocal()
      // Expected:  v1 = 1 / v1 = (0.333333333)
      cell_ops->reciprocal(cvindx[1], cvindx[1]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      double val_rec = 0.33333333333;
      if (!doubleDataSameAsValue(cvindx[1], val_rec, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #12: math::HierarchyCellDataOpsReal::reciprocal()\n"
         << "Expected: v1 = " << val_rec << std::endl;
         cell_ops->printData(cvindx[1], tbox::plog);
      }

      // Test #13: math::HierarchyCellDataOpsReal::abs()
      // Expected:  v3 = abs(v2) = 4.0
      cell_ops->abs(cvindx[3], cvindx[2]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      double val_abs = 4.0;
      if (!doubleDataSameAsValue(cvindx[3], val_abs, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #13: math::HierarchyCellDataOpsReal::abs()\n"
         << "Expected: v3 = " << val_abs << std::endl;
         cell_ops->printData(cvindx[3], tbox::plog);
      }

      // Test #14: Place some bogus values on coarse level
      std::shared_ptr<pdat::CellData<double> > cdata;

      // set values
      std::shared_ptr<hier::PatchLevel> level_zero(
         hierarchy->getPatchLevel(0));
      for (hier::PatchLevel::iterator ip(level_zero->begin());
           ip != level_zero->end(); ++ip) {
         patch = *ip;
         cdata = SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>,
                            hier::PatchData>(patch->getPatchData(cvindx[2]));
         TBOX_ASSERT(cdata);
         hier::Index index0(dim, 2);
         hier::Index index1(dim, 3);
         index1(0) = 5;
         if (patch->getBox().contains(index0)) {
            (*cdata)(pdat::CellIndex(index0), 0) = 100.0;
         }
         if (patch->getBox().contains(index1)) {
            (*cdata)(pdat::CellIndex(index1), 0) = -1000.0;
         }
      }

      // check values
      bool bogus_value_test_passed = true;
      for (hier::PatchLevel::iterator ipp(level_zero->begin());
           ipp != level_zero->end(); ++ipp) {
         patch = *ipp;
         cdata = SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>,
                            hier::PatchData>(patch->getPatchData(cvindx[2]));
         TBOX_ASSERT(cdata);
         hier::Index index0(dim, 2);
         hier::Index index1(dim, 3);
         index1(0) = 5;

         pdat::CellIterator cend(pdat::CellGeometry::end(cdata->getBox()));
         for (pdat::CellIterator c(pdat::CellGeometry::begin(cdata->getBox()));
              c != cend && bogus_value_test_passed; ++c) {
            pdat::CellIndex cell_index = *c;

            if (cell_index == pdat::CellIndex(index0)) {
               if (!tbox::MathUtilities<double>::equalEps((*cdata)(cell_index),
                      100.0)) {
                  bogus_value_test_passed = false;
               }
            } else {
               if (cell_index == pdat::CellIndex(index1)) {
                  if (!tbox::MathUtilities<double>::equalEps((*cdata)(
                            cell_index), -1000.0)) {
                     bogus_value_test_passed = false;
                  }
               } else {
                  if (!tbox::MathUtilities<double>::equalEps((*cdata)(
                            cell_index), 4.0)) {
                     bogus_value_test_passed = false;
                  }
               }
            }
         }
      }
      if (!bogus_value_test_passed) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #14:  Place some bogus values on coarse level"
         << std::endl;
         cell_ops->printData(cvindx[2], tbox::plog);
      }

      // Test #15: math::HierarchyCellDataOpsReal::L1Norm() - w/o control weight
      // Expected:  bogus_l1_norm = 1452 in 2d, 3732 in 3d
      double bogus_l1_norm = cell_ops->L1Norm(cvindx[2]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(bogus_l1_norm,
             ((dim == tbox::Dimension(2)) ? 1452 : 3732))) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #15: math::HierarchyCellDataOpsReal::L1Norm()"
         << " - w/o control weight\n"
         << "Expected value = "
         << ((dim == tbox::Dimension(2)) ? 1452 : 3732) << ", Computed value = "
         << std::setprecision(12) << bogus_l1_norm << std::endl;
      }

      // Test #16: math::HierarchyCellDataOpsReal::L1Norm() - w/control weight
      // Expected:  correct_l1_norm = 2.0
      double correct_l1_norm = cell_ops->L1Norm(cvindx[2], cwgt_id);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(correct_l1_norm, 2.0)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #16: math::HierarchyCellDataOpsReal::L1Norm()"
         << " - w/control weight\n"
         << "Expected value = 2.0, Computed value = "
         << correct_l1_norm << std::endl;
      }

      // Test #17: math::HierarchyCellDataOpsReal::L2Norm()
      // Expected:  l2_norm = 2.82842712475
      double l2_norm = cell_ops->L2Norm(cvindx[2], cwgt_id);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(l2_norm, 2.82842712475)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #17: math::HierarchyCellDataOpsReal::L2Norm()\n"
         << "Expected value = 2.82842712475, Computed value = "
         << l2_norm << std::endl;
      }

      // Test #18: math::HierarchyCellDataOpsReal::L2Norm() - w/o control weight
      // Expected:  bogus_max_norm = 1000.0
      double bogus_max_norm = cell_ops->maxNorm(cvindx[2]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(bogus_max_norm, 1000.0)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #18: math::HierarchyCellDataOpsReal::L2Norm()"
         << " - w/o control weight\n"
         << "Expected value = 1000.0, Computed value = "
         << bogus_max_norm << std::endl;
      }

      // Test #19: math::HierarchyCellDataOpsReal::L2Norm() - w/control weight
      // Expected:  max_norm = 4.0
      double max_norm = cell_ops->maxNorm(cvindx[2], cwgt_id);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(max_norm, 4.0)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #19: math::HierarchyCellDataOpsReal::L2Norm()"
         << " - w/control weight\n"
         << "Expected value = 4.0, Computed value = "
         << max_norm << std::endl;
      }

      // Reset data and test sums, axpy's
      cell_ops->setToScalar(cvindx[0], 1.00);
      cell_ops->setToScalar(cvindx[1], 2.5);
      cell_ops->setToScalar(cvindx[2], 7.0);

#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      // Test #20: math::HierarchyCellDataOpsReal::linearSum()
      // Expected:  v3 = 5.0
      cell_ops->linearSum(cvindx[3], 2.0, cvindx[1], 0.00, cvindx[0]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      double val_linearSum = 5.0;
      if (!doubleDataSameAsValue(cvindx[3], val_linearSum, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #20: math::HierarchyCellDataOpsReal::linearSum()\n"
         << "Expected: v3 = " << val_linearSum << std::endl;
         cell_ops->printData(cvindx[3], tbox::plog);
      }

      // Test #21: math::HierarchyCellDataOpsReal::axmy()
      // Expected:  v3 = 6.5
      cell_ops->axmy(cvindx[3], 3.0, cvindx[1], cvindx[0]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      double val_axmy = 6.5;
      if (!doubleDataSameAsValue(cvindx[3], val_axmy, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #21: math::HierarchyCellDataOpsReal::axmy()\n"
         << "Expected: v3 = " << val_axmy << std::endl;
         cell_ops->printData(cvindx[3], tbox::plog);
      }

      // Test #22a: math::HierarchyCellDataOpsReal::dot() - (ind2) * (ind1)
      // Expected:  cdot = 8.75
      double cdot = cell_ops->dot(cvindx[2], cvindx[1], cwgt_id);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(cdot, 8.75)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #22a: math::HierarchyCellDataOpsReal::dot() - (ind2) * (ind1)\n"
         << "Expected Value = 8.75, Computed Value = "
         << cdot << std::endl;
      }

      // Test #22b: math::HierarchyCellDataOpsReal::dot() - (ind1) * (ind2)
      // Expected:  cdot = 8.75
      cdot = cell_ops->dot(cvindx[1], cvindx[2], cwgt_id);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(cdot, 8.75)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #22b: math::HierarchyCellDataOpsReal::dot() - (ind1) * (ind2)\n"
         << "Expected Value = 8.75, Computed Value = "
         << cdot << std::endl;
      }

      // deallocate data on hierarchy
      for (ln = 0; ln < 2; ++ln) {
         hierarchy->getPatchLevel(ln)->deallocatePatchData(cwgt_id);
         for (iv = 0; iv < NVARS; ++iv) {
            hierarchy->getPatchLevel(ln)->deallocatePatchData(cvindx[iv]);
         }
      }

      for (iv = 0; iv < NVARS; ++iv) {
         cvar[iv].reset();
      }
      cwgt.reset();

      geometry.reset();
      hierarchy.reset();
      cell_ops.reset();
      cwgt_ops.reset();

      if (num_failures == 0) {
         tbox::pout << "\nPASSED:  cell hiertest" << std::endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return num_failures;
}

/*
 * Returns true if all the data in the hierarchy is equal to the specified
 * value.  Returns false otherwise.
 */
static bool
doubleDataSameAsValue(
   int desc_id,
   double value,
   std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   bool test_passed = true;

   int ln;
   std::shared_ptr<hier::Patch> patch;
   for (ln = 0; ln < 2; ++ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         patch = *ip;
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
      }
   }

   return test_passed;
}

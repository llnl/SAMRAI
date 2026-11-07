/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Main program to test cell-centered complex patch data ops
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"

#include "SAMRAI/tbox/SAMRAIManager.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/math/HierarchyDataOpsComplex.h"
#include "SAMRAI/math/HierarchyCellDataOpsComplex.h"
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
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/VariableContext.h"

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <memory>

using namespace SAMRAI;

/* Helper function declarations */
static bool
complexDataSameAsValue(
   int desc_id,
   dcomplex value,
   std::shared_ptr<hier::PatchHierarchy> hierarchy);
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

   if (argc < 2) {
      TBOX_ERROR("Usage: " << argv[0] << " [dimension]");
   }

   const unsigned short d = static_cast<unsigned short>(atoi(argv[1]));
   TBOX_ASSERT(d > 0);
   TBOX_ASSERT(d <= SAMRAI::MAX_DIM_VAL);
   const tbox::Dimension dim(d);

   const std::string log_fn = std::string("cell_cplxtest.")
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

      const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
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

      // Make some dummy variables and register them with hier::VariableDatabase
      std::shared_ptr<pdat::CellVariable<dcomplex> > cvar[NVARS];
      int cvindx[NVARS];
      cvar[0].reset(new pdat::CellVariable<dcomplex>(dim, "cvar0", 1));
      cvindx[0] = variable_db->registerVariableAndContext(
            cvar[0], dummy, no_ghosts);
      cvar[1].reset(new pdat::CellVariable<dcomplex>(dim, "cvar1", 1));
      cvindx[1] = variable_db->registerVariableAndContext(
            cvar[1], dummy, no_ghosts);
      cvar[2].reset(new pdat::CellVariable<dcomplex>(dim, "cvar2", 1));
      cvindx[2] = variable_db->registerVariableAndContext(
            cvar[2], dummy, no_ghosts);
      cvar[3].reset(new pdat::CellVariable<dcomplex>(dim, "cvar3", 1));
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

      std::shared_ptr<math::HierarchyDataOpsComplex> cell_ops(
         new math::HierarchyCellDataOpsComplex(
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
#if defined(HAVE_RAJA)
            tbox::parallel_synchronize();
#endif
         }
      }

      // Test #1: Print out control volume data and compute its integral

      // Test #1a: Check control volume data set properly
      // Expected: cwgt = 0.01 on coarse (except where finer patch exists) and
      // 0.0025 on fine level for 2d.  0.001 and 0.000125 for 3d.
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
                     double compare;
                     if ((dim == tbox::Dimension(2)))
                        compare = 0.01;
                     else
                        compare = 0.001;

                     if (!tbox::MathUtilities<double>::equalEps((*cvdata)(
                               cell_index), compare)) {

                        vol_test_passed = false;
                     }
                  }
               }

               if (ln == 1) {
                  double compare;
                  if ((dim == tbox::Dimension(2)))
                     compare = 0.0025;
                  else
                     compare = 0.000125;

                  if (!tbox::MathUtilities<double>::equalEps((*cvdata)(
                            cell_index), compare)) {
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

      // Test #1b: math::HierarchyCellDataOpsComplex::sumControlVolumes()
      // Expected: norm = 0.5
      double norm = cell_ops->sumControlVolumes(cvindx[0], cwgt_id);
#if defined(HAVE_RAJA)
            tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(norm, 0.5)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #1b: math::HierarchyCellDataOpsComplex::sumControlVolumes()\n"
         << "Expected value = 0.5 , Computed value = "
         << norm << std::endl;
      }

      // Test #2: math::HierarchyCellDataOpsComplex::numberOfEntries()
      // Expected: num_data_points = 90 in 2d, 660 in 3d
      size_t num_data_points = cell_ops->numberOfEntries(cvindx[0]);

      {
         size_t compare;
         if ((dim == tbox::Dimension(2)))
            compare = 90;
         else
            compare = 660;

         if (num_data_points != compare) {
            ++num_failures;
            tbox::perr
            << "FAILED: - Test #2: math::HierarchyCellDataOpsReal::numberOfEntries()\n"
            << "Expected value = " << compare
            << ", Computed value = "
            << num_data_points << std::endl;
         }
      }

      // Test #3a: math::HierarchyCellDataOpsComplex::setToScalar()
      // Expected: v0 = (2.0,1.5)
      dcomplex val0 = dcomplex(2.0, 1.5);
      cell_ops->setToScalar(cvindx[0], val0);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!complexDataSameAsValue(cvindx[0], val0, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #3a: math::HierarchyCellDataOpsComplex::setToScalar()\n"
         << "Expected: v0 = " << val0 << std::endl;
         cell_ops->printData(cvindx[0], tbox::plog);
      }

      // Test #3b: math::HierarchyCellDataOpsComplex::setToScalar()
      // Expected: v1 = (4.0,3.0)
      dcomplex val1(4.0, 3.0);
      cell_ops->setToScalar(cvindx[1], val1);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!complexDataSameAsValue(cvindx[1], val1, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #3b: math::HierarchyCellDataOpsComplex::setToScalar()\n"
         << "Expected: v1 = " << val1 << std::endl;
         cell_ops->printData(cvindx[1], tbox::plog);
      }

      // Test #4: math::HierarchyCellDataOpsComplex::copyData()
      // Expected: v2 = v1 = (4.0, 3.0)
      cell_ops->copyData(cvindx[2], cvindx[1]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!complexDataSameAsValue(cvindx[2], val1, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #4: math::HierarchyCellDataOpsComplex::copyData()\n"
         << "Expected: v2 = v1 = " << val1 << std::endl;
         cell_ops->printData(cvindx[2], tbox::plog);
      }

      // Test #5: math::HierarchyCellDataOpsComplex::swapData()
      // Expected: v0 = (4.0, 3.0), v1 = (2.0,1.5)
      cell_ops->swapData(cvindx[0], cvindx[1]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!complexDataSameAsValue(cvindx[0], val1, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #5a: math::HierarchyCellDataOpsComplex::swapData()\n"
         << "Expected: v0 = " << val1 << std::endl;
         cell_ops->printData(cvindx[0], tbox::plog);
      }
      if (!complexDataSameAsValue(cvindx[1], val0, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #5b: math::HierarchyCellDataOpsComplex::swapData()\n"
         << "Expected: v1 = " << val0 << std::endl;
         cell_ops->printData(cvindx[1], tbox::plog);
      }

      // Test #6: math::HierarchyCellDataOpsComplex::scale()
      // Expected: v2 = 0.25 * v2 = (1.0,0.75)
      cell_ops->scale(cvindx[2], 0.25, cvindx[2]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex val_scale(1.0, 0.75);
      if (!complexDataSameAsValue(cvindx[2], val_scale, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #6: math::HierarchyCellDataOpsComplex::scale()\n"
         << "Expected: v2 = " << val_scale << std::endl;
         cell_ops->printData(cvindx[2], tbox::plog);
      }

      // Test #7: math::HierarchyCellDataOpsComplex::add()
      // Expected: v3 = v0 + v1 = (6.0, 4.5)
      cell_ops->add(cvindx[3], cvindx[0], cvindx[1]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex val_add(6.0, 4.5);
      if (!complexDataSameAsValue(cvindx[3], val_add, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #7: math::HierarchyCellDataOpsComplex::add()\n"
         << "Expected: v3 = " << val_add << std::endl;
         cell_ops->printData(cvindx[3], tbox::plog);
      }

      // Reset v0:  v0 = (0.0,4.5)
      cell_ops->setToScalar(cvindx[0], dcomplex(0.0, 4.5));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif

      // Test #8: math::HierarchyCellDataOpsComplex::subtract()
      // Expected: v1 = v3 - v0 = (6.0,0.0)
      cell_ops->subtract(cvindx[1], cvindx[3], cvindx[0]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex val_sub(6.0, 0.0);
      if (!complexDataSameAsValue(cvindx[1], val_sub, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #8: math::HierarchyCellDataOpsComplex::subtract()\n"
         << "Expected: v1 = " << val_sub << std::endl;
         cell_ops->printData(cvindx[1], tbox::plog);
      }

      // Test #9a: math::HierarchyCellDataOpsComplex::addScalar()
      // Expected: v1 = v1 + (0.0,-4.0) = (6.0,-4.0)
      cell_ops->addScalar(cvindx[1], cvindx[1], dcomplex(0.0, -4.0));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex val_addScalar(6.0, -4.0);
      if (!complexDataSameAsValue(cvindx[1], val_addScalar, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #9a: math::HierarchyCellDataOpsComplex::addScalar()\n"
         << "Expected: v1 = " << val_addScalar << std::endl;
         cell_ops->printData(cvindx[1], tbox::plog);
      }

      // Test #9b: math::HierarchyCellDataOpsComplex::addScalar()
      // Expected:   v2 = v2 + (0.0,0.25) = (1.0,1.0)
      cell_ops->addScalar(cvindx[2], cvindx[2], dcomplex(0.0, 0.25));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      val_addScalar = dcomplex(1.0, 1.0);
      if (!complexDataSameAsValue(cvindx[2], val_addScalar, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #9c: math::HierarchyCellDataOpsComplex::addScalar()\n"
         << "Expected: v2 = " << val_addScalar << std::endl;
         cell_ops->printData(cvindx[2], tbox::plog);
      }

      // Test #9c: math::HierarchyCellDataOpsComplex::addScalar()
      // Expected:  v2 = v2 + (3.0,-4.0) = (4.0,-3.0)
      cell_ops->addScalar(cvindx[2], cvindx[2], dcomplex(3.0, -4.0));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      val_addScalar = dcomplex(4.0, -3.0);
      if (!complexDataSameAsValue(cvindx[2], val_addScalar, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #9d: math::HierarchyCellDataOpsComplex::addScalar()\n"
         << "Expected: v2 = " << val_addScalar << std::endl;
         cell_ops->printData(cvindx[2], tbox::plog);
      }

      // Reset v3: v3 = (0.5, 0.0)
      cell_ops->setToScalar(cvindx[3], dcomplex(0.5, 0.0));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif

      // Test #10: math::HierarchyCellDataOpsComplex::multiply()
      // Expected:  v1 = v3 * v1 = (3.0,-2.0)
      cell_ops->multiply(cvindx[1], cvindx[3], cvindx[1]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex val_mult(3.0, -2.0);
      if (!complexDataSameAsValue(cvindx[1], val_mult, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #10: math::HierarchyCellDataOpsComplex::multiply()\n"
         << "Expected: v1 = " << val_mult << std::endl;
         cell_ops->printData(cvindx[1], tbox::plog);
      }

      // Test #11: math::HierarchyCellDataOpsComplex::divide()
      // Expected:  v0 = v2 / v1 = (1.3846153846154,-0.076923076923077)
      cell_ops->divide(cvindx[0], cvindx[2], cvindx[1]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex val_div(1.3846153846154, -0.076923076923077);
      if (!complexDataSameAsValue(cvindx[0], val_div, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #11: math::HierarchyCellDataOpsComplex::divide()\n"
         << "Expected: v0 = " << val_div << std::endl;
         cell_ops->printData(cvindx[0], tbox::plog);
      }

      // Test #12: math::HierarchyCellDataOpsComplex::reciprocal()
      // Expected:  v1 = 1 / v1 = (0.23076923076923, 0.15384615384615)
      cell_ops->reciprocal(cvindx[1], cvindx[1]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex val_rec(0.23076923076923, 0.15384615384615);
      if (!complexDataSameAsValue(cvindx[1], val_rec, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #12: math::HierarchyCellDataOpsComplex::reciprocal()\n"
         << "Expected: v1 = " << val_rec << std::endl;
         cell_ops->printData(cvindx[1], tbox::plog);
      }

      // Test #13:  Place some bogus values on coarse level
      std::shared_ptr<pdat::CellData<dcomplex> > cdata;

      // set values
      std::shared_ptr<hier::PatchLevel> level_zero(
         hierarchy->getPatchLevel(0));
      for (hier::PatchLevel::iterator ip(level_zero->begin());
           ip != level_zero->end(); ++ip) {
         patch = *ip;
         cdata = SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>,
                            hier::PatchData>(patch->getPatchData(cvindx[2]));
         TBOX_ASSERT(cdata);
         hier::Index index0(dim, 2);
         hier::Index index1(dim, 3);
         index1(0) = 5;
         if (patch->getBox().contains(index0)) {
            (*cdata)(pdat::CellIndex(index0), 0) = dcomplex(100.0, -50.0);
         }
         if (patch->getBox().contains(index1)) {
            (*cdata)(pdat::CellIndex(index1), 0) = dcomplex(-1000.0, 20.0);
         }
      }

      // check values
      bool bogus_value_test_passed = true;
      for (hier::PatchLevel::iterator ipp(level_zero->begin());
           ipp != level_zero->end(); ++ipp) {
         patch = *ipp;
         cdata = SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>,
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
               if (!tbox::MathUtilities<dcomplex>::equalEps((*cdata)(cell_index),
                      dcomplex(100.0, -50.0))) {
                  bogus_value_test_passed = false;
               }
            } else {
               if (cell_index == pdat::CellIndex(index1)) {
                  if (!tbox::MathUtilities<dcomplex>::equalEps((*cdata)(
                            cell_index),
                         dcomplex(-1000.0, 20.0))) {
                     bogus_value_test_passed = false;
                  }
               } else {
                  if (!tbox::MathUtilities<dcomplex>::equalEps((*cdata)(
                            cell_index),
                         dcomplex(4.0, -3.0))) {
                     bogus_value_test_passed = false;
                  }
               }
            }
         }
      }
      if (!bogus_value_test_passed) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #13:  Place some bogus values on coarse level"
         << std::endl;
         cell_ops->printData(cvindx[2], tbox::plog);
      }

      // Test norms on patch data with cvindx[2] on hierarchy with bogus values

      // Test #14: math::HierarchyCellDataOpsComplex::L1Norm() - w/o control weight
      // Expected:  bogus_l1_norm = 1552.00337888 in 2d, 4402.00337888 in 3d.
      double bogus_l1_norm = cell_ops->L1Norm(cvindx[2]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif

      {
         double compare;
         if ((dim == tbox::Dimension(2)))
            compare = 1552.00337888;
         else
            compare = 4402.00337888;

         if (!tbox::MathUtilities<double>::equalEps(bogus_l1_norm, compare)) {
            ++num_failures;
            tbox::perr
            << "FAILED: - Test #14: math::HierarchyCellDataOpsComplex::L1Norm()"
            << " - w/o control weight\n"
            << "Expected value = " << compare
            << ", Computed value = "
            << std::setprecision(12) << bogus_l1_norm << std::endl;
         }
      }

      // Test #15: math::HierarchyCellDataOpsComplex::L1Norm() - w/control weight
      // Expected:  correct_l1_norm = 2.5
      double correct_l1_norm = cell_ops->L1Norm(cvindx[2], cwgt_id);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(correct_l1_norm, 2.5)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #15: math::HierarchyCellDataOpsComplex::L1Norm()"
         << " - w/control weight\n"
         << "Expected value = 2.5, Computed value = "
         << correct_l1_norm << std::endl;
      }

      // Test #16: math::HierarchyCellDataOpsComplex::L2Norm()
      // Expected:  l2_norm = 3.53553390593
      double l2_norm = cell_ops->L2Norm(cvindx[2], cwgt_id);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(l2_norm, 3.53553390593)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #16: math::HierarchyCellDataOpsComplex::L2Norm()\n"
         << "Expected value = 3.53553390593, Computed value = "
         << l2_norm << std::endl;
      }

      // Test #17: math::HierarchyCellDataOpsComplex::maxNorm() - w/o control weight
      // Expected:  bogus_max_norm = 1000.19998
      double bogus_max_norm = cell_ops->maxNorm(cvindx[2]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(bogus_max_norm, 1000.19998)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #17: math::HierarchyCellDataOpsComplex::maxNorm() "
         << "- w/o control weight\n"
         << "Expected value = 1000.19998, Computed value = "
         << bogus_max_norm << std::endl;
      }

      // Test #18: math::HierarchyCellDataOpsComplex::maxNorm() - w/control weight
      // Expected:  max_norm = 5.0
      double max_norm = cell_ops->maxNorm(cvindx[2], cwgt_id);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!tbox::MathUtilities<double>::equalEps(max_norm, 5.0)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #18: math::HierarchyCellDataOpsComplex::maxNorm() "
         << "- w/control weight\n"
         << "Expected value = 5.0, Computed value = "
         << max_norm << std::endl;
      }

      // Reset data and test sums, axpy's
      cell_ops->setToScalar(cvindx[0], dcomplex(1.0, -3.0));
      cell_ops->setToScalar(cvindx[1], dcomplex(2.5, 3.0));
      cell_ops->setToScalar(cvindx[2], dcomplex(7.0, 0.0));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif

      // Test #19: math::HierarchyCellDataOpsComplex::linearSum()
      // Expected:  v3 = (2.0,5.0)
      cell_ops->linearSum(cvindx[3],
         dcomplex(2.0, 0.0), cvindx[1], dcomplex(0.0, -1.0), cvindx[0]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex val_linearSum(2.0, 5.0);
      if (!complexDataSameAsValue(cvindx[3], val_linearSum, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #19: math::HierarchyCellDataOpsComplex::linearSum()\n"
         << "Expected: v3 = " << val_linearSum << std::endl;
         cell_ops->printData(cvindx[3], tbox::plog);
      }

      // Test #20: math::HierarchyCellDataOpsComplex::axmy()
      // Expected:  v3 = (6.5,12.0)
      cell_ops->axmy(cvindx[3], 3.0, cvindx[1], cvindx[0]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex val_axmy(6.5, 12.0);
      if (!complexDataSameAsValue(cvindx[3], val_axmy, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #20: math::HierarchyCellDataOpsComplex::axmy()\n"
         << "Expected: v3 = " << val_axmy << std::endl;
         cell_ops->printData(cvindx[3], tbox::plog);
      }

      // Test #21a: math::HierarchyCellDataOpsComplex::dot()
      // Expected:  cdot = (8.75,-10.5)
      dcomplex cdot = cell_ops->dot(cvindx[2], cvindx[1], cwgt_id);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex ans_2_dot_1(8.75, -10.5);
      if (!tbox::MathUtilities<dcomplex>::equalEps(cdot, ans_2_dot_1)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #21a: math::HierarchyCellDataOpsComplex::dot()\n"
         << "Expected value = (8.75,-10.5), Computed value = "
         << cdot << std::endl;
      }

      // Test #21b: math::HierarchyCellDataOpsComplex::dot()
      // Expected:  cdot = (8.75,10.5)
      dcomplex cdot2 = cell_ops->dot(cvindx[1], cvindx[2], cwgt_id);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      dcomplex ans_1_dot_2(8.75, 10.5);
      if (!tbox::MathUtilities<dcomplex>::equalEps(cdot2, ans_1_dot_2)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #21b: math::HierarchyCellDataOpsComplex::dot()\n"
         << "Expected value = (8.75,10.5), Computed value = "
         << cdot2 << std::endl;
      }

      // Test #22: math::HierarchyCellDataOpsComplex::abs()
      // Expected:  abs(v0) = 5.0
      cell_ops->setToScalar(cvindx[0], dcomplex(4.0, -3.0));
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      cell_ops->abs(cwgt_id, cvindx[0]);
#if defined(HAVE_RAJA)
      tbox::parallel_synchronize();
#endif
      if (!doubleDataSameAsValue(cwgt_id, 5.0, hierarchy)) {
         ++num_failures;
         tbox::perr
         << "FAILED: - Test #22: math::HierarchyCellDataOpsComplex::abs()\n"
         << "Expected: abs(v0) = 5.0" << std::endl;
         cwgt_ops->printData(cwgt_id, tbox::plog);
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
         tbox::pout << "\nPASSED:  cell cplxtest" << std::endl;
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
complexDataSameAsValue(
   int desc_id,
   dcomplex value,
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
         std::shared_ptr<pdat::CellData<dcomplex> > cvdata(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<dcomplex>, hier::PatchData>(
               patch->getPatchData(desc_id)));

         TBOX_ASSERT(cvdata);

         pdat::CellIterator cend(pdat::CellGeometry::end(cvdata->getBox()));
         for (pdat::CellIterator c(pdat::CellGeometry::begin(cvdata->getBox()));
              c != cend && test_passed; ++c) {
            pdat::CellIndex cell_index = *c;
            if (!tbox::MathUtilities<dcomplex>::equalEps((*cvdata)(cell_index),
                   value)) {
               test_passed = false;
            }
         }
      }
   }

   return test_passed;
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

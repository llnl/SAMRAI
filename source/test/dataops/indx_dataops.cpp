/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Main program to test index data operations
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

// class holding information stored in index data
#include "SampleIndexData.h"

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/IndexData.h"
#include "SAMRAI/pdat/IndexVariable.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/IOStream.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/VariableContext.h"

#include <memory>

using namespace SAMRAI;

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

   const std::string log_fn = std::string("indx_dataops.")
      + tbox::Utilities::intToString(dim.getValue(), 1) + "d.log";
   tbox::PIO::logAllNodes(log_fn);

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

/*
 ************************************************************************
 *
 *   Create a simple 2-level hierarchy to test.
 *   (NOTE: it is setup to work on at most 2 processors)
 *
 ************************************************************************
 */
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
      hier::BoxContainer fine_domain;
      coarse_domain.pushBack(coarse0);
      coarse_domain.pushBack(coarse1);
      fine_domain.pushBack(fine0);
      fine_domain.pushBack(fine1);

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
      const int n_fine_boxes = fine_domain.size();

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

      hier::BoxContainer::iterator fine_itr = fine_domain.begin();
      for (int ib = 0; ib < n_fine_boxes; ++ib) {
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

      /*
       * Create an IndexData<SampleIndexData> variable and register it with
       * the variable database.
       */
      hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
      std::shared_ptr<hier::VariableContext> cxt(
         variable_db->getContext("dummy"));
      const hier::IntVector no_ghosts(dim, 0);

      std::shared_ptr<pdat::IndexVariable<SampleIndexData,
                                            pdat::CellGeometry> > data(
         new pdat::IndexVariable<SampleIndexData, pdat::CellGeometry>(
            dim, "sample"));
      int data_id = variable_db->registerVariableAndContext(
            data, cxt, no_ghosts);

/*
 ************************************************************************
 *
 *   Set index data.
 *
 ************************************************************************
 */

      /*
       * Loop over hierarchy levels and set index data on cells of patches
       */
      int counter = 0;
      std::ostream& os = tbox::plog;
      for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
         std::shared_ptr<hier::PatchLevel> level(
            hierarchy->getPatchLevel(ln));

         // allocate "sample" data
         level->allocatePatchData(data_id);
         os << "\nLevel: " << level->getLevelNumber() << " ";

         // loop over patches on level
         for (hier::PatchLevel::iterator ip(level->begin());
              ip != level->end(); ++ip) {
            std::shared_ptr<hier::Patch> patch(*ip);
            os << "Patch: " << patch->getLocalId() << std::endl;

            // access sample data from patch
            std::shared_ptr<pdat::IndexData<SampleIndexData,
                                              pdat::CellGeometry> > sample(
               SAMRAI_SHARED_PTR_CAST<pdat::IndexData<SampleIndexData, pdat::CellGeometry>,
                          hier::PatchData>(
                  patch->getPatchData(data_id)));
            TBOX_ASSERT(sample);

            // iterate over cells of patch and invoke one "SampleIndexData"
            // instance on each cell (its possible to do more).
            pdat::CellIterator icend(pdat::CellGeometry::end(patch->getBox()));
            for (pdat::CellIterator ic(pdat::CellGeometry::begin(patch->getBox()));
                 ic != icend; ++ic) {
               SampleIndexData sd;
               sd.setInt(counter);
               sample->appendItem(*ic, sd);
               ++counter;
            }

            // iterate over the "SampleIndexData" index data stored on the patch
            // and dump the integer stored on it.
            int currData = counter - 1;
            pdat::IndexData<SampleIndexData, pdat::CellGeometry>::iterator idend(*sample, false);
            for (pdat::IndexData<SampleIndexData,
                                 pdat::CellGeometry>::iterator id(*sample, true);
                 id != idend;
                 ++id) {
               os << "      SampleIndexData data: " << id->getInt()
                  << std::endl;
               if (id->getInt() != currData) {
                  ++num_failures;
                  tbox::perr
                  << "FAILED: - Index data set incorrectly" << std::endl;
               }
               --currData;
            }

         }
      }

      geometry.reset();
      hierarchy.reset();

      if (num_failures == 0) {
         tbox::pout << "\nPASSED:  indx dataops" << std::endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return num_failures;
}

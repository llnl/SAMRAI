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

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/Index.h"

#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/hier/PatchDataRestartManager.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/tbox/RestartManager.h"

#include "SAMRAI/pdat/SparseData.h"
#include "SAMRAI/pdat/SparseDataVariable.h"

#include <memory>

using namespace SAMRAI;

bool
checkIterators(
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   const int data_id1);

bool
checkCopyOps(
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   const int data_id1,
   const int data_id2);

bool
checkRemoveOps(
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   const int data_id1);

std::shared_ptr<geom::CartesianGridGeometry>
getGeometry(
   hier::BoxContainer& coarse_domain,
   hier::BoxContainer& fine_domain,
   const tbox::Dimension& dim);

void
_getDblKeys(
   std::vector<std::string>& keys);
void
_getDblValues(
   double* values,
   int mult = 1);
void
_getIntKeys(
   std::vector<std::string>& keys);
void
_getIntValues(
   int* values,
   int mult = 1);

const int DSIZE = 10;
const int ISIZE = 4;

int main(
   int argc,
   char* argv[]) {

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   // Note: For these simple tests we allow at most 2 processors.
   tbox::SAMRAI_MPI mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   const int nproc = mpi.getSize();

   // Currently this test only works for 2 dimensions.
   const tbox::Dimension dim(2);

   const std::string log_fn = std::string("sparse_dataops.")
      + tbox::Utilities::intToString(dim.getValue(), 1) + "d.log";
   tbox::PIO::logAllNodes(log_fn);

   bool check_it = true;
   bool copy_ops = true;
   bool remove_ops = true;

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {
      std::ostream& os = tbox::plog;

      /********************************************************************
      *
      *   Create a simple 2-level hierarchy to test.
      *   (NOTE: it is setup to work on at most 2 processors)
      *
      ********************************************************************/
      hier::BoxContainer coarse_domain;
      hier::BoxContainer fine_domain;

      std::shared_ptr<geom::CartesianGridGeometry> geometry(
         getGeometry(coarse_domain, fine_domain, dim));

      std::shared_ptr<hier::PatchHierarchy> hierarchy(
         new hier::PatchHierarchy("PatchHierarchy", geometry));

      hier::IntVector ratio(dim, 2);
      hierarchy->setMaxNumberOfLevels(2);
      hierarchy->setRatioToCoarserLevel(ratio, 1);

      const int n_coarse_boxes = coarse_domain.size();
      const int n_fine_boxes = fine_domain.size();

      std::shared_ptr<hier::BoxLevel> layer0(
         std::make_shared<hier::BoxLevel>(
            hier::IntVector(dim, 1), geometry));
      std::shared_ptr<hier::BoxLevel> layer1(
         std::make_shared<hier::BoxLevel>(ratio, geometry));

      hier::BoxContainer::iterator coarse_domain_itr = coarse_domain.begin();
      for (int ib = 0; ib < n_coarse_boxes; ++ib, ++coarse_domain_itr) {
         if (nproc > 1) {
            if (ib == layer0->getMPI().getRank()) {
               layer0->addBox(hier::Box(*coarse_domain_itr,
                     hier::LocalId(ib), layer0->getMPI().getRank()));
            }
         } else {
            layer0->addBox(hier::Box(*coarse_domain_itr,
                  hier::LocalId(ib), 0));
         }
      }

      hier::BoxContainer::iterator fine_domain_itr = fine_domain.begin();
      for (int ib = 0; ib < n_fine_boxes; ++ib, ++fine_domain_itr) {
         if (nproc > 1) {
            if (ib == layer1->getMPI().getRank()) {
               layer1->addBox(hier::Box(*fine_domain_itr,
                     hier::LocalId(ib), layer1->getMPI().getRank()));
            }
         } else {
            layer1->addBox(hier::Box(*fine_domain_itr,
                  hier::LocalId(ib), 0));
         }
      }

      hierarchy->makeNewPatchLevel(0, layer0);
      hierarchy->makeNewPatchLevel(1, layer1);

      /*
       * Create an SparseData<BOX_GEOMETRY> variable and register it with
       * the variable database.
       */
      hier::VariableDatabase* variable_db =
         hier::VariableDatabase::getDatabase();
      hier::PatchDataRestartManager* pdrm =
         hier::PatchDataRestartManager::getManager();
      std::shared_ptr<hier::VariableContext> cxt(
         variable_db->getContext("dummy"));
      const hier::IntVector no_ghosts(dim, 0);

      typedef pdat::SparseData<pdat::CellGeometry> LSparseData;
      typedef pdat::SparseDataVariable<pdat::CellGeometry> LSparseDataVar;

      std::vector<std::string> dkeys;
      _getDblKeys(dkeys);
      std::vector<std::string> ikeys;
      _getIntKeys(ikeys);

      std::shared_ptr<LSparseDataVar> data1(
         new LSparseDataVar(dim, "sample1", dkeys, ikeys));
      int data_id1 = variable_db->registerVariableAndContext(
            data1, cxt, no_ghosts);

      std::shared_ptr<LSparseDataVar> data2(
         new LSparseDataVar(dim, "sample2", dkeys, ikeys));
      int data_id2 = variable_db->registerVariableAndContext(
            data2, cxt, no_ghosts);

      // set us up for restart.
      pdrm->registerPatchDataForRestart(data_id1);
      pdrm->registerPatchDataForRestart(data_id2);

      for (int i = 0; i < 2; ++i) {
         // allocate "sample" data
         hierarchy->getPatchLevel(i)->allocatePatchData(data_id1);
         hierarchy->getPatchLevel(i)->allocatePatchData(data_id2);
      }

      /*
       * Loop over hierarchy levels and populate data.
       */
      for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
         std::shared_ptr<hier::PatchLevel> level(
            hierarchy->getPatchLevel(ln));

         // loop over patches on level
         for (hier::PatchLevel::iterator ip(level->begin());
              ip != level->end(); ++ip) {
            const std::shared_ptr<hier::Patch>& patch = *ip;

            // access sample data from patch
            std::shared_ptr<LSparseData> sample1(
               SAMRAI_SHARED_PTR_CAST<LSparseData, hier::PatchData>(
                  patch->getPatchData(data_id1)));
            std::shared_ptr<LSparseData> sample2(
               SAMRAI_SHARED_PTR_CAST<LSparseData, hier::PatchData>(
                  patch->getPatchData(data_id2)));
            TBOX_ASSERT(sample1);
            TBOX_ASSERT(sample2);

            // add items to the sparse data objects.
            pdat::CellIterator ic(pdat::CellGeometry::begin(patch->getBox()));
            pdat::CellIterator icend(pdat::CellGeometry::end(patch->getBox()));
            for ( ; ic != icend; ++ic) {
               const hier::Index* idx = &(*ic);
               LSparseData::iterator iter1 = sample1->registerIndex(*idx);
               LSparseData::iterator iter2 = sample2->registerIndex(*idx);

               double dvals1[DSIZE], dvals2[DSIZE];
               for (int i = 0; i < DSIZE; ++i) {
                  dvals1[i] = i;
                  dvals2[i] = i * 2;
               }
               int ivals1[ISIZE], ivals2[ISIZE];
               for (int j = 0; j < ISIZE; ++j) {
                  ivals1[j] = j;
                  ivals2[j] = j * 2;
               }
               iter1.insert(dvals1, ivals1);
               iter2.insert(dvals2, ivals2);
            }

            LSparseData::iterator iter1(sample1.get());
            LSparseData::iterator iter2(sample2.get());

            for ( ; iter1 != sample1->end() && iter2 != sample2->end();
                  ++iter1, ++iter2) {
               os << iter1;
               os << iter2;
            }
         }
      }

      /********************************************************************
      *
      *   Run the tests
      *
      ********************************************************************/
      check_it = checkIterators(hierarchy, data_id1);
      // Test 1 check iterators
      os << (check_it ? "PASSED: " : "FAILED: ")
         << "Test 1: Iterator test"
         << std::endl;

      // Test 2: copying items
      copy_ops = checkCopyOps(hierarchy, data_id1, data_id2);
      os << (copy_ops ? "PASSED: " : "FAILED: ")
         << "Test 2: Copy test"
         << std::endl;

      // Test 3:  removing and erasing items.
      remove_ops = checkRemoveOps(hierarchy, data_id1);
      os << (remove_ops ? "PASSED: " : "FAILED: ")
         << "Test 3: Remove test"
         << std::endl;

      /*
       * Tests Completed.
       */
      geometry.reset();
      hierarchy.reset();
   }

   if (check_it && copy_ops && remove_ops) {
      tbox::pout << "\nPASSED: sparse data ops" << std::endl;
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return 0;
}

void
_getDblKeys(std::vector<std::string>& keys) {
   for (int i = 0; i < DSIZE; ++i) {
      std::stringstream key_name;
      key_name << "key_" << i;
      keys.push_back(key_name.str());
   }
}

void
_getDblValues(double* values, int mult)
{
   for (int i = 0; i < DSIZE; ++i) {
      values[i] = (double)(i * mult);
   }
}

void
_getIntKeys(std::vector<std::string>& keys) {
   for (int i = 0; i < ISIZE; ++i) {
      std::stringstream key_name;
      key_name << "key_" << i;
      keys.push_back(key_name.str());
   }
}

void
_getIntValues(int* values, int mult)
{
   for (int i = 0; i < ISIZE; ++i) {
      values[i] = i * mult;
   }
}

bool
checkIterators(
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   const int data_id1)
{

   // Test 1 - check the functionality of the SparseData API
   //
   int num_failures(0);
   typedef pdat::SparseData<pdat::CellGeometry> LSparseData;
   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& patch = *ip;

         std::shared_ptr<LSparseData> sample(
            SAMRAI_SHARED_PTR_CAST<LSparseData, hier::PatchData>(
               patch->getPatchData(data_id1)));
         TBOX_ASSERT(sample);

         // Test #1a: check empty.  This should be false.
         if (sample->empty()) {
            ++num_failures;
            tbox::perr
            << "FAILED: - sparse data structure reports empty. "
            << std::endl;
         }

         pdat::CellIterator ic(pdat::CellGeometry::begin(patch->getBox()));
         pdat::CellIterator icend(pdat::CellGeometry::end(patch->getBox()));
         for ( ; ic != icend; ++ic) {
            const hier::Index& idx = *ic;
            LSparseData::AttributeIterator it = sample->begin(idx);
            for ( ; it != sample->end(idx); ++it) {

               // check element access.
               for (int i = 0; i < DSIZE; ++i) {
                  if (it[pdat::DoubleAttributeId(i)] != i) {
                     ++num_failures;
                  }
               }

               for (int j = 0; j < ISIZE; ++j) {
                  if (it[pdat::IntegerAttributeId(j)] != j) {
                     ++num_failures;
                  }
               }
            } // for (; it != ... (attribute iterator)
         } // for (; ic; ++ic) ... (cell iterator)
      } // for (hier::PatchLevel::iterator...
   } // hierarchy iteration

   bool it_passed = true;
   if (num_failures > 0) {
      it_passed = false;
      tbox::perr
      << "FAILED: - Test #1: Check iterator functionality."
      << std::endl;
   }

   return it_passed;
}

bool checkCopyOps(
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   const int data_id1, const int data_id2)
{

   bool copy_passed = true;

   typedef pdat::SparseData<pdat::CellGeometry> LSparseData;
   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& patch = *ip;
         std::shared_ptr<LSparseData> control(
            SAMRAI_SHARED_PTR_CAST<LSparseData, hier::PatchData>(
               patch->getPatchData(data_id1)));
         std::shared_ptr<LSparseData> copiedTo(
            SAMRAI_SHARED_PTR_CAST<LSparseData, hier::PatchData>(
               patch->getPatchData(data_id1)));
         std::shared_ptr<LSparseData> copiedFrom(
            SAMRAI_SHARED_PTR_CAST<LSparseData, hier::PatchData>(
               patch->getPatchData(data_id2)));
         TBOX_ASSERT(control);
         TBOX_ASSERT(copiedTo);
         TBOX_ASSERT(copiedFrom);

         int edit = copiedTo->size() / 2;
         LSparseData::iterator ct_it(copiedTo.get());
         for ( ; ct_it != copiedTo->end() && edit > 0; ++ct_it, --edit) {
         }

         while (ct_it != copiedTo->end()) {
            copiedTo->remove(ct_it);
         }
         copiedTo->copy(*copiedFrom);

         edit = copiedTo->size() / 2;
         ct_it = copiedTo->begin();

         LSparseData::iterator ctrl_it(control.get());
         bool first_passed = true;
         for (int i = 0; i < edit; ++i) {
            if (!ct_it.equals(ctrl_it)) {
               first_passed = false;
            } else {
               ++ct_it;
               ++ctrl_it;
            }
         }
         LSparseData::iterator cf_it(copiedFrom.get());
         edit = (copiedTo->size() / 2);
         for (int i = 0; i < edit; ++i) {
            ++cf_it;
         }

         bool second_passed = true;
         if (*copiedTo != *copiedFrom) {
            second_passed = false;
         }

         if (!first_passed) {
            tbox::perr
            << "FAILED: (first) copy for level " << ln << std::endl;
            copy_passed = false;
         }
         if (!second_passed) {
            tbox::perr
            << "FAILED: (second) copy for level " << ln << std::endl;
            copy_passed = false;
         }
      }
   }
   return copy_passed;
}

bool checkRemoveOps(
   const std::shared_ptr<hier::PatchHierarchy> hierarchy,
   const int data_id1)
{
   bool remove_passed = true;

   typedef pdat::SparseData<pdat::CellGeometry> LSparseData;
   int num_failures(0);
   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const std::shared_ptr<hier::Patch>& patch = *ip;

         std::shared_ptr<LSparseData> sample(
            SAMRAI_SHARED_PTR_CAST<LSparseData, hier::PatchData>(
               patch->getPatchData(data_id1)));
         TBOX_ASSERT(sample);

         LSparseData::iterator it;
         int stop = sample->size() / 2;
         for (it = sample->begin(); it != sample->end() && stop != 0;
              ++it, --stop) {
         }

         it = sample->end();
         if (it != sample->end()) {
            tbox::perr
            << "FAILED: - Test 3a: iterator fastforward"
            << std::endl;
         }

         it = sample->begin();
         if (it != sample->begin()) {
            tbox::perr
            << "FAILED: - Test 3b: iterator rewind."
            << std::endl;
         }

         // do a clear of the elements
         sample->clear();

         if (!sample->empty()) {
            ++num_failures;
            remove_passed = false;
            tbox::perr << "sample size is " << sample->size() << std::endl;
         }
      }
   }

   if (!remove_passed) {
      tbox::perr
      << "FAILED: the container is not empty and it should be."
      << std::endl;
   }
   return remove_passed;
}

std::shared_ptr<geom::CartesianGridGeometry>
getGeometry(
   hier::BoxContainer& coarse_domain,
   hier::BoxContainer& fine_domain,
   const tbox::Dimension& dim)
{
   std::vector<double> lo(dim.getValue());
   lo[0] = 0.0;
   lo[1] = 0.0;
   std::vector<double> hi(dim.getValue());
   hi[0] = 1.0;
   hi[1] = 0.5;

   // Sparse data sample 1 info
   hier::Box coarse0(hier::Index(0, 0), hier::Index(9, 2), hier::BlockId(0));
   hier::Box coarse1(hier::Index(0, 3), hier::Index(9, 4), hier::BlockId(0));
   hier::Box fine0(hier::Index(4, 4), hier::Index(7, 7), hier::BlockId(0));
   hier::Box fine1(hier::Index(8, 4), hier::Index(13, 7), hier::BlockId(0));

   coarse_domain.pushBack(coarse0);
   coarse_domain.pushBack(coarse1);
   fine_domain.pushBack(fine0);
   fine_domain.pushBack(fine1);

   std::shared_ptr<geom::CartesianGridGeometry> geometry(
      new geom::CartesianGridGeometry(
         "CartesianGeometry",
         &lo[0],
         &hi[0],
         coarse_domain));
   return geometry;
}

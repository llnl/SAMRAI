/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Main program for testing BoxContainer iterators
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxContainerSingleBlockIterator.h"
#include "SAMRAI/hier/BoxContainerSingleOwnerIterator.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"

using namespace SAMRAI;

int main(
   int argc,
   char* argv[])
{
   int fail_count = 0;

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      /*
       * Iterator tests.
       */

      const tbox::Dimension dim(2);

      hier::BoxContainer mboxes;

      const int num_boxes = 100;
      const int num_blocks = 5;
      const int num_owners = 10;

      // Build the BoxContainer.
      for (int i = 0; i < num_boxes; ++i) {

         int owner(i % num_owners);
         hier::BlockId bid(i / num_blocks);
         hier::LocalId lid(i);
         hier::BoxId mbid(lid, owner);

         hier::Box mb(dim, mbid);
         mb.setBlockId(bid);
         mboxes.insert(mb);

      }

      // Test 1: Block iterator.

      for (int b = 0; b < num_blocks; ++b) {

         const hier::BlockId bid(b);

         for (hier::BoxContainerSingleBlockIterator bi(mboxes.begin(hier::BlockId(b)));
              bi != mboxes.end(hier::BlockId(b)); ++bi) {

            if (bi->getBlockId() != bid) {
               tbox::perr << "FAILED: - Test #1: box id " << bi->getBlockId()
                          << " should have BlockId " << bid << std::endl;
               ++fail_count;
            }

         }
      }

      // Test 2: Owner iterator.

      for (int owner_rank = 0; owner_rank < num_owners; ++owner_rank) {

         for (hier::BoxContainerSingleOwnerIterator bi(mboxes.begin(owner_rank));
              bi != mboxes.end(owner_rank); ++bi) {

            if (bi->getOwnerRank() != owner_rank) {
               tbox::perr << "FAILED: - Test #2: box id " << bi->getBlockId()
                          << " should have rank " << owner_rank << std::endl;
               ++fail_count;
            }

         }
      }

      if (fail_count == 0) {
         tbox::pout << "\nPASSED:  testboxcontaineriterator" << std::endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return fail_count;
}

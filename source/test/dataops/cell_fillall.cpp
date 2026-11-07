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

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"

#include "SAMRAI/tbox/SAMRAIManager.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include <string>

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

using namespace SAMRAI;

int main(
   int argc,
   char* argv[]) {

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   const tbox::Dimension dim(2);

   int num_failures = 0;

   hier::Index lo(0,0);
   hier::Index hi(9,17);
   hier::IntVector ghosts(dim, 0);
   hier::Box box(lo, hi, hier::BlockId(0));
   pdat::CellData<int> cdata(box, 1, ghosts);

   cdata.fillAll(-1);

   pdat::CellIterator icend(pdat::CellGeometry::end(cdata.getGhostBox()));
   for (pdat::CellIterator c(pdat::CellGeometry::begin(cdata.getGhostBox()));
        c != icend; ++c) {
      if (cdata(*c) != -1) {
         ++num_failures;
      }
   }

   cdata.fillAll(5, cdata.getGhostBox());

   for (pdat::CellIterator c(pdat::CellGeometry::begin(cdata.getGhostBox()));
        c != icend; ++c) {
      if (cdata(*c) != 5) {
         ++num_failures;
      }
   }

   hier::Index lo_2(3,4);
   hier::Index hi_2(7,12);
   hier::Box box_2(lo_2, hi_2, hier::BlockId(0));

   cdata.fillAll(12, box_2);

   for (pdat::CellIterator c(pdat::CellGeometry::begin(cdata.getGhostBox()));
        c != icend; ++c) {
      if (box_2.contains(*c)) {
         if (cdata(*c) != 12) {
            ++num_failures;
         }
      } else {
         if (cdata(*c) != 5) {
            ++num_failures;
         }
      }
   }

   if (num_failures == 0) {
      tbox::pout << "\nPASSED:  cell fillall" << std::endl;
   } else {
      tbox::perr << "\nFAILED:  cell fillall" << std::endl;
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return num_failures;
}


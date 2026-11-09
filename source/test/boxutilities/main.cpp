/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Test program to demonstrate/test the Dimension class
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BlockId.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"

#include <list>
#include <vector>

using namespace SAMRAI;

namespace {

int
testCutPoints(
   const char* test_name,
   const int cell_count,
   std::list<int> cuts,
   const std::vector<int>& bad_cut_positions,
   const int min_size,
   const int cut_factor,
   const std::list<int>& expected_cuts)
{
   const tbox::Dimension dim(1);
   const hier::Box box(
      hier::Index(dim, 0),
      hier::Index(dim, cell_count - 1),
      hier::BlockId(0));
   std::vector<bool> bad_cuts(cell_count, false);
   for (const int cut_position : bad_cut_positions) {
      bad_cuts[cut_position] = true;
   }

   hier::BoxUtilities::fixBadCutPointsForDirection(
      0, cuts, bad_cuts, box, min_size, cut_factor);

   if (cuts != expected_cuts) {
      tbox::perr << "FAILED: " << test_name << "\nExpected cuts:";
      for (const int cut_position : expected_cuts) {
         tbox::perr << ' ' << cut_position;
      }
      tbox::perr << "\nActual cuts:";
      for (const int cut_position : cuts) {
         tbox::perr << ' ' << cut_position;
      }
      tbox::perr << std::endl;
      return 1;
   }

   return 0;
}

}

int main(
   int argc,
   char* argv[])
{
   int fail_count = 0;

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   {
      fail_count += testCutPoints(
         "empty cut list", 15, {}, {6}, 6, 1, {});

      fail_count += testCutPoints(
         "no bad cuts", 18, {6, 12}, {}, 6, 1, {6, 12});

      fail_count += testCutPoints(
         "single good cut with bad points elsewhere",
         15, {6}, {5}, 6, 1, {6});

      fail_count += testCutPoints(
         "single bad cut replaced below",
         15, {6}, {6}, 1, 1, {5});

      fail_count += testCutPoints(
         "single bad cut replaced above",
         15, {6}, {6}, 6, 1, {7});

      fail_count += testCutPoints(
         "single bad cut removed",
         15, {6}, {6, 7, 8, 9}, 6, 1, {});

      fail_count += testCutPoints(
         "first bad cut replaced below",
         18, {6, 12}, {6}, 1, 1, {5, 12});

      fail_count += testCutPoints(
         "first bad cut removed with upper cut retained",
         18, {6, 12}, {6}, 6, 1, {12});

      fail_count += testCutPoints(
         "upper bad cut removed after lower iterator advances onto it",
         18, {6, 12}, {12}, 6, 1, {6});

      fail_count += testCutPoints(
         "all cuts removed",
         18, {6, 12}, {6, 12}, 6, 1, {});

      fail_count += testCutPoints(
         "upper replacement becomes first remaining cut",
         19, {6, 12}, {6, 12}, 6, 1, {13});

      fail_count += testCutPoints(
         "upper successor becomes first remaining cut",
         30, {6, 12, 18, 24},
         {6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18},
         6, 1, {24});

      fail_count += testCutPoints(
         "shared bad cut replaced below with predecessor",
         12, {3, 6, 9}, {6}, 1, 1, {3, 5, 9});

      fail_count += testCutPoints(
         "shared bad cut replaced above with successor",
         12, {3, 6, 9}, {4, 5, 6}, 1, 1, {3, 7, 9});

      fail_count += testCutPoints(
         "replacement respects rounded minimum size and cut factor",
         16, {6}, {6}, 3, 2, {4});
   }

   if (fail_count == 0) {
      tbox::pout << "\nPASSED:  boxutilities" << std::endl;
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return fail_count;
}

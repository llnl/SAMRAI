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

// Headers for basic SAMRAI objects used in this code.
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>
#include <cassert>


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
      /* This test assumes that the maximum dim is at least 3 which is the default */
      assert(SAMRAI::MAX_DIM_VAL >= 3);

      tbox::Dimension dim1(1);
      tbox::Dimension dim2(2);
      tbox::Dimension dim3(3);

      if (dim1.getValue() != 1) {
         ++fail_count;
         TBOX_ERROR("Failed dim check; dim = 1");
      }

      if (dim2.getValue() != 2) {
         ++fail_count;
         TBOX_ERROR("Failed dim check; dim = 2");
      }

      if (dim3.getValue() != 3) {
         ++fail_count;
         TBOX_ERROR("Failed dim check; dim = 3");
      }

      tbox::Dimension a(2), b(2);
      if (!(a == b)) {
         ++fail_count;
         TBOX_ERROR("Failed dim comparison check: ==");
      }

      if (!(a != dim1)) {
         ++fail_count;
         TBOX_ERROR("Failed dim comparison check: !=");
      }

      if (dim1 > dim2) {
         ++fail_count;
         TBOX_ERROR("Failed dim comparison check: >");
      }

      if (!(dim3 > dim2)) {
         ++fail_count;
         TBOX_ERROR("Failed dim comparison check: >");
      }

      if (dim1 >= dim2) {
         ++fail_count;
         TBOX_ERROR("Failed dim comparison check: >=");
      }

      if (!(dim3 >= dim2)) {
         ++fail_count;
         TBOX_ERROR("Failed dim comparison check: >=");
      }

      if (!(dim2 >= dim2)) {
         ++fail_count;
         TBOX_ERROR("Failed dim comparison check: >=");
      }

      if (dim2 < dim1) {
         ++fail_count;
         TBOX_ERROR("Failed dim comparison check: <");
      }

      if (!(dim2 < dim3)) {
         ++fail_count;
         TBOX_ERROR("Failed dim comparison check: <");
      }

      if (dim2 <= dim1) {
         ++fail_count;
         TBOX_ERROR("Failed dim comparison check: <=");
      }

      if (!(dim2 <= dim3)) {
         ++fail_count;
         TBOX_ERROR("Failed dim comparison check: <=");
      }

      if (!(dim2 <= dim2)) {
         ++fail_count;
         TBOX_ERROR("Failed dim comparison check: <=");
      }

#if 0
      // Currently not allowed.
      a = dim3;
      if (a != dim3) {
         ++fail_count;
         TBOX_ERROR("Failed dim assignment check");
      }
#endif

   }

   if (fail_count == 0) {
      tbox::pout << "\nPASSED:  dimension" << std::endl;
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();
   return fail_count;
}

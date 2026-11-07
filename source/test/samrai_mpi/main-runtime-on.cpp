/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Test program for SAMRAI_MPI with run-time MPI on.
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "mpi-interface-tests.h"

using namespace SAMRAI;
using namespace tbox;

/*
 ************************************************************************
 * Test run-time MPI mode on.
 ************************************************************************
 */

int main(
   int argc,
   char* argv[])
{
   int fail_count = 0;

   bool runtime_mpi = true;
   bool disable_mpi = false;
   mpiInterfaceTests(argc, argv, fail_count, runtime_mpi, disable_mpi);

   if (fail_count == 0) {
      if (tbox::SAMRAI_MPI::usingMPI()) {
         tbox::pout << "\nPASSED:  " << argv[0] << std::endl;
      } else {
         std::cout << "\nPASSED:  " << argv[0] << std::endl;
      }
   }

   return fail_count;
}

/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Routines for tracking memory use in SAMRAI.
 *
 ************************************************************************/

#include "SAMRAI/tbox/MemoryUtilities.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/IOStream.h"

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace tbox {

double MemoryUtilities::s_max_memory = 0.;

/*
 *************************************************************************
 *
 * Prints memory usage to specified output stream.  Each time this
 * method is called, it prints in the format:
 *
 *    253.0MB (265334688) in 615 allocs, 253.9MB reserved (871952 unused)
 *
 * where
 *
 *    253.0MB is how many megabytes your current allocation has malloced.
 *    2653346688 is the precise size (in bytes) of your current alloc.
 *    615 is the number of items allocated with malloc.
 *    253.9MB is the current memory reserved by the system for mallocs.
 *    871952 is the bytes currently not used in this reserved memory.
 *
 *************************************************************************
 */
void
MemoryUtilities::printMemoryInfo(
   std::ostream& os)
{
   NULL_USE(os);

#ifdef HAVE_MALLINFO
   /*
    * NOTE: This was taken directly from John Gyllenhal...
    */

   /* Get malloc info structure */
   struct mallinfo my_mallinfo = mallinfo();

   /* Get total memory reserved by the system for malloc currently*/
   double reserved_mem = my_mallinfo.arena;

   /* Get all the memory currently allocated to user by malloc, etc. */
   double used_mem = my_mallinfo.hblkhd + my_mallinfo.usmblks
      + my_mallinfo.uordblks;

   /* Get memory not currently allocated to user but malloc controls */
   double free_mem = my_mallinfo.fsmblks + my_mallinfo.fordblks;

   /* Get number of items currently allocated */
   double number_allocated = my_mallinfo.ordblks + my_mallinfo.smblks;

   /* Record high-water mark for memory used. */
   s_max_memory = MathUtilities<double>::Max(s_max_memory, used_mem);

   /* Print out concise malloc info line */
   os << used_mem / (1024.0 * 1024.0) << "MB ("
      << used_mem << ") in "
      << number_allocated << " allocs, "
      << reserved_mem / (1024.0 * 1024.0) << "MB reserved ("
      << free_mem << " unused)" << std::endl;

#endif
}

/*
 *************************************************************************
 *
 * Prints maximum memory used (i.e. high-water mark).  The max is
 * determined each time the "printMemoryInfo" or "recordMemoryInfo"
 * functions are called.
 *
 *************************************************************************
 */
void
MemoryUtilities::printMaxMemory(
   std::ostream& os)
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   /*
    * Step through all nodes (>0) and send max memory to processor 0,
    * which subsequently writes it out.
    */
   int maxmem = 0;
   int len = 1;
   SAMRAI_MPI::Status status;
   for (int p = 0; p < mpi.getSize(); ++p) {
      if (mpi.getSize() > 1) {
         if (mpi.getRank() == p) {
            maxmem = static_cast<int>(s_max_memory);
            mpi.Send(&maxmem, len, MPI_INT, 0, 0);
         }
         if (mpi.getRank() == 0) {
            mpi.Recv(&maxmem, len, MPI_INT, p, 0, &status);
         }
      }
      os << "Maximum memory used on processor " << p
         << ": " << maxmem / (1024. * 1024.) << " MB" << std::endl;
   }

}

size_t
MemoryUtilities::align(
   const size_t bytes)
{
   size_t aligned = bytes + ArenaAllocationAlignment - 1;
   aligned -= aligned % ArenaAllocationAlignment;
   return aligned;
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Unsuppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

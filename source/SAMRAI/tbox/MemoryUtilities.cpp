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
#include <vector>
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

#if defined(HAVE_MALLINFO2) || defined(HAVE_MALLINFO)
   /*
    * NOTE: This was taken directly from John Gyllenhal...
    */

   /*
    * Get malloc info structure.  mallinfo2() has size_t fields; the
    * older mallinfo() has int fields, which wrap around above 2 GiB.
    */
#ifdef HAVE_MALLINFO2
   struct mallinfo2 my_mallinfo = mallinfo2();
#else
   struct mallinfo my_mallinfo = mallinfo();
#endif

   /*
    * Each field is converted to double before summing so that the sums
    * cannot overflow the integer type of the fields.
    */

   /* Get total memory reserved by the system for malloc currently*/
   double reserved_mem = static_cast<double>(my_mallinfo.arena);

   /* Get all the memory currently allocated to user by malloc, etc. */
   double used_mem = static_cast<double>(my_mallinfo.hblkhd)
      + static_cast<double>(my_mallinfo.usmblks)
      + static_cast<double>(my_mallinfo.uordblks);

   /* Get memory not currently allocated to user but malloc controls */
   double free_mem = static_cast<double>(my_mallinfo.fsmblks)
      + static_cast<double>(my_mallinfo.fordblks);

   /* Get number of items currently allocated */
   double number_allocated = static_cast<double>(my_mallinfo.ordblks)
      + static_cast<double>(my_mallinfo.smblks);

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
    * Gather the high-water mark of every processor on processor 0,
    * which subsequently writes them out.  The values are in bytes and
    * are kept as double since they can exceed the range of int.
    */
   const int nprocs = mpi.getSize();
   std::vector<double> maxmem(nprocs, 0.);
   double local_maxmem = s_max_memory;
   if (SAMRAI_MPI::usingMPI() && nprocs > 1) {
      mpi.Gather(&local_maxmem, 1, MPI_DOUBLE,
         maxmem.data(), 1, MPI_DOUBLE, 0);
   } else {
      maxmem[0] = local_maxmem;
   }

   if (mpi.getRank() == 0) {
      for (int p = 0; p < nprocs; ++p) {
         os << "Maximum memory used on processor " << p
            << ": " << maxmem[p] / (1024. * 1024.) << " MB" << std::endl;
      }
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

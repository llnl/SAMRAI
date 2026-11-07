/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Routine for tracking memory use in applications.
 *
 ************************************************************************/

#ifndef included_tbox_MemoryUtilities
#define included_tbox_MemoryUtilities

#include "SAMRAI/SAMRAI_config.h"

#ifndef included_iostream
#include <iostream>
#endif

namespace SAMRAI {
namespace tbox {

/*!
 * @brief Class MemoryUtilities provides utility methods for memory usage and
 * misc memory utilities.
 *
 * "printMemoryInfo()" which does a simple dump of the current memory
 * usage on a processor, and "recordMemoryInfo()" which records the
 * memory for post-process analysis.
 *
 * Calls to these methods may be placed at various points in an application
 * to track memory usage characteristics.  For applications running on a
 * single processor, the  call the print method is likely sufficient. The
 * information can simply be printed to a log file or output stream. 
 *
 * Note that all member functions of this class are static so it is not
 * necessary to instantiate the class.  Simply call the functions as
 * static functions; e.g.,MemoryUtilities::function(...).
 */
struct MemoryUtilities {
   /*!
    * Print memory information to the supplied output stream.
    */
   static void
   printMemoryInfo(
      std::ostream& os);

   /*!
    * Print maximum memory used (i.e. high-water mark) to the
    * supplied output stream.
    */
   static void
   printMaxMemory(
      std::ostream& os);

   /**
    * Static function to compute alignment for memory allocation.
    * Data allocations less than the alignment size are rounded up to
    * the next multiple of the allocation size.  All data allocations
    * are aligned on 16 byte boundaries.  Thus, a memory allocation of
    * only 9 bytes will actually return a 16 byte chunk of memory.
    */
   static size_t
   align(
      const size_t bytes);

private:
   /*
    * Keep track of maximum memory used (updated each time print function
    * called).
    */
   static double s_max_memory;

   enum { ArenaAllocationAlignment = 16 };
};

}
}

#endif

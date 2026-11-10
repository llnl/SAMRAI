/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Test of MemoryUtilities memory reporting.
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/MemoryUtilities.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"

#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

using namespace SAMRAI;

/*
 * Each processor allocates a buffer (64 MiB by default, or the number of
 * MiB given as first argument), records its memory with printMemoryInfo
 * and calls printMaxMemory.  Processor 0 must then write one line per
 * processor, each reporting at least the size of the buffer, and the
 * other processors must write nothing.
 *
 * Passing a size above 2048 MiB checks that values beyond the range of
 * int are reported correctly.
 */
int main(
   int argc,
   char* argv[])
{
   int fail_count = 0;

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   {
      const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

      const double buffer_mib = argc > 1 ? std::atof(argv[1]) : 64.;
      const size_t buffer_bytes =
         static_cast<size_t>(buffer_mib * 1024. * 1024.);
      std::vector<char> buffer(buffer_bytes, 1);
      // Let the buffer escape so that its allocation cannot be elided.
      tbox::plog << "buffer at " << static_cast<const void *>(buffer.data())
                 << ", " << buffer_bytes << " bytes" << std::endl;

      std::ostringstream info;
      tbox::MemoryUtilities::printMemoryInfo(info);
      tbox::plog << info.str();

      std::ostringstream report;
      tbox::MemoryUtilities::printMaxMemory(report);
      tbox::pout << report.str();

#if defined(HAVE_MALLINFO2) || defined(HAVE_MALLINFO)
      if (mpi.getRank() == 0) {
         std::istringstream lines(report.str());
         std::string line;
         int nlines = 0;
         while (std::getline(lines, line)) {
            const std::string prefix = "Maximum memory used on processor "
               + std::to_string(nlines) + ": ";
            if (line.compare(0, prefix.size(), prefix) != 0) {
               tbox::perr << "FAILED: unexpected line: " << line << std::endl;
               ++fail_count;
            } else {
               const double mb = std::atof(line.c_str() + prefix.size());
               if (mb < buffer_mib) {
                  tbox::perr << "FAILED: processor " << nlines
                             << " reports " << mb << " MB, expected at least "
                             << buffer_mib << " MB" << std::endl;
                  ++fail_count;
               }
            }
            ++nlines;
         }
         if (nlines != mpi.getSize()) {
            tbox::perr << "FAILED: " << nlines << " lines for "
                       << mpi.getSize() << " processors" << std::endl;
            ++fail_count;
         }
      } else if (!report.str().empty()) {
         tbox::perr << "FAILED: processor " << mpi.getRank()
                    << " wrote the report" << std::endl;
         ++fail_count;
      }
#else
      tbox::pout << "mallinfo not available, memory values not checked"
                 << std::endl;
#endif

      if (fail_count == 0) {
         tbox::pout << "\nPASSED:  memory_utilities" << std::endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();
   return fail_count;
}

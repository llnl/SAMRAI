/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Utility functions for error reporting, file manipulation, etc.
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/Logger.h"
#include "SAMRAI/tbox/PIO.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace tbox {

/*
 * Routine to recursively construct directories based on a relative path name.
 */
void
Utilities::recursiveMkdir(
   const std::string& path,
   mode_t mode,
   bool only_node_zero_creates)
{

#ifdef _MSC_VER
   const char seperator = '/';
#define mkdir(path, mode) mkdir(path)
#else
   const char seperator = '/';
#endif

   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   if ((!only_node_zero_creates) || (mpi.getRank() == 0)) {
      int length = static_cast<int>(path.length());
      char* path_buf = new char[length + 1];
      snprintf(path_buf, length + 1, "%s", path.c_str());
      struct stat status;
      int pos = length - 1;

      /* find part of path that has not yet been created */
      while ((stat(path_buf, &status) != 0) && (pos >= 0)) {

         /* slide backwards in string until next slash found */
         bool slash_found = false;
         while ((!slash_found) && (pos >= 0)) {
            if (path_buf[pos] == seperator) {
               slash_found = true;
               if (pos >= 0) path_buf[pos] = '\0';
            } else --pos;
         }
      }

      /*
       * if there is a part of the path that already exists make sure
       * it is really a directory
       */
      if (pos >= 0) {
         if (!S_ISDIR(status.st_mode)) {
            TBOX_ERROR("Error in Utilities::recursiveMkdir...\n"
               << "    Cannot create directories in path = " << path
               << "\n    because some intermediate item in path exists and"
               << "is NOT a directory" << std::endl);
         }
      }

      /* make all directories that do not already exist */

      /*
       * if (pos < 0), then there is no part of the path that
       * already exists.  Need to make the first part of the
       * path before sliding along path_buf.
       */
      if (pos < 0) {
         if (mkdir(path_buf, mode) != 0) {
            TBOX_ERROR("Error in Utilities::recursiveMkdir...\n"
               << "    Cannot create directory  = "
               << path_buf << std::endl);
         }
         pos = 0;
      }

      /* make rest of directories */
      do {

         /* slide forward in string until next '\0' found */
         bool null_found = false;
         while ((!null_found) && (pos < length)) {
            if (path_buf[pos] == '\0') {
               null_found = true;
               path_buf[pos] = seperator;
            }
            ++pos;
         }

         /* make directory if not at end of path */
         if (pos < length) {
            if (mkdir(path_buf, mode) != 0) {
               TBOX_ERROR("Error in Utilities::recursiveMkdir...\n"
                  << "    Cannot create directory  = "
                  << path_buf << std::endl);
            }
         }
      } while (pos < length);

      delete[] path_buf;
   }

   /*
    * Make sure all processors wait until node zero creates
    * the directory structure.
    */
   if (only_node_zero_creates) {
      SAMRAI_MPI::getSAMRAIWorld().Barrier();
   }
}

/*
 * Routine to convert an integer to a string.
 */
std::string
Utilities::intToString(
   int num,
   int min_width)
{
   int tmp_width = (min_width > 0 ? min_width : 1);
   std::ostringstream os;
   if (num < 0) {
      os << '-' << std::setw(tmp_width - 1) << std::setfill('0') << -num;
   } else {
      os << std::setw(tmp_width) << std::setfill('0') << num;
   }
   os << std::flush;

   return os.str();  //returns the string form of the stringstream object
}

/*
 * Routine to convert a size_t to a string.
 */
std::string
Utilities::sizetToString(
   size_t num,
   int min_width)
{
   int tmp_width = (min_width > 0 ? min_width : 1);
   std::ostringstream os;
   os << std::setw(tmp_width) << std::setfill('0') << num;
   os << std::flush;

   return os.str();  //returns the string form of the stringstream object
}

/*
 * Routine that calls abort and prints calling location to error stream.
 */
void
Utilities::abort(
   const std::string& message,
   const std::string& filename,
   const int line)
{
   Logger::getInstance()->logAbort(message, filename, line);

   SAMRAI_MPI::abort();
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Parallel I/O classes pout, perr, and plog and control class
 *
 ************************************************************************/

#include <string>

#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/ParallelBuffer.h"

namespace SAMRAI {
namespace tbox {

int PIO::s_rank = -1;
std::ofstream * PIO::s_filestream = 0;

/*
 *************************************************************************
 *
 * Define the parallel buffers and the associated ostream objects.
 *
 *************************************************************************
 */

static ParallelBuffer pout_buffer;
static ParallelBuffer perr_buffer;
static ParallelBuffer plog_buffer;

std::ostream pout(&pout_buffer);
std::ostream perr(&perr_buffer);
std::ostream plog(&plog_buffer);

/*
 *************************************************************************
 *
 * Initialize the parallel I/O streams.  This routine must be called
 * before pout, perr, and plog are used for output but after SAMRAI_MPI
 * has been initialized.  By default, logging is disabled.
 *
 *************************************************************************
 */

void
PIO::initialize()
{
   const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
   mpi.Comm_rank(&s_rank);
   s_filestream = 0;

   /*
    * Initialize the standard parallel output stream
    */

   pout_buffer.setActive(s_rank == 0);
   pout_buffer.setPrefixString(std::string());
   pout_buffer.setOutputStream1(&std::cout);
   pout_buffer.setOutputStream2(0);

   /*
    * Initialize the error parallel output stream
    */

   std::string buffer = "P=" + Utilities::processorToString(s_rank) + ":";

   perr_buffer.setActive(true);
   perr_buffer.setPrefixString(buffer);
   perr_buffer.setOutputStream1(&std::cerr);
   perr_buffer.setOutputStream2(0);

   /*
    * Initialize the parallel log file (disabled by default)
    */

   plog_buffer.setActive(false);
   plog_buffer.setPrefixString(std::string());
   plog_buffer.setOutputStream1(0);
   plog_buffer.setOutputStream2(0);
}

/*
 *************************************************************************
 *
 * Close the output streams.  Flush both cout and cerr.  If logging,
 * then flush and close the log stream.
 *
 *************************************************************************
 */

void
PIO::finalize()
{
   std::cout.flush();
   std::cerr.flush();
   shutdownFilestream();
}

/*
 *************************************************************************
 *
 * If the log file stream is open, then shut down the filestream.  Close
 * and flush the channel and disconnect the output stream buffers.
 *
 *************************************************************************
 */

void
PIO::shutdownFilestream()
{
   if (s_filestream) {
      s_filestream->flush();
      s_filestream->close();

      delete s_filestream;
      s_filestream = 0;

      pout_buffer.setOutputStream2(0);
      perr_buffer.setOutputStream2(0);
      plog_buffer.setOutputStream1(0);
      plog_buffer.setActive(false);
   }
}

/*
 *************************************************************************
 *
 * Log messages for node zero only.  If a log stream was open, close
 * it.  If this is node zero, then open a new log stream and set the
 * appropriate buffer streams to point to the log file.
 *
 *************************************************************************
 */

void
PIO::logOnlyNodeZero(
   const std::string& filename)
{
   /*
    * If the filestream was open, then close it and reset streams
    */

   shutdownFilestream();

   /*
    * If this is node zero, then open the log stream and redirect output
    */

   if (s_rank == 0) {
      s_filestream = new std::ofstream(filename.c_str());
      if (!(*s_filestream)) {
         delete s_filestream;
         s_filestream = 0;
         perr << "PIO: Could not open log file ``" << filename.c_str()
              << "''\n";
      } else {
         pout_buffer.setOutputStream2(s_filestream);
         perr_buffer.setOutputStream2(s_filestream);
         plog_buffer.setOutputStream1(s_filestream);
         plog_buffer.setActive(true);
      }
   }
}

/*
 *************************************************************************
 *
 * Log messages for all nodes.  If a log stream was open, the close it.
 * Open a log stream on every processor.  The filename for the log file
 * will be appended with the processor number.
 *
 *************************************************************************
 */

void
PIO::logAllNodes(
   const std::string& filename)
{
   /*
    * If the filestream was open, then close it and reset streams
    */

   shutdownFilestream();

   /*
    * Open the log stream and redirect output
    */

   std::string full_filename = filename + "."
      + Utilities::processorToString(s_rank);
   s_filestream = new std::ofstream(full_filename.c_str());

   if (!(*s_filestream)) {
      delete s_filestream;
      s_filestream = 0;
      perr << "PIO: Could not open log file ``" << full_filename << "''\n";
   } else {
      pout_buffer.setOutputStream2(s_filestream);
      perr_buffer.setOutputStream2(s_filestream);
      plog_buffer.setOutputStream1(s_filestream);
      plog_buffer.setActive(true);
   }
}

/*
 *************************************************************************
 *
 * Suspend logging of data to the file stream.  This does not close the
 * filestream (assuming it is open) but just disables logging.
 *
 *************************************************************************
 */
void
PIO::suspendLogging()
{
   pout_buffer.setOutputStream2(0);
   perr_buffer.setOutputStream2(0);
   plog_buffer.setOutputStream1(0);
   plog_buffer.setActive(false);
}

/*
 *************************************************************************
 *
 * Resume logging of the file stream (assuming it was open).  If the
 * file stream is 0, then do nothing.
 *
 *************************************************************************
 */
void
PIO::resumeLogging()
{
   if (s_filestream) {
      pout_buffer.setOutputStream2(s_filestream);
      perr_buffer.setOutputStream2(s_filestream);
      plog_buffer.setOutputStream1(s_filestream);
      plog_buffer.setActive(true);
   }
}

}
}

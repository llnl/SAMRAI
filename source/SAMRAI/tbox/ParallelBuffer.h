/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Parallel I/O class buffer to manage parallel ostreams output
 *
 ************************************************************************/

#ifndef included_tbox_ParallelBuffer
#define included_tbox_ParallelBuffer

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/OpenMPUtilities.h"

#include <iostream>
#include <string>

namespace SAMRAI {
namespace tbox {

/**
 * Class ParallelBuffer is a simple I/O stream utility that
 * intercepts output from an ostream and redirects the output as necessary
 * for parallel I/O.  This class defines a stream buffer class for an
 * ostream class.
 */

class ParallelBuffer:public std::streambuf
{
public:
   /**
    * Create a parallel buffer class.  The object will require further
    * initialization to set up the I/O streams and prefix string.
    */
   ParallelBuffer();

   /**
    * The destructor simply deallocates any internal data
    * buffers.  It does not modify the output streams.
    */
   ~ParallelBuffer();

   /**
    * Set whether the output stream will be active.  If the parallel buffer
    * stream is disabled, then no data is forwarded to the output streams.
    * The internal data buffer is deallocated and pointers are reset
    * whenever the parallel buffer is deactivated.
    */
   void
   setActive(
      bool active);

   /**
    * Set the prefix that begins every new line to the output stream.
    * A sample prefix is "P=XXXXX: ", where XXXXX represents the node
    * number.
    */
   void
   setPrefixString(
      const std::string& text)
   {
      d_prefix = text;
   }

   /**
    * Set the primary output stream.  If not NULL, then output data is
    * sent to this stream.  The primary output stream is typically stderr
    * or stdout or perhaps a log file.
    */
   void
   setOutputStream1(
      std::ostream* stream)
   {
      d_ostream1 = stream;
   }

   /**
    * Set the secondary output stream.  If not NULL, then output data is sent
    * to this stream.  The secondary output stream is typically NULL or a log
    * file that mirrors the primary output stream.
    */
   void
   setOutputStream2(
      std::ostream* stream)
   {
      d_ostream2 = stream;
   }

   /**
    * Write a text string to the output stream.  Note that the string is
    * not actually written until an end-of-line is detected.
    */
   void
   outputString(
      const std::string& text)
   {
      outputString(text, static_cast<int>(text.length()));
   }

   /**
    * Write a text string of the specified length to the output file.  Note
    * that the string is not actually written until an end-of-line is detected.
    */
   void
   outputString(
      const std::string& text,
      const int length,
      const bool recursive = false);

   /**
    * Synchronize the parallel buffer (called from streambuf).
    */
   int
   sync() override;

#if !defined(__INTEL_COMPILER) && (defined(__GNUG__))
   /**
    * Write the specified number of characters into the output stream (called
    * from streambuf).
    */
   std::streamsize
   xsputn(
      const char* text,
      std::streamsize n) override;
#endif

   /**
    * Write an overflow character into the parallel buffer (called from
    * streambuf).
    */
   int
   overflow(
      int ch) override;

#ifdef _MSC_VER

   /**
    * Read an overflow character from the parallel buffer (called from
    * streambuf).  This is not implemented.  It is needed by the
    * MSVC++ stream implementation.
    */
   int
   underflow() override
   {
      return EOF;
   }
#endif

private:
   // Unimplemented copy constructor.
   ParallelBuffer(
      const ParallelBuffer& other);

   // Unimplemented assignment operator.
   ParallelBuffer&
   operator = (
      const ParallelBuffer& rhs);

   void
   copyToBuffer(
      const std::string& text,
      const int length);
   void
   outputBuffer();              // output internal buffer data to streams

   bool d_active;               // whether this output stream is active
   std::string d_prefix;        // string prefix to prepend output strings
   std::ostream* d_ostream1;    // primary output stream for buffer
   std::ostream* d_ostream2;    // secondary output stream (e.g., for log file)
   char* d_buffer;              // internal buffer to store accumulated string
   int d_buffer_size;           // size of the internal output buffer
   int d_buffer_ptr;            // number of charcters in the output buffer
   TBOX_omp_lock_t l_buffer;    // OpenMP lock for buffer operations.

   static const int DEFAULT_BUFFER_SIZE;
};

}
}

#endif

/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Parallel I/O class buffer to manage parallel ostreams output
 *
 ************************************************************************/

#include "SAMRAI/tbox/ParallelBuffer.h"

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <string>
#include <cstring>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace tbox {

const int ParallelBuffer::DEFAULT_BUFFER_SIZE = 128;

/*
 *************************************************************************
 *
 * Construct a parallel buffer object.  The object will require further
 * initialization to set up I/O streams and the prefix string.
 *
 *************************************************************************
 */

ParallelBuffer::ParallelBuffer()
{
   d_active = true;
   d_prefix = std::string();
   d_ostream1 = 0;
   d_ostream2 = 0;
   d_buffer = 0;
   d_buffer_size = 0;
   d_buffer_ptr = 0;
   TBOX_omp_init_lock(&l_buffer);
}

/*
 *************************************************************************
 *
 * The destructor deallocates internal data buffer.  It does not modify
 * the output streams.
 *
 *************************************************************************
 */

ParallelBuffer::~ParallelBuffer()
{
   if (d_buffer) {
      delete[] d_buffer;
   }
   TBOX_omp_destroy_lock(&l_buffer);
}

/*
 *************************************************************************
 *
 * Activate or deactivate the output stream.  If the stream has been
 * deactivated, then deallocate the internal data buffer.
 *
 *************************************************************************
 */

void
ParallelBuffer::setActive(
   bool active)
{
   TBOX_omp_set_lock(&l_buffer);
   if (!active && d_buffer) {
      delete[] d_buffer;
      d_buffer = 0;
      d_buffer_size = 0;
      d_buffer_ptr = 0;
   }
   d_active = active;
   TBOX_omp_unset_lock(&l_buffer);
}

/*
 *************************************************************************
 *
 * Write a text string of the specified length to the output stream.
 * Note that the string data is accumulated into the internal output
 * buffer until an end-of-line is detected.
 *
 *************************************************************************
 */

void
ParallelBuffer::outputString(
   const std::string& text,
   const int length,
   const bool recursive)
{
#ifndef HAVE_OPENMP
   NULL_USE(recursive);
#endif

   if ((length > 0) && d_active) {

#ifdef HAVE_OPENMP 
      if (!recursive) {
         TBOX_omp_set_lock(&l_buffer);
      }
#endif

      /*
       * If we need to allocate the internal buffer, then do so
       */

      if (!d_buffer) {
         d_buffer = new char[DEFAULT_BUFFER_SIZE];
         d_buffer_size = DEFAULT_BUFFER_SIZE;
         d_buffer_ptr = 0;
      }

      /*
       * If the buffer pointer is zero, then prepend the prefix if not empty
       */

      if ((d_buffer_ptr == 0) && !d_prefix.empty()) {
         copyToBuffer(d_prefix, static_cast<int>(d_prefix.length()));
      }

      /*
       * Search for an end-of-line in the string
       */

      int eol_ptr = 0;
      for ( ; (eol_ptr < length) && (text[eol_ptr] != '\n'); ++eol_ptr)
         NULL_STATEMENT;

      /*
       * If no end-of-line found, copy the entire text string but no output
       */

      if (eol_ptr == length) {
         copyToBuffer(text, length);

         /*
          * If we found end-of-line, copy and output; recurse if more chars
          */

      } else {
         const int ncopy = eol_ptr + 1;
         copyToBuffer(text, ncopy);
         outputBuffer();
         if (ncopy < length) {
            outputString(text.substr(ncopy), length - ncopy, true);
         }
      }

#ifdef HAVE_OPENMP
      if (!recursive) {
         TBOX_omp_unset_lock(&l_buffer);
      }
#endif

   }
}

/*
 *************************************************************************
 *
 * Copy data from the text string into the internal output buffer.
 * If the internal buffer is not large enough to hold all of the string
 * data, then allocate a new internal buffer.
 *
 * This method is not thread-safe, but it is only called from
 * outputString(), which prevents multiple thread access to the
 * buffer.
 *
 *************************************************************************
 */

void
ParallelBuffer::copyToBuffer(
   const std::string& text,
   const int length)
{
   /*
    * First check whether we need to increase the size of the buffer
    */

   if (d_buffer_ptr + length > d_buffer_size) {
      const int new_size =
         MathUtilities<int>::Max(d_buffer_ptr + length, 2 * d_buffer_size);
      char* new_buffer = new char[new_size];

      if (d_buffer_ptr > 0) {
         (void)strncpy(new_buffer, d_buffer, d_buffer_ptr);
      }
      delete[] d_buffer;

      d_buffer = new_buffer;
      d_buffer_size = new_size;
   }

   /*
    * Copy data from the input into the internal buffer and increment pointer
    */

   TBOX_ASSERT(d_buffer_ptr + length <= d_buffer_size);

   strncpy(d_buffer + d_buffer_ptr, text.c_str(), length);
   d_buffer_ptr += length;
}

/*
 *************************************************************************
 *
 * Output buffered stream data to the active output streams and reset
 * the buffer pointer to its empty state.
 *
 *************************************************************************
 */

void
ParallelBuffer::outputBuffer()
{
   if (d_buffer_ptr > 0) {
      if (d_ostream1) {
         d_ostream1->write(d_buffer, d_buffer_ptr);
         d_ostream1->flush();
      }
      if (d_ostream2) {
         d_ostream2->write(d_buffer, d_buffer_ptr);
         d_ostream2->flush();
      }
      d_buffer_ptr = 0;
   }
}

/*
 *************************************************************************
 *
 * Synchronize the parallel buffer and write string data.  This routine
 * is called from streambuf.
 *
 *************************************************************************
 */

int
ParallelBuffer::sync()
{
   const int n = static_cast<int>(pptr() - pbase());
   if (n > 0) {
      outputString(pbase(), n);
   }
   return 0;
}

/*
 *************************************************************************
 *
 * Write the specified number of characters into the output stream.
 * This routine is called from streambuf.  If this routine is not
 * provided, then overflow() is called instead for each character.
 *
 * Note that this routine is not required; it only
 * offers some efficiency over overflow().
 *
 *************************************************************************
 */

#if !defined(__INTEL_COMPILER) && (defined(__GNUG__))
std::streamsize
ParallelBuffer::xsputn(
   const char* text,
   std::streamsize n)
{
   sync();
   if (n > 0) outputString(std::string(text, n), static_cast<int>(n));
   return n;
}
#endif

/*
 *************************************************************************
 *
 * Write a single character into the parallel buffer.  This routine is
 * called from streambuf.
 *
 *************************************************************************
 */

int
ParallelBuffer::overflow(
   int ch)
{
   const int n = static_cast<int>(pptr() - pbase());
   if (n && sync()) {
      return EOF;
   }
   if (ch != EOF) {
      char character[2];
      character[0] = (char)ch;
      character[1] = 0;
      outputString(character, 1);
   }
   pbump(-n);
   return 0;
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

/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Fixed-size message buffer used in interprocessor communication
 *
 ************************************************************************/
#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/tbox/Utilities.h"

#include "SAMRAI/tbox/AllocatorDatabase.h"

namespace SAMRAI {
namespace tbox {

/*
 *************************************************************************
 *
 * The constructor and destructor for MessageStream.
 *
 *************************************************************************
 */

MessageStream::MessageStream(
   const size_t num_bytes,
   const StreamMode mode,
   const void* data_to_read,
   bool deep_copy
#ifdef HAVE_UMPIRE
   , umpire::TypedAllocator<char> allocator
#endif
   ):
   d_mode(mode),
#ifdef HAVE_UMPIRE
   d_write_buffer(allocator),
#else
   d_write_buffer(),
#endif
   d_read_buffer(0),
#ifdef HAVE_UMPIRE
   d_allocator(allocator),
#endif
   d_buffer_size(0),
   d_buffer_index(0),
   d_grow_as_needed(false),
   d_deep_copy_read(deep_copy)
{
   TBOX_ASSERT(num_bytes >= 1);

   if (mode == Read) {
      if (data_to_read == 0) {
         TBOX_ERROR("MessageStream::MessageStream: error:\n"
            << "No data_to_read was given to a Read-mode MessageStream.\n");
      }
      if (deep_copy) {
#ifdef HAVE_UMPIRE
         d_read_buffer = allocator.allocate(num_bytes);
#else
         d_read_buffer = new char[num_bytes];
#endif
         memcpy(const_cast<char *>(d_read_buffer), data_to_read, num_bytes);
      } else {
         d_read_buffer = static_cast<const char *>(data_to_read);
      }
      d_buffer_size = num_bytes;
   } else {
      d_write_buffer.reserve(num_bytes);
   }
}

#ifdef HAVE_UMPIRE
MessageStream::MessageStream(umpire::TypedAllocator<char> allocator):
   d_mode(Write),
   d_write_buffer(allocator),
   d_read_buffer(0),
   d_allocator(allocator),
   d_buffer_size(0),
   d_buffer_index(0),
   d_grow_as_needed(true),
   d_deep_copy_read(false)
{
   d_write_buffer.reserve(10);
}
#endif

MessageStream::MessageStream():
   d_mode(Write),
#ifdef HAVE_UMPIRE
   d_write_buffer(AllocatorDatabase::getDatabase()->getInternalHostAllocator()),
#else
   d_write_buffer(),
#endif
   d_read_buffer(0),
#ifdef HAVE_UMPIRE
   d_allocator(AllocatorDatabase::getDatabase()->getInternalHostAllocator()),
#endif
   d_buffer_size(0),
   d_buffer_index(0),
   d_grow_as_needed(true),
   d_deep_copy_read(false)
{
   d_write_buffer.reserve(10);
}

MessageStream::~MessageStream()
{
   if (d_mode == Read && d_deep_copy_read) {
#ifdef HAVE_UMPIRE
      d_allocator.deallocate((char*)d_read_buffer, d_buffer_size);   
#else
      delete[] d_read_buffer;
#endif
   }
   d_read_buffer = 0;
}

/*
 *************************************************************************
 *
 * Print out class data if an assertion is thrown.
 *
 *************************************************************************
 */

void
MessageStream::printClassData(
   std::ostream& os) const
{
   os << "Maximum buffer size = " << d_buffer_size << std::endl;
   os << "Current buffer index = " << d_buffer_index << std::endl;
   os << "Pointer to buffer data = " << static_cast<const void *>(d_read_buffer) << std::endl;
}

}
}

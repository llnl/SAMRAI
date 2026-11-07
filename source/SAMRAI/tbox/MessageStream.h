/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Fixed-size message buffer used in interprocessor communication
 *
 ************************************************************************/

#ifndef included_tbox_MessageStream
#define included_tbox_MessageStream

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/AllocatorDatabase.h"
#include "SAMRAI/tbox/Utilities.h"

#ifdef HAVE_UMPIRE
#include "umpire/ResourceManager.hpp"
#include "umpire/TypedAllocator.hpp"
#endif

#include <cstring>
#include <iostream>
#include <vector>

namespace SAMRAI {
namespace tbox {

/*!
 * @brief Class to provide buffers for communication of data.
 *
 * MessageStream provides a message buffer that can hold data of any
 * type.  It is used by communication routines in the Schedule class.
 *
 * TODO: Because this class supports both read and write modes, it has
 * extra data and methods that don't make sense, depending on the
 * mode.  It should be rewritten as two classes, like std::cin and
 * std::cout are.  BTNG.
 *
 * @see Schedule
 */

class MessageStream
{
public:
   /*!
    * @brief Enumeration to identify if a buffer is being used to read or
    * write data.
    */
   enum StreamMode { Read, Write };

   /*!
    * @brief Create a message stream of the specified size and mode
    *
    * @param[in] num_bytes   Number of bytes in the stream.
    *
    * @param[in] mode    MessageStream::Read or MessageStream::Write.
    *
    * @param[in] data_to_read    Data for unpacking, should be num_bytes bytes
    *   long.  This is used when mode == MessageStream::Read, ignored in write
    *   mode.
    *
    * @param[in] deep_copy Whether to make deep copy of data_to_read.
    * The default is to make a deep copy, which is safer but slower
    * than a shallow (pointer) copy.  This is used when mode ==
    * MessageStream::Read, ignored in write mode.  In shallow copy mode,
    * you cannot call growBufferAsNeeded().
    *
    * @param[in] allocator  Optional argument only available when
    * configured with Umpire library.  This allocator will be used for
    * allocations of internal buffers inside this object.
    *
    * @pre num_bytes >= 0
    * @pre mode != Read || data_to_read != 0
    */
   MessageStream(
      const size_t num_bytes,
      const StreamMode mode,
      const void* data_to_read = 0,
      bool deep_copy = true
#ifdef HAVE_UMPIRE
      , umpire::TypedAllocator<char> allocator = 
           AllocatorDatabase::getDatabase()->getInternalHostAllocator()
#endif
   );

#ifdef HAVE_UMPIRE
   /*!
    * @brief Creates a message stream that grows as needed for writing
    * which useds a specified Umpire allocator
    *
    * @param[in] allocator  This allocator will be used for allocations
    * of internal buffers inside this object.
    */
   MessageStream(umpire::TypedAllocator<char> allocator);
#endif

   /*!
    * @brief Default constructor creates a message stream with a
    * buffer that automatically grows as needed, for writing.
    */
   MessageStream();

   /*!
    * Destructor for a message stream.
    */
   ~MessageStream();

   /*!
    * @brief Static method to get amount of message stream space needed to
    * communicate data type indicated by template parameter.
    *
    * IMPORTANT:  All size information given to the message stream should
    * be based on values returned by this method.
    *
    * TODO:  Implementation should be moved out of header?  If we do this,
    * Then we need to create another implementation file to include in this
    * header.  I don't think it's worth it. RDH
    *
    * @return The number of bytes for num_items of type DATA_TYPE.
    *
    * @param[in] num_items
    */
   template<typename DATA_TYPE>
   static size_t getSizeof(size_t num_items = 1)
   {
      return num_items * static_cast<unsigned int>(sizeof(DATA_TYPE));
   }

   /*!
    * @brief Return a pointer to the start of the message buffer.
    */
   const void *
   getBufferStart() const
   {
      if (d_mode == Read) {
         return static_cast<const void *>(d_read_buffer);
      } else {
         return &d_write_buffer[0];
      }
   }

   /*!
    * @brief Return the current size of the buffer in bytes.
    */
   size_t
   getCurrentSize() const
   {
      return d_buffer_index;
   }

   /*!
    * @brief Tell a Write-mode stream to allocate more buffer
    * as needed for data.
    *
    * It is an error to use this method for a Read-mode stream.
    *
    * @pre writeMode()
    */
   void
   growBufferAsNeeded()
   {
      TBOX_ASSERT(writeMode());
      d_grow_as_needed = true;
   }

   /*!
    * @brief Whether a Read-mode MessageStream has reached the end of
    * its data.
    *
    * @pre readMode()
    */
   bool endOfData() const
   {
      TBOX_ASSERT(readMode());
      return d_buffer_index >= d_buffer_size;
   }

   /*!
    * @brief Returns a pointer into the message stream valid for
    * num_bytes.
    *
    * @param[in] num_bytes  Number of bytes requested for window.
    *
    * @pre readMode()
    */
   template<typename DATA_TYPE>
   const DATA_TYPE *
   getReadBuffer(
      size_t num_entries)
   {
      TBOX_ASSERT(readMode());
      const size_t num_bytes = getSizeof<DATA_TYPE>(num_entries);
      TBOX_ASSERT(canCopyOut(num_bytes));
      const DATA_TYPE *buffer =
         reinterpret_cast<const DATA_TYPE *>(&d_read_buffer[getCurrentSize()]);
      d_buffer_index += num_bytes;
      return buffer;
   }

   /*!
    * @brief Returns a pointer into the message stream valid for
    * num_bytes.
    *
    * @param[in] num_bytes  Number of bytes requested for window.
    *
    * @pre writeMode()
    */
   template<typename DATA_TYPE>
   DATA_TYPE *
   getWriteBuffer(
      size_t num_entries)
   {
      TBOX_ASSERT(writeMode());
      const size_t num_bytes = getSizeof<DATA_TYPE>(num_entries);
      if (num_bytes > 0) {
         d_write_buffer.resize(getCurrentSize() + num_bytes);
         d_buffer_size = d_write_buffer.size();
      }
      DATA_TYPE *buffer =
         reinterpret_cast<DATA_TYPE *>(&d_write_buffer[getCurrentSize()]);
      d_buffer_index += num_bytes;
      return buffer;
   }

   /*!
    * @brief Pack a single data item into message stream.
    *
    * @param[in] data  Single item of type DATA_TYPE to be copied
    * into the stream.
    *
    * @pre writeMode()
    */
   template<typename DATA_TYPE>
   MessageStream&
   operator << (
      const DATA_TYPE& data)
   {
      TBOX_ASSERT(writeMode());
      static const size_t nbytes =
         MessageStream::getSizeof<DATA_TYPE>(1);
      copyDataIn(static_cast<const void *>(&data), nbytes);
      return *this;
   }

   /*!
    * @brief Pack an array of data items into message stream.
    *
    * @param[in] data  Pointer to an array of data of type DATA_TYPE
    *                  to be copied into the stream.
    * @param[in] size  Number of items to pack.
    *
    * @pre writeMode()
    */
   template<typename DATA_TYPE>
   void
   pack(
      const DATA_TYPE* data,
      size_t size = 1)
   {
      TBOX_ASSERT(writeMode());
      if (data && (size > 0)) {
         const size_t nbytes = MessageStream::getSizeof<DATA_TYPE>(size);
         copyDataIn(static_cast<const void *>(data), nbytes);
      }
   }

   /*!
    * @brief Pack content of another data stream into this one.
    *
    * @param[in] other  The other data stream.
    *
    * @pre writeMode()
    */
   void
   pack(
      const MessageStream& other)
   {
      TBOX_ASSERT(writeMode());
      if (other.getCurrentSize() > 0) {
         copyDataIn(other.getBufferStart(), other.getCurrentSize());
      }
   }

   /*!
    * @brief Unpack a single data item from message stream.
    *
    * @param[out] data  Single item of type DATA_TYPE that will be
    *                   copied from the stream.
    *
    * @pre readMode()
    */
   template<typename DATA_TYPE>
   MessageStream&
   operator >> (
      DATA_TYPE& data)
   {
      TBOX_ASSERT(readMode());
      static const size_t nbytes =
         MessageStream::getSizeof<DATA_TYPE>(1);
      copyDataOut(static_cast<void *>(&data), nbytes);
      return *this;
   }

   /*!
    * @brief Unpack an array of data items from message stream.
    *
    * @param[out] data  Pointer to an array of data of type DATA_TYPE
    *                   that will receive data copied from
    *                   the stream.
    * @param[out] size  Number of items that will be copied.
    *
    * @pre readMode()
    */
   template<typename DATA_TYPE>
   void
   unpack(
      DATA_TYPE* data,
      size_t size = 1)
   {
      TBOX_ASSERT(readMode());
      if (data) {
         const size_t nbytes = MessageStream::getSizeof<DATA_TYPE>(size);
         copyDataOut(static_cast<void *>(data), nbytes);
      }
   }

   /*!
    * @brief Print out internal object data.
    *
    * @param[out] os  Output stream.
    */
   void
   printClassData(
      std::ostream& os) const;

   /*!
    * @brief Returns true if stream is in read mode.
    */
   bool
   readMode() const
   {
      return d_mode == Read;
   }

   /*!
    * @brief Returns true if stream is in write mode.
    */
   bool
   writeMode() const
   {
      return d_mode == Write;
   }

   /*!
    * @brief Returns true if the buffer is grown as needed in Write mode.
    */
   bool
   growAsNeeded() const
   {
      return d_grow_as_needed;
   }

   /*!
    * @brief Returns true if num_bytes can be copied into the stream.
    */
   bool
   canCopyIn(
      size_t num_bytes) const
   {
      return d_buffer_index + num_bytes <= d_write_buffer.capacity();
   }

   /*!
    * @brief Returns true if num_bytes can be copied out of the stream.
    */
   bool
   canCopyOut(
      size_t num_bytes) const
   {
      return d_buffer_index + num_bytes <= d_buffer_size;
   }

private:
   /*!
    * @brief Copy data into the stream, advancing the stream pointer.
    *
    * @param[in]  input_data
    * @param[in]  num_bytes
    *
    * @pre growAsNeeded() || canCopyIn(num_bytes)
    */
   void copyDataIn(
      const void* input_data,
      const size_t num_bytes)
   {
      if (!growAsNeeded()) {
         TBOX_ASSERT(canCopyIn(num_bytes));
      }
      if (num_bytes > 0) {
         d_write_buffer.insert(d_write_buffer.end(),
            static_cast<const char *>(input_data),
            static_cast<const char *>(input_data) + num_bytes);
         d_buffer_size = d_write_buffer.size();
         d_buffer_index += num_bytes;
      }
   }

   /*!
    * @brief Copy data out of the stream, advancing the stream pointer.
    *
    * @param[in]  output_data
    * @param[in]  num_bytes
    *
    * @pre canCopyOut(num_bytes)
    */
   void copyDataOut(
      void* output_data,
      const size_t num_bytes)
   {
      TBOX_ASSERT(canCopyOut(num_bytes));
      memcpy(output_data, &d_read_buffer[d_buffer_index], num_bytes);
      d_buffer_index += num_bytes;
   }

   MessageStream(
      const MessageStream&);            // not implemented
   MessageStream&
   operator = (
      const MessageStream&);            // not implemented

   /*!
    * @brief  Read/write mode of the stream.
    */
   const StreamMode d_mode;

   /*!
    * The buffer for the streamed data to be written.
    */
#ifdef HAVE_UMPIRE
   std::vector<char, umpire::TypedAllocator<char> > d_write_buffer;
#else
   std::vector<char> d_write_buffer;
#endif

   /*!
    * @brief Pointer to the externally supplied memory to read from in
    * shallow-copy Read mode, or the internal copy of the externally supplied
    * memory to read from in deep-copy Read mode.
    */
   const char* d_read_buffer;

   /*!
    * The allocator used for internal buffers
    */
#ifdef HAVE_UMPIRE
   umpire::TypedAllocator<char> d_allocator;
#endif 

   /*!
    * @brief Number of bytes in the buffer.
    *
    * Equal to d_write_buffer.size() in write mode, size of supplied external
    * buffer size in read mode.
    */
   size_t d_buffer_size;

   /*!
    * Current index into the buffer used when traversing.
    */
   size_t d_buffer_index;

   /*!
    * @brief Whether to grow buffer as needed in a Write-mode stream.
    */
   bool d_grow_as_needed;

   /*!
    * @brief True if d_read_buffer is a deep copy (locally allocated copy) of
    * externally supplied memory.
    */
   bool d_deep_copy_read;

};

}
}

#endif

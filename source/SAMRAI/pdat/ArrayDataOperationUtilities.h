/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated array data looping operations supporting patch data types
 *
 ************************************************************************/

#ifndef included_pdat_ArrayDataOperationUtilities
#define included_pdat_ArrayDataOperationUtilities

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"

namespace SAMRAI {
namespace pdat {

template<class TYPE>
class ArrayData;

/*!
 * @brief Struct ArrayDataOperationUtilities<TYPE, OP> provides
 * generic looping operations for all array-based patch data types.
 * The operations are templated on data type
 * and the operation that will be performed on individual array
 * elements in the innermost loop.
 *
 * @see ArrayData
 */

template<class TYPE, class OP>
class ArrayDataOperationUtilities
{
public:
   /*!
    * Perform operation on a subset of data components of source and
    * destination array data objects and put results in destination array data
    * object.
    *
    * @param dst    Reference to destination array data object.
    * @param src    Const reference to source array data object.
    * @param opbox  Const reference to Box indicating index space region of
    *               operation.
    * @param src_shift  Const reference to IntVector indicating shift required
    *                   to put source index space region into destination index
    *                   space region.
    * @param dst_start_depth  Integer specifying starting depth component of
    *                         operation in destination array.
    * @param src_start_depth  Integer specifying starting depth component of
    *                         operation in source array.
    * @param num_depth  Integer number of depth components on which to perform
    *                   operation.
    * @param op  Const reference to object that performs operations on
    *            individual data array elements.
    *
    * @pre (dst.getDim() == src.getDim()) &&
    *      (dst.getDim() == opbox.getDim()) &&
    *      (dst.getDim() == src_shift.getDim())
    * @pre num_depth >= 0
    * @pre (0 <= dst_start_depth) &&
    *      (dst_start_depth + num_depth <= dst.getDepth())
    * @pre (0 <= src_start_depth) &&
    *      (src_start_depth + num_depth <= src.getDepth())
    */
   static void
   doArrayDataOperationOnBox(
      ArrayData<TYPE>& dst,
      const ArrayData<TYPE>& src,
      const hier::Box& opbox,
      const hier::IntVector& src_shift,
      unsigned int dst_start_depth,
      unsigned int src_start_depth,
      unsigned int num_depth,
      const OP& op);

   /*!
    * Perform operation on all data components of array data object and
    * corresponding buffer data, putting results in either the array data
    * object or buffer.
    *
    * @param arraydata   Const reference to array data object.
    * @param buffer      Const pointer to first element in buffer.
    * @param opbox       Const reference to Box indicating operation region in
    *                    index space of array data object.
    * @param src_is_buffer  Boolean value indicating whether buffer is source
    *                       data for operation; if true results will be placed
    *                       in array data object, otherwise results will go in
    *                       buffer.
    * @param op  Const reference to object that performs operations on
    *            individual data array elements.
    *
    * @pre arraydata.getDim() == opbox.getDim()
    * @pre buffer != 0
    * @pre opbox.isSpatiallyEqual(opbox * arraydata.getBox())
    */
   static void
   doArrayDataBufferOperationOnBox(
      const ArrayData<TYPE>& arraydata,
      const TYPE * buffer,
      const hier::Box& opbox,
      bool src_is_buffer,
      const OP& op);
private:
   // the following are not implemented:
   ArrayDataOperationUtilities();
   ~ArrayDataOperationUtilities();
   ArrayDataOperationUtilities(
      const ArrayDataOperationUtilities&);
   ArrayDataOperationUtilities&
   operator = (
      const ArrayDataOperationUtilities&);
};

}
}

#include "SAMRAI/pdat/ArrayDataOperationUtilities.cpp"

#endif

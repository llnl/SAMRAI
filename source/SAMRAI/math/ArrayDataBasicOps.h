/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Basic templated operations for array data.
 *
 ************************************************************************/

#ifndef included_math_ArrayDataBasicOps
#define included_math_ArrayDataBasicOps

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/ArrayData.h"

namespace SAMRAI {
namespace math {

/**
 * Class ArrayDataBasicOps implements a set of basic operations
 * that apply to numerical data maintained as pdat::ArrayData<TYPE> objects.
 * These operations include simple arithmetic operations as well as min
 * and max, etc.  This class provides a single implementation of these
 * operations that may used to manipulate any of the standard array-based
 * patch data types defined on a patch.  Note that each member function
 * accepts a box argument which specifies the portion of the array data
 * on which the associated operation is performed.   The actual index
 * region on which the operation occurs is the intersection of this box
 * and the boxes of all the pdat::ArrayData<TYPE> objects involved.
 *
 * These operations typically apply only to the numerical standard built-in
 * types, such as double, float, and int, and the complex type (which may or
 * may not be a built-in type depending on the C++ compiler).  Thus, this
 * templated class should only be used to instantiate objects with those
 * types as the template parameter.  Those operations whose implementations
 * depend of the data type are specialized for each numerical type.  To use
 * this class with other standard types or user-defined types (which may or
 * may not make sense), the member functions must be specialized so that the
 * correct operations are performed.
 *
 * @see pdat::ArrayData
 */

template<class TYPE>
class ArrayDataBasicOps
{
public:
   /**
    * Empty constructor and destructor.
    */
   ArrayDataBasicOps();

   ~ArrayDataBasicOps();

   /**
    * Set dst = alpha * src, elementwise.
    *
    * @pre (dst.getDim() == src.getDim()) && (dst.getDim() == box.getDim())
    * @pre (alpha == tbox::MathUtilities<TYPE>::getZero()) ||
    *      (alpha == tbox::MathUtilities<TYPE>::getOne()) ||
    *      (dst.getDepth() == src.getDepth())
    */
   void
   scale(
      pdat::ArrayData<TYPE>& dst,
      const TYPE& alpha,
      const pdat::ArrayData<TYPE>& src,
      const hier::Box& box) const;

   /**
    * Set dst = src + alpha, elementwise.
    *
    * @pre (dst.getDim() == src.getDim()) && (dst.getDim() == box.getDim())
    * @pre (alpha == tbox::MathUtilities<TYPE>::getZero()) ||
    *      (dst.getDepth() == src.getDepth())
    */
   void
   addScalar(
      pdat::ArrayData<TYPE>& dst,
      const pdat::ArrayData<TYPE>& src,
      const TYPE& alpha,
      const hier::Box& box) const;

   /**
    * Set dst = src1 + src2, elementwise.
    *
    * @pre (dst.getDim() == src1.getDim()) &&
    *      (dst.getDim() == src2.getDim()) && (dst.getDim() == box.getDim())
    * @pre (dst.getDepth() == src1.getDepth()) &&
    *      (dst.getDepth() == src2.getDepth())
    */
   void
   add(
      pdat::ArrayData<TYPE>& dst,
      const pdat::ArrayData<TYPE>& src1,
      const pdat::ArrayData<TYPE>& src2,
      const hier::Box& box) const;

   /**
    * Set dst = src1 - src2, elementwise.
    *
    * @pre (dst.getDim() == src1.getDim()) &&
    *      (dst.getDim() == src2.getDim()) && (dst.getDim() == box.getDim())
    * @pre (dst.getDepth() == src1.getDepth()) &&
    *      (dst.getDepth() == src2.getDepth())
    */
   void
   subtract(
      pdat::ArrayData<TYPE>& dst,
      const pdat::ArrayData<TYPE>& src1,
      const pdat::ArrayData<TYPE>& src2,
      const hier::Box& box) const;

   /**
    * Set dst = src1 * src2, elementwise.
    *
    * @pre (dst.getDim() == src1.getDim()) &&
    *      (dst.getDim() == src2.getDim()) && (dst.getDim() == box.getDim())
    * @pre (dst.getDepth() == src1.getDepth()) &&
    *      (dst.getDepth() == src2.getDepth())
    */
   void
   multiply(
      pdat::ArrayData<TYPE>& dst,
      const pdat::ArrayData<TYPE>& src1,
      const pdat::ArrayData<TYPE>& src2,
      const hier::Box& box) const;

   /**
    * Set dst = src1 / src2, elementwise.  No check for division by zero.
    *
    * @pre (dst.getDim() == src1.getDim()) &&
    *      (dst.getDim() == src2.getDim()) && (dst.getDim() == box.getDim())
    * @pre (dst.getDepth() == src1.getDepth()) &&
    *      (dst.getDepth() == src2.getDepth())
    */
   void
   divide(
      pdat::ArrayData<TYPE>& dst,
      const pdat::ArrayData<TYPE>& src1,
      const pdat::ArrayData<TYPE>& src2,
      const hier::Box& box) const;

   /**
    * Set dst = 1 / src, elementwise.  No check for division by zero.
    *
    * @pre (dst.getDim() == src.getDim()) && (dst.getDim() == box.getDim())
    * @pre dst.getDepth() == src.getDepth()
    */
   void
   reciprocal(
      pdat::ArrayData<TYPE>& dst,
      const pdat::ArrayData<TYPE>& src,
      const hier::Box& box) const;

   /**
    * Set dst = alpha * src1 + beta * src2, elementwise.
    *
    * @pre (dst.getDim() == src1.getDim()) &&
    *      (dst.getDim() == src2.getDim()) && (dst.getDim() == box.getDim())
    * @pre (dst.getDepth() == src1.getDepth()) &&
    *      (dst.getDepth() == src2.getDepth())
    */
   void
   linearSum(
      pdat::ArrayData<TYPE>& dst,
      const TYPE& alpha,
      const pdat::ArrayData<TYPE>& src1,
      const TYPE& beta,
      const pdat::ArrayData<TYPE>& src2,
      const hier::Box& box) const;

   /**
    * Set dst = alpha * src1 + src2, elementwise.
    *
    * @pre (dst.getDim() == src1.getDim()) &&
    *      (dst.getDim() == src2.getDim()) && (dst.getDim() == box.getDim())
    * @pre (alpha == tbox::MathUtilities<TYPE>::getZero()) ||
    *      (alpha == tbox::MathUtilities<TYPE>::getOne()) ||
    *      (alpha == -tbox::MathUtilities<TYPE>::getOne()) ||
    *      ((dst.getDepth() == src1.getDepth()) &&
    *       (dst.getDepth() == src2.getDepth()))
    */
   void
   axpy(
      pdat::ArrayData<TYPE>& dst,
      const TYPE& alpha,
      const pdat::ArrayData<TYPE>& src1,
      const pdat::ArrayData<TYPE>& src2,
      const hier::Box& box) const;

   /**
    * Set dst = alpha * src1 - src2, elementwise.
    *
    * @pre (dst.getDim() == src1.getDim()) &&
    *      (dst.getDim() == src2.getDim()) && (dst.getDim() == box.getDim())
    * @pre (alpha == tbox::MathUtilities<TYPE>::getZero()) ||
    *      (alpha == tbox::MathUtilities<TYPE>::getOne()) ||
    *      ((dst.getDepth() == src1.getDepth()) &&
    *       (dst.getDepth() == src2.getDepth()))
    */
   void
   axmy(
      pdat::ArrayData<TYPE>& dst,
      const TYPE& alpha,
      const pdat::ArrayData<TYPE>& src1,
      const pdat::ArrayData<TYPE>& src2,
      const hier::Box& box) const;

   /**
    * Return the minimum array data entry.  If data is complex, return the
    * array data entry with the minimum norm.
    *
    * @pre data.getDim() == box.getDim()
    */
   TYPE
   min(
      const pdat::ArrayData<TYPE>& data,
      const hier::Box& box) const;

   /**
    * Return the maximum array data entry.  If data is complex, return the
    * array data entry with the maximum norm.
    *
    * @pre data.getDim() == box.getDim()
    */
   TYPE
   max(
      const pdat::ArrayData<TYPE>& data,
      const hier::Box& box) const;

   /**
    * Set dst to random values.  If the data is int, each element of dst
    * is set as dst = mrand48(). If the data is double or float, each
    * element of dst is set as dst = width * drand48() + low.  If the
    * data is complex, each element of dst is set as dst = dcomplex(rval, ival),
    * where rval = real(width) * drand48() + real(low), and
    * ival = imag(width) * drand48() + imag(low).
    *
    * @pre dst.getDim() == box.getDim()
    */
   void
   setRandomValues(
      pdat::ArrayData<TYPE>& dst,
      const TYPE& width,
      const TYPE& low,
      const hier::Box& box) const;

private:
   // The following are not implemented:
   ArrayDataBasicOps(
      const ArrayDataBasicOps&);
   ArrayDataBasicOps&
   operator = (
      const ArrayDataBasicOps&);

};

}
}

#include "SAMRAI/math/ArrayDataBasicOps.cpp"

#endif

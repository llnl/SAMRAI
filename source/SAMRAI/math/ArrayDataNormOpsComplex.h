/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Norm operations for complex data arrays.
 *
 ************************************************************************/

#ifndef included_math_ArrayDataNormOpsComplex
#define included_math_ArrayDataNormOpsComplex

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/Utilities.h"

#include <cmath>

namespace SAMRAI {
namespace math {

/**
 * Class ArrayDataNormOpsComplex provides a set of common norm
 * operations that may be applied to arrays of complex data values
 * maintained as pdat::ArrayData<TYPE> objects.  The intent of this class is to
 * provide a single implementation of these operations as they are needed
 * by objects that perform these operations on the standard array-based patch
 * data types (i.e., cell-centered, face-centered, node-centered).  Each of
 * the norm operations is implemented in two different ways.  The choice of
 * operation is based on whether control volume information is to be used to
 * weight the contribution of each data entry to the norm calculation.  The
 * use of control volumes is important when these operations are used in
 * vector kernels where the data resides over multiple levels in an AMR
 * hierarchy.  Note also that each operation will be performed on the
 * intersection of the box in the function argument list and the boxes
 * associated with all pdat::ArrayData<TYPE> objects.
 *
 * Note that a similar set of norm operations is implemented for real array
 * data (double and float) in the class ArrayDataNormOpsReal.
 *
 * @see pdat::ArrayData
 */

class ArrayDataNormOpsComplex
{
public:
   /**
    * Empty constructor and destructor.
    */
   ArrayDataNormOpsComplex();

   ~ArrayDataNormOpsComplex();

   /**
    * Set destination component to norm of source component.  That is,
    * each destination entry is set to
    * \f$d_i = \sqrt{ {real(s_i)}^2 + {imag(s_i)}^2 }\f$.
    *
    * @pre (dst.getDim() == src.getDim()) && (dst.getDim() == box.getDim())
    * @pre dst.getDepth() == src.getDepth()
    */
   void
   abs(
      pdat::ArrayData<double>& dst,
      const pdat::ArrayData<dcomplex>& src,
      const hier::Box& box) const;

   /**
    * Return sum of entries in control volume array.
    *
    * @pre (data.getDim() == cvol.getDim()) && (data.getDim() == box.getDim())
    */
   double
   sumControlVolumes(
      const pdat::ArrayData<dcomplex>& data,
      const pdat::ArrayData<double>& cvol,
      const hier::Box& box) const;

   /**
    * Return discrete \f$L_1\f$-norm of the data using the control volume to
    * weight the contribution of each data entry to the sum.  That is, the
    * return value is the sum \f$\sum_i ( \sqrt{data_i * \bar{data_i}} cvol_i )\f$.
    *
    * @pre (data.getDim() == cvol.getDim()) && (data.getDim() == box.getDim())
    */
   double
   L1NormWithControlVolume(
      const pdat::ArrayData<dcomplex>& data,
      const pdat::ArrayData<double>& cvol,
      const hier::Box& box) const;

   /**
    * Return discrete \f$L_1\f$-norm of the data.  That is, the return value is
    * the sum \f$\sum_i ( \sqrt{data_i * \bar{data_i}} )\f$.
    *
    * @pre (data.getDim() == box.getDim()
    */
   double
   L1Norm(
      const pdat::ArrayData<dcomplex>& data,
      const hier::Box& box) const;

   /**
    * Return discrete \f$L_2\f$-norm of the data using the control volume to
    * weight the contribution of each data entry to the sum.  That is, the
    * return value is the sum \f$\sqrt{ \sum_i (
    * data_i * \bar{data_i} cvol_i ) }\f$.
    *
    * @pre (data.getDim() == cvol.getDim()) && (data.getDim() == box.getDim())
    */
   double
   L2NormWithControlVolume(
      const pdat::ArrayData<dcomplex>& data,
      const pdat::ArrayData<double>& cvol,
      const hier::Box& box) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY3(data, cvol, box);
      return sqrt(real(dotWithControlVolume(data, data, cvol, box)));
   }

   /**
    * Return discrete \f$L_2\f$-norm of the data using the control volume to
    * weight the contribution of each data entry to the sum.  That is, the
    * return value is the sum \f$\sqrt{ \sum_i ( data_i * \bar{data_i} ) }\f$.
    *
    * @pre data.getDim() == box.getDim()
    */
   double
   L2Norm(
      const pdat::ArrayData<dcomplex>& data,
      const hier::Box& box) const
   {
      TBOX_ASSERT_OBJDIM_EQUALITY2(data, box);
      return sqrt(real(dot(data, data, box)));
   }

   /**
    * Return discrete weighted \f$L_2\f$-norm of the data using the control
    * volume to weight the contribution of the data and weight entries to
    * the sum.  That is, the return value is the sum \f$\sqrt{ \sum_i (
    * (data_i * wgt_i) * \bar{(data_i * wgt_i)} cvol_i ) }\f$.
    *
    * @pre (data.getDim() == wgt.getDim()) &&
    *      (data.getDim() == cvol.getDim()) && (data.getDim() == box.getDim())
    * @pre data.getDepth() == wgt.getDepth()
    */
   double
   weightedL2NormWithControlVolume(
      const pdat::ArrayData<dcomplex>& data,
      const pdat::ArrayData<dcomplex>& wgt,
      const pdat::ArrayData<double>& cvol,
      const hier::Box& box) const;

   /**
    * Return discrete weighted \f$L_2\f$-norm of the data.  That is, the return
    * value is the sum \f$\sqrt{ \sum_i ( (data_i * wgt_i)
    * \bar{(data_i * wgt_i)} cvol_i ) }\f$.
    *
    * @pre (data.getDim() == wgt.getDim()) && (data.getDim() == box.getDim())
    * @pre data.getDepth() == wgt.getDepth()
    */
   double
   weightedL2Norm(
      const pdat::ArrayData<dcomplex>& data,
      const pdat::ArrayData<dcomplex>& wgt,
      const hier::Box& box) const;

   /**
    * Return the \f$\max\f$-norm of the data using the control volume to weight
    * the contribution of each data entry to the maximum.  That is, the return
    * value is \f$\max_i ( \sqrt{data_i * \bar{data_i}} )\f$, where the max is
    * over the data elements where \f$cvol_i > 0\f$.
    *
    * @pre (data.getDim() == cvol.getDim()) && (data.getDim() == box.getDim())
    */
   double
   maxNormWithControlVolume(
      const pdat::ArrayData<dcomplex>& data,
      const pdat::ArrayData<double>& cvol,
      const hier::Box& box) const;

   /**
    * Return the \f$\max\f$-norm of the data.  That is, the return value is
    * \f$\max_i ( \sqrt{data_i * \bar{data_i}} )\f$.
    *
    * @pre data.getDim() == box.getDim()
    */
   double
   maxNorm(
      const pdat::ArrayData<dcomplex>& data,
      const hier::Box& box) const;

   /**
    * Return the dot product of the two data arrays using the control volume
    * to weight the contribution of each product to the sum.  That is, the
    * return value is the sum \f$\sum_i ( data1_i * \bar{data2_i} * cvol_i )\f$.
    *
    * @pre (data1.getDim() == data2.getDim()) &&
    *      (data1.getDim() == cvol.getDim()) &&
    *      (data1.getDim() == box.getDim())
    * @pre data1.getDepth() == data2.getDepth()
    */
   dcomplex
   dotWithControlVolume(
      const pdat::ArrayData<dcomplex>& data1,
      const pdat::ArrayData<dcomplex>& data2,
      const pdat::ArrayData<double>& cvol,
      const hier::Box& box) const;

   /**
    * Return the dot product of the two data arrays.  That is, the
    * return value is the sum \f$\sum_i ( data1_i * \bar{data2_i} )\f$.
    *
    * @pre (data1.getDim() == data2.getDim()) &&
    *      (data1.getDim() == box.getDim())
    * @pre data1.getDepth() == data2.getDepth()
    */
   dcomplex
   dot(
      const pdat::ArrayData<dcomplex>& data1,
      const pdat::ArrayData<dcomplex>& data2,
      const hier::Box& box) const;

   /**
    * Return the integral of the function based on the data array.
    * The return value is the sum \f$\sum_i ( data_i * vol_i )\f$.
    *
    * @pre (data.getDim() == vol.getDim()) && (data.getDim() == box.getDim())
    */
   dcomplex
   integral(
      const pdat::ArrayData<dcomplex>& data,
      const pdat::ArrayData<double>& vol,
      const hier::Box& box) const;

private:
   // The following are not implemented:
   ArrayDataNormOpsComplex(
      const ArrayDataNormOpsComplex&);
   ArrayDataNormOpsComplex&
   operator = (
      const ArrayDataNormOpsComplex&);
};

}
}

#endif

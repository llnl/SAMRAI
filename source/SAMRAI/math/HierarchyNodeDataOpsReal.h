/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Templated operations for real node data on multiple levels.
 *
 ************************************************************************/

#ifndef included_math_HierarchyNodeDataOpsReal
#define included_math_HierarchyNodeDataOpsReal

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/math/HierarchyDataOpsReal.h"
#include "SAMRAI/math/PatchNodeDataOpsReal.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/PatchHierarchy.h"

#include <iostream>
#include <memory>

namespace SAMRAI {
namespace math {

/**
 * Class HierarchyNodeDataOpsReal provides a collection of operations
 * that manipulate real (double and float) node-centered patch data components
 * over multiple levels in an AMR hierarchy.  It is derived from the abstract
 * base class HierarchyDataOpsReal which defines the interface to
 * similar operations for node-centered, face-centered, node-centered patch
 * data objects where the data is either double or float.  The operations
 * include basic arithmetic, norms and ordering, and assorted miscellaneous
 * operations.  On each patch, the operations are performed by the
 * PatchNodeDataOpsReal data member.
 *
 * The patch hierarchy and set of levels within that hierarcy over which the
 * operations will be performed are set in the constructor.  However, note
 * that the constructor accepts default arguments for the coarsest and finest
 * level numbers.  If the level numbers are not specified when calling the
 * constructor the levels which exist in the hierarchy will be assumed in
 * all operations.  The hierarchy and levels may be changed at any time using
 * the proper member functions.
 *
 * Note that, when it makes sense, an operation accept a boolean argument
 * which indicates whether the operation should be performed on all of the
 * data or just those data elements corresponding to the patch interiors.
 * If no boolean argument is provided, the default behavior is to treat only
 * the patch interiors.  Also, a similar set of operations for complex and
 * integer node-centered data is provided in the classes
 * HierarchyNodeDataOpsComplex and HierarchyNodeDataOpsInteger,
 * respectively.
 *
 * @see PatchNodeDataOpsReal
 */

template<class TYPE>
class HierarchyNodeDataOpsReal:public HierarchyDataOpsReal<TYPE>
{
public:
   /**
    * The constructor for the HierarchyNodeDataOpsReal class sets
    * the default patch hierarchy and coarsest and finest patch levels
    * in that hierarchy over which operations will be performed.  The
    * hierarchy and operations may be reset using the member fuctions
    * setPatchHierarchy() and resetLevels() below.  If no level number
    * arguments are given here, the levels over which the operations will
    * be performed are those already existing in the hierarchy.  If the
    * hierarchy level configuration changes, the operations must be explicitly
    * reset by calling the resetLevels() function.
    *
    * @pre hierarchy
    */
   explicit HierarchyNodeDataOpsReal(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int coarsest_level = -1,
      const int finest_level = -1);

   /**
    * Virtual destructor for the HierarchyNodeDataOpsReal class.
    */
   virtual ~HierarchyNodeDataOpsReal();

   /**
    * Reset patch hierarchy over which operations occur.
    *
    * @pre hierarchy
    */
   void
   setPatchHierarchy(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy);

   /**
    * Reset range of patch levels over which operations occur.
    *
    * @pre getPatchHierarchy()
    * @pre (coarsest_level >= 0) && (finest_level >= coarsest_level) &&
    *      (finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   resetLevels(
      const int coarsest_level,
      const int finest_level);

   /**
    * Return const pointer to patch hierarchy associated with operations.
    */
   const std::shared_ptr<hier::PatchHierarchy>
   getPatchHierarchy() const;

   /**
    * Copy source data to destination data.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   copyData(
      const int dst_id,
      const int src_id,
      const bool interior_only = true) const;

   /**
    * Swap data pointers (i.e., storage) between two data components.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   swapData(
      const int data1_id,
      const int data2_id) const;

   /**
    * Print data over multiple levels to specified output stream.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   printData(
      const int data_id,
      std::ostream& s,
      const bool interior_only = true) const;

   /**
    * Set data component to given scalar.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   setToScalar(
      const int data_id,
      const TYPE& alpha,
      const bool interior_only = true) const;

   /**
    * Set destination to source multiplied by given scalar, pointwise.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   scale(
      const int dst_id,
      const TYPE& alpha,
      const int src_id,
      const bool interior_only = true) const;

   /**
    * Add scalar to each entry in source data and set destination to result.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   addScalar(
      const int dst_id,
      const int src_id,
      const TYPE& alpha,
      const bool interior_only = true) const;

   /**
    * Set destination to sum of two source components, pointwise.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   add(
      const int dst_id,
      const int src1_id,
      const int src2_id,
      const bool interior_only = true) const;

   /**
    * Subtract second source component from first source component pointwise
    * and set destination data component to result.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   subtract(
      const int dst_id,
      const int src1_id,
      const int src2_id,
      const bool interior_only = true) const;

   /**
    * Set destination component to product of two source components, pointwise.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   multiply(
      const int dst_id,
      const int src1_id,
      const int src2_id,
      const bool interior_only = true) const;

   /**
    * Divide first data component by second source component pointwise
    * and set destination data component to result.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   divide(
      const int dst_id,
      const int src1_id,
      const int src2_id,
      const bool interior_only = true) const;

   /**
    * Set each entry of destination component to reciprocal of corresponding
    * source data component entry.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   reciprocal(
      const int dst_id,
      const int src_id,
      const bool interior_only = true) const;

   /**
    * Set \f$d = \alpha s_1 + \beta s_2\f$, where \f$d\f$ is the destination patch
    * data component and \f$s_1, s_2\f$ are the first and second source components,
    * respectively.  Here \f$\alpha, \beta\f$ are scalar values.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   linearSum(
      const int dst_id,
      const TYPE& alpha,
      const int src1_id,
      const TYPE& beta,
      const int src2_id,
      const bool interior_only = true) const;

   /**
    * Set \f$d = \alpha s_1 + s_2\f$, where \f$d\f$ is the destination patch data
    * component and \f$s_1, s_2\f$ are the first and second source components,
    * respectively.  Here \f$\alpha\f$ is a scalar.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   axpy(
      const int dst_id,
      const TYPE& alpha,
      const int src1_id,
      const int src2_id,
      const bool interior_only = true) const;

   /**
    * Set \f$d = \alpha s_1 - s_2\f$, where \f$d\f$ is the destination patch data
    * component and \f$s_1, s_2\f$ are the first and second source components,
    * respectively.  Here \f$\alpha\f$ is a scalar.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   axmy(
      const int dst_id,
      const TYPE& alpha,
      const int src1_id,
      const int src2_id,
      const bool interior_only = true) const;

   /**
    * Set destination data to absolute value of source data, pointwise.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   abs(
      const int dst_id,
      const int src_id,
      const bool interior_only = true) const;

   /**
    * Return minimum data value over all patches in the collection of levels.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   TYPE
   min(
      const int data_id,
      const bool interior_only = true) const;

   /**
    * Return maximum data value over all patches in the collection of levels.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   TYPE
   max(
      const int data_id,
      const bool interior_only = true) const;

   /**
    * Set data entries to random values.  See the operations in the
    * array data operation classes for each data type for details on the
    * generation of the random values.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   setRandomValues(
      const int data_id,
      const TYPE& width,
      const TYPE& low,
      const bool interior_only = true) const;

   /**
    * Return the total number of data values for the component on the set
    * of hierarchy levels.  If the boolean argument is true, the number of
    * elements will be summed over patch interiors in a unique way which
    * avoids multiple counting of redundant values (recall the definition
    * of node points on a patch interior).  If the boolean argument is false,
    * all elements will be counted (including ghost values) over all patches.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   size_t
   numberOfEntries(
      const int data_id,
      const bool interior_only = true) const;

   /**
    * Return sum of the control volumes associated with the data component.
    * Note that if the control volumes are set properly, this is equivalent to
    * integrating a data component containing all ones over the collection of
    * hierarchy levels.
    *
    * @pre vol_id >= 0
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   double
   sumControlVolumes(
      const int data_id,
      const int vol_id) const;

   /**
    * Return discrete \f$L_1\f$-norm of the data using the control volume to
    * weight the contribution of each data entry to the sum.  That is, the
    * return value is the sum \f$\sum_i ( \| data_i \| cvol_i )\f$.  If the
    * control volume is not defined (vol_id < 0), the return value is
    * \f$\sum_i ( \| data_i \| )\f$.
    * If local_only is true, the global reduction is not
    * performed (thus each process will get only local results).
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   double
   L1Norm(
      const int data_id,
      const int vol_id = -1,
      bool local_only = false) const;

   /**
    * Return discrete \f$L_2\f$-norm of the data using the control volume to
    * weight the contribution of each data entry to the sum.  That is, the
    * return value is the sum \f$\sqrt{ \sum_i ( (data_i)^2 cvol_i ) }\f$.
    * If the control volume is  not defined (vol_id < 0), the return value
    * is \f$\sqrt{ \sum_i ( (data_i)^2 cvol_i ) }\f$.
    * If local_only is true, the global reduction is not
    * performed (thus each process will get only local results).
    */
   double
   L2Norm(
      const int data_id,
      const int vol_id = -1,
      bool local_only = false) const;

   /**
    * Return discrete weighted \f$L_2\f$-norm of the data using the control
    * volume to weight the contribution of the data and weight entries to
    * the sum.  That is, the return value is the sum \f$\sqrt{ \sum_i (
    * (data_i * weight_i)^2 cvol_i ) }\f$.  If the control volume is not defined
    * (vol_id < 0), the return value is
    * \f$\sqrt{ \sum_i ( (data_i * weight_i)^2 ) }\f$.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   double
   weightedL2Norm(
      const int data_id,
      const int weight_id,
      const int vol_id = -1) const;

   /**
    * Return discrete root mean squared norm of the data.  If the control
    * volume is defined, the return value is the \f$L_2\f$-norm divided by
    * the square root of the sum of the control volumes.  Otherwise, the
    * return value is the \f$L_2\f$-norm divided by the square root of the
    * number of data entries.
    */
   double
   RMSNorm(
      const int data_id,
      const int vol_id = -1) const;

   /**
    * Return discrete weighted root mean squared norm of the data.  If the
    * control volume is define, the return value is the weighted \f$L_2\f$-norm
    * divided by the square root of the sum of the control volumes.  Otherwise,
    * the return value is the weighted \f$L_2\f$-norm divided by the square root
    * of the number of data entries.
    */
   double
   weightedRMSNorm(
      const int data_id,
      const int weight_id,
      const int vol_id = -1) const;

   /**
    * Return the \f$\max\f$-norm of the data using the control volume to weight
    * the contribution of each data entry to the maximum.  That is, the return
    * value is \f$\max_i ( \| data_i \| )\f$, where the max is over the data
    * elements where \f$cvol_i > 0\f$.  If the control volume is undefined
    * (vol_id < 0), it is ignored during the computation of the maximum.
    * If local_only is true, the global reduction is not
    * performed (thus each process will get only local results).
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   double
   maxNorm(
      const int data_id,
      const int vol_id = -1,
      bool local_only = false) const;

   /**
    * Return the dot product of the two data arrays using the control volume
    * to weight the contribution of each product to the sum.  That is, the
    * return value is the sum \f$\sum_i ( data1_i * data2_i * cvol_i )\f$.
    * If the control volume is undefined (vol_id < 0), it is ignored during
    * the summation.  If local_only is true, the global reduction is not
    * performed (thus each process will get only local results).
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   TYPE
   dot(
      const int data1_id,
      const int data2_id,
      const int vol_id = -1,
      bool local_only = false) const;

   /**
    * Return the integral of the function represented by the data array.
    * The return value is the sum \f$\sum_i ( data_i * vol_i )\f$.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   TYPE
   integral(
      const int data_id,
      const int vol_id) const;

   /**
    * Return 1 if \f$\|data2_i\| > 0\f$ and \f$data1_i * data2_i \leq 0\f$, for
    * any \f$i\f$ in the set of patch data indices, where \f$cvol_i > 0\f$.  Otherwise,
    * return 0.  If the control volume is undefined (vol_id < 0), all values
    * on the patch interiors are considered.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   int
   computeConstrProdPos(
      const int data1_id,
      const int data2_id,
      const int vol_id = -1) const;

   /**
    * Wherever \f$cvol_i > 0\f$ in the set of patch data indices, set \f$dst_i = 1\f$
    * if \f$\|src_i\| > \alpha\f$, and \f$dst_i = 0\f$ otherwise.  If the control
    * volume is undefined (vol_id < 0), all values on the patch interiors
    * are considered.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   compareToScalar(
      const int dst_id,
      const int src_id,
      const TYPE& alpha,
      const int vol_id = -1) const;

   /**
    * Wherever \f$cvol_i > 0\f$ in the set of patch data indices, set
    * \f$dst_i = 1/src_i\f$ if \f$src_i \neq 0\f$, and \f$dst_i = 0\f$ otherwise.
    * If \f$dst_i = 0\f$ anywhere, 0 is the return value.  Otherwise 1 is returned.
    * If the control volume is undefined (vol_id < 0), all values on the
    * patch interiors are considered.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   int
   testReciprocal(
      const int dst_id,
      const int src_id,
      const int vol_id = -1) const;

   /*!
    * @brief Compute max of "conditional" quotients of two arrays.
    *
    * Return the maximum of pointwise "conditional" quotients of the numerator
    * and denominator.
    *
    * The "conditional" quotient is defined as |numerator/denominator|
    * where the denominator is nonzero.  Otherwise, it is defined as
    * |numerator|.
    *
    * @b Note: This method is currently intended to support the
    * PETSc-2.1.6 vector wrapper only.  Please do not use it!
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   TYPE
   maxPointwiseDivide(
      const int numer_id,
      const int denom_id,
      bool local_only = false) const;

   /*!
    * @brief Compute min of quotients of two arrays.
    *
    * Return the minimum of pointwise quotients of the numerator
    * and denominator.
    *
    * The quotient is defined as (numerator/denominator)
    * where the denominator is nonzero.  When the denominator is zero, the
    * entry is skipped.  If the denominator is always zero, the value of
    * tbox::IEEE::getFLT_MAX() is returned (see @ref SAMRAI::tbox::IEEE).
    *
    * @b Note: This method is currently intended to support the
    * SUNDIALS vector wrapper only.  Please do not use it!
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   TYPE
   minPointwiseDivide(
      const int numer_id,
      const int denom_id,
      bool local_only = false) const;

   /*!
    * @brief Get the length of the data with the given ID
    *
    * @param data_id         Patch data ID
    * @param interior_only   If true, only count the elements in interior
    *                        of patches.
    * @return The length (number of elements in the underlying data)
    */
   int64_t getLength(
      const int data_id,
      const bool interior_only = true) const;
 
private:
   // The following are not implemented
   HierarchyNodeDataOpsReal(
      const HierarchyNodeDataOpsReal&);
   HierarchyNodeDataOpsReal&
   operator = (
      const HierarchyNodeDataOpsReal&);

   std::shared_ptr<hier::PatchHierarchy> d_hierarchy;
   int d_coarsest_level;
   int d_finest_level;
   std::vector<std::vector<hier::BoxContainer> > d_nonoverlapping_node_boxes;

   PatchNodeDataOpsReal<TYPE> d_patch_ops;

};

}
}

#include "SAMRAI/math/HierarchyNodeDataOpsReal.cpp"

#endif

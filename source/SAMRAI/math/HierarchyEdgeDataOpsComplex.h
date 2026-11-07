/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Operations for complex edge data on multiple levels.
 *
 ************************************************************************/

#ifndef included_math_HierarchyEdgeDataOpsComplex
#define included_math_HierarchyEdgeDataOpsComplex

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/math/HierarchyDataOpsComplex.h"
#include "SAMRAI/math/PatchEdgeDataOpsComplex.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/tbox/Complex.h"

#include <iostream>
#include <vector>
#include <memory>

namespace SAMRAI {
namespace math {

/**
 * Class HierarchyEdgeDataOpsComplex provides a collection of
 * operations that manipulate complex edge-centered patch data components over
 * multiple levels in an AMR hierarchy.  It is derived from the abstract
 * base class HierarchyDataOpsComplex which defines the interface to
 * similar operations for edge-centered, edge-centered, edge-centered patch
 * data objects where the data is complex.   The operations include basic
 * arithmetic and norms.  On each patch, the operations are performed by the
 * PatchEdgeDataOpsComplex data member.
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
 * the patch interiors.  Also, a similar set of operations for real (double
 * and float) and integer edge-centered data is provided in the classes
 * HierarchyEdgeDataOpsReal and HierarchyEdgeDataOpsInteger,
 * respectively.
 *
 * @see PatchEdgeDataOpsComplex
 */

class HierarchyEdgeDataOpsComplex:public HierarchyDataOpsComplex
{
public:
   /**
    * The constructor for the HierarchyEdgeDataOpsComplex class sets
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
   explicit HierarchyEdgeDataOpsComplex(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int coarsest_level = -1,
      const int finest_level = -1);

   /**
    * Virtual destructor for the HierarchyEdgeDataOpsComplex class.
    */
   virtual ~HierarchyEdgeDataOpsComplex();

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
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
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
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
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
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
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
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   setToScalar(
      const int data_id,
      const dcomplex& alpha,
      const bool interior_only = true) const;

   /**
    * Set destination to source multiplied by given scalar, pointwise.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   scale(
      const int dst_id,
      const dcomplex& alpha,
      const int src_id,
      const bool interior_only = true) const;

   /**
    * Add scalar to each entry in source data and set destination to result.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   addScalar(
      const int dst_id,
      const int src_id,
      const dcomplex& alpha,
      const bool interior_only = true) const;

   /**
    * Set destination to sum of two source components, pointwise.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
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
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
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
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
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
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
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
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
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
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   linearSum(
      const int dst_id,
      const dcomplex& alpha,
      const int src1_id,
      const dcomplex& beta,
      const int src2_id,
      const bool interior_only = true) const;

   /**
    * Set \f$d = \alpha s_1 + s_2\f$, where \f$d\f$ is the destination patch data
    * component and \f$s_1, s_2\f$ are the first and second source components,
    * respectively.  Here \f$\alpha\f$ is a scalar.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   axpy(
      const int dst_id,
      const dcomplex& alpha,
      const int src1_id,
      const int src2_id,
      const bool interior_only = true) const;

   /**
    * Set \f$d = \alpha s_1 - s_2\f$, where \f$d\f$ is the destination patch data
    * component and \f$s_1, s_2\f$ are the first and second source components,
    * respectively.  Here \f$\alpha\f$ is a scalar.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   axmy(
      const int dst_id,
      const dcomplex& alpha,
      const int src1_id,
      const int src2_id,
      const bool interior_only = true) const;

   /**
    * Set destination data to absolute value of source data, pointwise.
    * Note that the source data must be dcomplex and the destination must
    * be double.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   abs(
      const int dst_id,
      const int src_id,
      const bool interior_only = true) const;

   /**
    * Set data entries to random values.  See the operations in the
    * array data operation classes for details on the generation of
    * the random values.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   setRandomValues(
      const int data_id,
      const dcomplex& width,
      const dcomplex& low,
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
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
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
    *
    * @pre vol_id >= 0
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    * hierarchy levels.
    */
   double
   sumControlVolumes(
      const int data_id,
      const int vol_id) const;

   /**
    * Return discrete \f$L_1\f$-norm of the data using the control volume to
    * weight the contribution of each data entry to the sum.  That is, the
    * return value is the sum \f$\sum_i ( \sqrt{data_i * \bar{data_i}}*cvol_i )\f$.
    * If the control volume is undefined (vol_id < 0), the
    * return value is \f$\sum_i ( \sqrt{data_i * \bar{data_i}} )\f$.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   double
   L1Norm(
      const int data_id,
      const int vol_id = -1) const;

   /**
    * Return discrete \f$L_2\f$-norm of the data using the control volume to
    * weight the contribution of each data entry to the sum.  That is, the
    * return value is the sum
    * \f$\sqrt{ \sum_i ( data_i * \bar{data_i} cvol_i ) }\f$.
    * If the control volume is undefined (vol_id < 0), the return value is
    * \f$\sqrt{ \sum_i ( data_i * \bar{data_i} ) }\f$.
    */
   double
   L2Norm(
      const int data_id,
      const int vol_id = -1) const;

   /**
    * Return discrete weighted \f$L_2\f$-norm of the data using the control
    * volume to weight the contribution of the data and weight entries to
    * the sum.  That is, the return value is the sum \f$\sqrt{ \sum_i (
    * (data_i * wgt_i) * \bar{(data_i * wgt_i)} cvol_i ) }\f$.  If the control
    * volume is undefined (vol_id < 0), the return value is
    * \f$\sqrt{ \sum_i ( (data_i * wgt_i) * \bar{(data_i * wgt_i)} cvol_i ) }\f$.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   double
   weightedL2Norm(
      const int data_id,
      const int weight_id,
      const int vol_id = -1) const;

   /**
    * Return discrete root mean squared norm of the data.  If the control
    * volume is specified (vol_id >= 0), the return value is the \f$L_2\f$-norm
    * divided by the square root of the sum of the control volumes.  Otherwise,
    * the return value is the \f$L_2\f$-norm divided by the square root of the
    * number of data entries.
    */
   double
   RMSNorm(
      const int data_id,
      const int vol_id = -1) const;

   /**
    * Return discrete weighted root mean squared norm of the data.  If the
    * control volume is specified (vol_id >= 0), the return value is the
    * weighted \f$L_2\f$-norm divided by the square root of the sum of the
    * control volumes.  Otherwise, the return value is the weighted \f$L_2\f$-norm
    * divided by the square root of the number of data entries.
    */
   double
   weightedRMSNorm(
      const int data_id,
      const int weight_id,
      const int vol_id = -1) const;

   /**
    * Return the \f$\max\f$-norm of the data using the control volume to weight
    * the contribution of each data entry to the maximum.  That is, the return
    * value is \f$\max_i ( \sqrt{data_i * \bar{data_i}} )\f$, where the max is
    * over the data elements where \f$cvol_i > 0\f$.  If the control volume is
    * undefined (vol_id < 0), it is ignored during the computation of the
    * maximum.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   double
   maxNorm(
      const int data_id,
      const int vol_id = -1) const;

   /**
    * Return the dot product of the two data arrays using the control volume
    * to weight the contribution of each product to the sum.  That is, the
    * return value is the sum \f$\sum_i ( data1_i * \bar{data2_i} * cvol_i )\f$.
    * If the control volume is undefined (vol_id < 0), it is ignored during
    * the summation.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   dcomplex
   dot(
      const int data1_id,
      const int data2_id,
      const int vol_id = -1) const;

   /**
    * Return the integral of the function represented by the data array.
    * The return value is the sum \f$\sum_i ( data_i * vol_i )\f$.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   dcomplex
   integral(
      const int data_id,
      const int vol_id) const;

private:
   // The following are not implemented
   HierarchyEdgeDataOpsComplex(
      const HierarchyEdgeDataOpsComplex&);
   HierarchyEdgeDataOpsComplex&
   operator = (
      const HierarchyEdgeDataOpsComplex&);

   std::shared_ptr<hier::PatchHierarchy> d_hierarchy;
   int d_coarsest_level;
   int d_finest_level;
   std::vector<std::vector<hier::BoxContainer> >
   d_nonoverlapping_edge_boxes[SAMRAI::MAX_DIM_VAL];

   PatchEdgeDataOpsComplex d_patch_ops;

};

}
}
#endif

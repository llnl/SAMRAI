/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Operations for integer node data on multiple levels.
 *
 ************************************************************************/

#ifndef included_math_HierarchyNodeDataOpsInteger
#define included_math_HierarchyNodeDataOpsInteger

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/math/HierarchyDataOpsInteger.h"
#include "SAMRAI/math/PatchNodeDataOpsInteger.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/PatchHierarchy.h"

#include <iostream>
#include <memory>

namespace SAMRAI {
namespace math {

/**
 * Class HierarchyNodeDataOpsInteger provides a collection of
 * operations that manipulate integer node-centered patch data components over
 * multiple levels in an AMR hierarchy.  It is derived from the abstract
 * base class HierarchyDataOpsInteger which defines the interface to
 * similar operations for node-centered, face-centered, node-centered patch
 * data objects where the data is of type integer.  The operations include
 * basic arithmetic and some ordering operations.  On each patch, the
 * operations are performed by the PatchNodeDataOpsInteger data member.
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
 * and float) and complex node-centered data is provided in the classes
 * HierarchyNodeDataOpsReal and HierarchyNodeDataOpsComplex,
 * respectively.
 *
 * @see PatchNodeDataOpsInteger
 */

class HierarchyNodeDataOpsInteger:public HierarchyDataOpsInteger
{
public:
   /**
    * The constructor for the HierarchyNodeDataOpsInteger class sets
    * the default patch hierarchy and coarsest and finest patch levels
    * in that hierarchy over which operations will be performed.  The
    * hierarchy and operations may be reset using the member fuctions
    * setPatchHierarchy() and resetLevels() below. If no level number
    * arguments are given here, the levels over which the operations will
    * be performed are those already existing in the hierarchy.  If the
    * hierarchy level configuration changes, the operations must be explicitly
    * reset by calling the resetLevels() function.
    *
    * @pre hierarchy
    */
   explicit HierarchyNodeDataOpsInteger(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int coarsest_level = -1,
      const int finest_level = -1);

   /**
    * Virtual destructor for the HierarchyNodeDataOpsInteger class.
    */
   virtual ~HierarchyNodeDataOpsInteger();

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
      const int& alpha,
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
      const int& alpha,
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
      const int& alpha,
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
      const int& alpha,
      const int src1_id,
      const int& beta,
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
      const int& alpha,
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
      const int& alpha,
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
   int
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
   int
   max(
      const int data_id,
      const bool interior_only = true) const;

   /**
    * Set data entries to random values.  See the operations in the
    * array data operation classes for details on the generation of
    * the random values.
    *
    * @pre getPatchHierarchy()
    * @pre (d_coarsest_level >= 0) && (d_finest_level >= d_coarsest_level) &&
    *      (d_finest_level <= getPatchHierarchy()->getFinestLevelNumber())
    */
   void
   setRandomValues(
      const int data_id,
      const int& width,
      const int& low,
      const bool interior_only = true) const;

private:
   // The following are not implemented
   HierarchyNodeDataOpsInteger(
      const HierarchyNodeDataOpsInteger&);
   HierarchyNodeDataOpsInteger&
   operator = (
      const HierarchyNodeDataOpsInteger&);

   std::shared_ptr<hier::PatchHierarchy> d_hierarchy;
   int d_coarsest_level;
   int d_finest_level;
   std::vector<std::vector<hier::BoxContainer> > d_nonoverlapping_node_boxes;

   PatchNodeDataOpsInteger d_patch_ops;

};

}
}
#endif

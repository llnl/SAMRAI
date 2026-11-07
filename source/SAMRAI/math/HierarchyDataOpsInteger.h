/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Interface to operations for integer data on hierarchy.
 *
 ************************************************************************/

#ifndef included_math_HierarchyDataOpsInteger
#define included_math_HierarchyDataOpsInteger

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/PatchHierarchy.h"

#include <iostream>
#include <memory>

namespace SAMRAI {
namespace math {

/**
 * Class HierarchyDataOpsInteger defines the interface to a collection
 * of operations that may be used to manipulate integer numerical patch data
 * components over multiple levels in an AMR hierarchy.  It serves as a base
 * class for subclasses which implement the operations for cell-centered,
 * face-centered, or node-centered data types.  The patch hierarchy and set
 * of levels within that hierarcy over which the operations will be performed
 * are set in the constructor of the subclass.  However, these data members may
 * be changed at any time via the virtual access functions setPatchHierarchy()
 * and resetLevels() below.  The operations include basic arithmetic and
 * some ordering operations.
 *
 * Note that, when it makes sense, an operation accept a boolean argument
 * which indicates whether the operation should be performed on all of the
 * data or just those data elements corresponding to the patch interiors.
 * If no boolean argument is provided, the default behavior is to treat only
 * the patch interiors.  Also, interfaces for similar sets of operations
 * for real (double and float) and complex hierarchy data are defined in the
 * classes HierarchyDataOpsReal and HierarchyDataOpsComplex,
 * respectively.
 */

class HierarchyDataOpsInteger
{
public:
   /**
    * The constructor for the HierarchyDataOpsInteger class.
    */
   HierarchyDataOpsInteger();

   /**
    * Virtual destructor for the HierarchyDataOpsInteger class.
    */
   virtual ~HierarchyDataOpsInteger();

   /**
    * Reset patch hierarchy over which operations occur.
    */
   virtual
   void
   setPatchHierarchy(
      const std::shared_ptr<hier::PatchHierarchy>& hierarchy) = 0;

   /**
    * Reset range of patch levels over which operations occur.
    */
   virtual void
   resetLevels(
      const int coarsest_level,
      const int finest_level) = 0;

   /**
    * Return const pointer to patch hierarchy associated with operations.
    */
   virtual
   const std::shared_ptr<hier::PatchHierarchy>
   getPatchHierarchy() const = 0;

   /**
    * Return the total number of data values for the component on the set
    * of hierarchy levels.  If the boolean argument is true, the number of
    * elements will be summed over patch interiors.  If the boolean argument
    * is false, all elements will be counted (including ghost values)
    * over all patches.
    */
   virtual size_t
   numberOfEntries(
      const int data_id,
      const bool interior_only = true) const = 0;

   /**
    * Copy source data to destination data.
    */
   virtual void
   copyData(
      const int dst_id,
      const int src_id,
      const bool interior_only = true) const = 0;

   /**
    * Swap data pointers (i.e., storage) between two data components.
    */
   virtual void
   swapData(
      const int data1_id,
      const int data2_id) const = 0;

   /**
    * Print data over multiple levels to specified output stream.
    */
   virtual void
   printData(
      const int data_id,
      std::ostream& s,
      const bool interior_only = true) const = 0;

   /**
    * Set data component to given scalar.
    */
   virtual void
   setToScalar(
      const int data_id,
      const int& alpha,
      const bool interior_only = true) const = 0;

   /**
    * Set destination to source multiplied by given scalar, pointwise.
    */
   virtual void
   scale(
      const int dst_id,
      const int& alpha,
      const int src_id,
      const bool interior_only = true) const = 0;

   /**
    * Add scalar to each entry in source data and set destination to result.
    */
   virtual void
   addScalar(
      const int dst_id,
      const int src_id,
      const int& alpha,
      const bool interior_only = true) const = 0;

   /**
    * Set destination to sum of two source components, pointwise.
    */
   virtual void
   add(
      const int dst_id,
      const int src1_id,
      const int src2_id,
      const bool interior_only = true) const = 0;

   /**
    * Subtract second source component from first source component pointwise
    * and set destination data component to result.
    */
   virtual void
   subtract(
      const int dst_id,
      const int src1_id,
      const int src2_id,
      const bool interior_only = true) const = 0;

   /**
    * Set destination component to product of two source components, pointwise.
    */
   virtual void
   multiply(
      const int dst_id,
      const int src1_id,
      const int src2_id,
      const bool interior_only = true) const = 0;

   /**
    * Divide first data component by second source component pointwise
    * and set destination data component to result.
    */
   virtual void
   divide(
      const int dst_id,
      const int src1_id,
      const int src2_id,
      const bool interior_only = true) const = 0;

   /**
    * Set each entry of destination component to reciprocal of corresponding
    * source data component entry.
    */
   virtual void
   reciprocal(
      const int dst_id,
      const int src_id,
      const bool interior_only = true) const = 0;

   /**
    * Set \f$d = \alpha s_1 + \beta s_2\f$, where \f$d\f$ is the destination patch
    * data component and \f$s_1, s_2\f$ are the first and second source components,
    * respectively.  Here \f$\alpha, \beta\f$ are scalar values.
    */
   virtual void
   linearSum(
      const int dst_id,
      const int& alpha,
      const int src1_id,
      const int& beta,
      const int src2_id,
      const bool interior_only = true) const = 0;

   /**
    * Set \f$d = \alpha s_1 + s_2\f$, where \f$d\f$ is the destination patch data
    * component and \f$s_1, s_2\f$ are the first and second source components,
    * respectively.  Here \f$\alpha\f$ is a scalar.
    */
   virtual void
   axpy(
      const int dst_id,
      const int& alpha,
      const int src1_id,
      const int src2_id,
      const bool interior_only = true) const = 0;

   /**
    * Set \f$d = \alpha s_1 - s_2\f$, where \f$d\f$ is the destination patch data
    * component and \f$s_1, s_2\f$ are the first and second source components,
    * respectively.  Here \f$\alpha\f$ is a scalar.
    */
   virtual void
   axmy(
      const int dst_id,
      const int& alpha,
      const int src1_id,
      const int src2_id,
      const bool interior_only = true) const = 0;

   /**
    * Set destination data to absolute value of source data, pointwise.
    */
   virtual void
   abs(
      const int dst_id,
      const int src_id,
      const bool interior_only = true) const = 0;

   /**
    * Return minimum data value over all patches in the collection of levels.
    */
   virtual int
   min(
      const int data_id,
      const bool interior_only = true) const = 0;

   /**
    * Return maximum data value over all patches in the collection of levels.
    */
   virtual int
   max(
      const int data_id,
      const bool interior_only = true) const = 0;

   /**
    * Set data entries to random values.  See the operations in the
    * array data operation classes for each data type for details on the
    * generation of the random values.
    */
   virtual void
   setRandomValues(
      const int data_id,
      const int& width,
      const int& low,
      const bool interior_only = true) const = 0;

private:
   // The following are not implemented
   HierarchyDataOpsInteger(
      const HierarchyDataOpsInteger&);
   HierarchyDataOpsInteger&
   operator = (
      const HierarchyDataOpsInteger&);

};

}
}
#endif

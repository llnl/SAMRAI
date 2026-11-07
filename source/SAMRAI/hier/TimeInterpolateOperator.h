/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Abstract base class for time interpolation operators.
 *
 ************************************************************************/

#ifndef included_hier_TimeInterpolateOperator
#define included_hier_TimeInterpolateOperator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/Variable.h"

#include <string>
#include <memory>

namespace SAMRAI {
namespace hier {

/**
 * Class TimeInterpolateOperator is an abstract base class for each
 * time interpolation operator used in the SAMRAI framework.  This class
 * defines the interface between numerical interpolation routines and the
 * rest of the framework.  Each concrete time interpolation operator subclass
 * provide two operations:
 *
 *
 *
 * - @b (1) an implementation of the time interpolation operation
 *            appropriate for its corresponding patch data type.
 * - @b (2) a function that determines whether or not the operator
 *             matches an arbitrary request for a time interpolation operator.
 *
 *
 *
 * To add a new time interpolation operator (either for a new patch data type
 * or for a new time interpolation routine on an existing type), define
 * the operator by inheriting from this abstract base class.  The operator
 * subclass must implement the interpolation operation in the timeInterpolate()
 * function.  Then, the new operator must be added to the operator list
 * for the appropriate transfer geometry object using the
 * BaseGridGeometry::addTimeInterpolateOperator() function.
 *
 * Although time interpolation operators usually depend only on patch data
 * centering and data type and not the mesh coordinate system, they are
 * defined in the @em geometry package.
 *
 * @see TransferOperatorRegistry
 */

class TimeInterpolateOperator
{
public:
   /**
    * The default constructor for the coarsening operator does
    * nothing interesting.
    */
   TimeInterpolateOperator(
      const std::string& name = "STD_LINEAR_TIME_INTERPOLATE");

   /**
    * The virtual destructor for the coarsening operator does
    * nothing interesting.
    */
   virtual ~TimeInterpolateOperator() = default;

   /**
    * Return name std::string identifier of the time interpolate operation.
    */
   const std::string&
   getOperatorName() const
   {
      return d_name;
   }

   /**
    * Perform time interpolation between two patch data sources
    * and place result in the destination patch data.  Time interpolation
    * is performed on the intersection of the destination patch data and
    * the input box.  The time to which data is interpolated is provided
    * by the destination data.
    */
   virtual void
   timeInterpolate(
      PatchData& dst_data,
      const Box& where,
      const BoxOverlap& overlap, 
      const PatchData& src_data_old,
      const PatchData& src_data_new) const = 0;

private:
   // Neither of these is implemented.
   TimeInterpolateOperator(
      const TimeInterpolateOperator&) = delete;
   TimeInterpolateOperator&
   operator = (
      const TimeInterpolateOperator&) = delete;

   const std::string d_name;
};

}
}
#endif

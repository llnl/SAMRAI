/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Linear time interp operator for side-centered complex data.
 *
 ************************************************************************/

#ifndef included_pdat_SideComplexLinearTimeInterpolateOp
#define included_pdat_SideComplexLinearTimeInterpolateOp

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/TimeInterpolateOperator.h"

#include <string>

namespace SAMRAI {
namespace pdat {

/**
 * Class SideComplexLinearTimeInterpolateOp implements standard
 * linear time interpolation for side-centered complex patch data.
 * It is derived from the hier::TimeInterpolateOperator base class.
 * The interpolation uses FORTRAN numerical routines.
 *
 * @see hier::TimeInterpolateOperator
 */

class SideComplexLinearTimeInterpolateOp:
   public hier::TimeInterpolateOperator
{
public:
   /**
    * Uninteresting default constructor.
    */
   SideComplexLinearTimeInterpolateOp();

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~SideComplexLinearTimeInterpolateOp();

   /**
    * Perform linear time interpolation between two side-centered complex
    * patch data sources and place result in the destination patch data.
    * Time interpolation is performed on the intersection of the destination
    * patch data and the input box.  The time to which data is interpolated
    * is provided by the destination data.
    *
    * @pre dynamic_cast<const SideData<float> *>(&src_data_old) != 0
    * @pre dynamic_cast<const SideData<float> *>(&src_data_new) != 0
    * @pre dynamic_cast<SideData<float> *>(&dst_data) != 0
    * @pre (where * src_data_old->getGhostBox()).isSpatiallyEqual(where)
    * @pre (where * src_data_new->getGhostBox()).isSpatiallyEqual(where)
    * @pre (where * dst_data->getGhostBox()).isSpatiallyEqual(where)
    * @pre (dst_data.getDim() == where.getDim()) &&
    *      (dst_data.getDim() == src_data_old.getDim()) &&
    *      (dst_data.getDim() == src_data_new.getDim())
    * @pre direction vector of dst_data is smaller than that of src_data_old or
    *      src_data_new
    */
   void
   timeInterpolate(
      hier::PatchData& dst_data,
      const hier::Box& where,
      const hier::BoxOverlap& overlap,
      const hier::PatchData& src_data_old,
      const hier::PatchData& src_data_new) const;

private:
};

}
}
#endif

/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Conservative linear refine operator for face-centered
 *                float data on a Cartesian mesh.
 *
 ************************************************************************/

#ifndef included_geom_CartesianFaceFloatConservativeLinearRefine
#define included_geom_CartesianFaceFloatConservativeLinearRefine

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"

#include <string>
#include <memory>

namespace SAMRAI {
namespace geom {

/**
 * Class CartesianFaceFloatConservativeLinearRefine implements
 * conservative linear interpolation for face-centered float patch data
 * defined over a Cartesian mesh.  It is derived from the base class
 * hier::RefineOperator.  The numerical operations for the interpolation
 * use FORTRAN numerical routines.
 *
 * @see hier::RefineOperator
 */

class CartesianFaceFloatConservativeLinearRefine:
   public hier::RefineOperator
{
public:
   /**
    * Uninteresting default constructor.
    */
   CartesianFaceFloatConservativeLinearRefine();

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~CartesianFaceFloatConservativeLinearRefine();

   /**
    * The priority of face-centered float conservative linear is 0.
    * It will be performed before any user-defined interpolation operations.
    */
   int
   getOperatorPriority() const;

   /**
    * The stencil width of the conservative linear interpolation operator is
    * the vector of ones.
    */
   hier::IntVector
   getStencilWidth(
      const tbox::Dimension& dim) const;

   /**
    * Refine the source component on the coarse patch to the destination
    * component on the fine patch using the face-centered float conservative
    * linear interpolation operator.  Interpolation is performed on the
    * intersection of the destination patch and the boxes contained in
    * fine_overlap.  It is assumed that the coarse patch contains sufficient
    * data for the stencil width of the refinement operator.
    *
    * @pre (fine.getDim() == coarse.getDim()) &&
    *      (fine.getDim() == ratio.getDim())
    * @pre dynamic_cast<const pdat::FaceOverlap *>(&fine_overlap) != 0
    * @pre coarse.getPatchData(src_component) is actually a std::shared_ptr<pdat::FaceData<float> >
    * @pre fine.getPatchData(dst_component) is actually a std::shared_ptr<pdat::FaceData<float> >
    * @pre coarse.getPatchData(src_component)->getDepth() == fine.getPatchData(dst_component)->getDepth()
    * @pre (fine.getDim().getValue() == 1) ||
    *      (fine.getDim().getValue() == 2) || (fine.getDim().getValue() == 3)
    */
   void
   refine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const int dst_component,
      const int src_component,
      const hier::BoxOverlap& fine_overlap,
      const hier::IntVector& ratio) const;
};

}
}
#endif

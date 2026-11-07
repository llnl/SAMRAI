/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Constant refine operator for side-centered integer data on
 *                a  mesh.
 *
 ************************************************************************/

#ifndef included_pdat_SideIntegerConstantRefine
#define included_pdat_SideIntegerConstantRefine

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/RefineOperator.h"

#include <string>

namespace SAMRAI {
namespace pdat {

/**
 * Class SideIntegerConstantRefine implements constant
 * interpolation for side-centered integer patch data defined over a
 * mesh.  It is derived from the hier::RefineOperator base class.
 * The numerical operations for interpolation use FORTRAN numerical routines.
 *
 * @see hier::RefineOperator
 */

class SideIntegerConstantRefine:
   public hier::RefineOperator
{
public:
   /**
    * Uninteresting default constructor.
    */
   SideIntegerConstantRefine();

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~SideIntegerConstantRefine();

   /**
    * The priority of side-centered integer constant interpolation is 0.
    * It will be performed before any user-defined interpolation operations.
    */
   int
   getOperatorPriority() const;

   /**
    * The stencil width of the constant interpolation operator is the vector
    * of zeros.  That is, its stencil does not extend outside the fine box.
    */
   hier::IntVector
   getStencilWidth(
      const tbox::Dimension& dim) const;

   /**
    * Refine the source component on the coarse patch to the destination
    * component on the fine patch using the side-centered integer constant
    * interpolation operator.  Interpolation is performed on the intersection
    * of the destination patch and the boxes contained in fine_overlap.
    * It is assumed that the coarse patch contains sufficient data for the
    * stencil width of the refinement operator.
    *
    * @pre dynamic_cast<const SideOverlap *>(&fine_overlap) != 0
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

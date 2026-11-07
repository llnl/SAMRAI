/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Constant refine operator for cell-centered integer data on
 *                a mesh.
 *
 ************************************************************************/

#ifndef included_pdat_CellIntegerConstantRefine
#define included_pdat_CellIntegerConstantRefine

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/RefineOperator.h"

#include <string>
#include <memory>

namespace SAMRAI {
namespace pdat {

/**
 * Class CellIntegerConstantRefine implements constant
 * interpolation for cell-centered integer patch data defined over a
 * mesh.  It is derived from the hier::RefineOperator base class.
 * The numerical operations for interpolation use FORTRAN numerical routines.
 *
 * @see hier::RefineOperator
 */

class CellIntegerConstantRefine:
   public hier::RefineOperator
{
public:
   /**
    * Uninteresting default constructor.
    */
   CellIntegerConstantRefine();

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~CellIntegerConstantRefine();

   /**
    * The priority of cell-centered integer constant interpolation is 0.
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
    * component on the fine patch using the cell-centered integer constant
    * interpolation operator.  Interpolation is performed on the intersection
    * of the destination patch and the boxes contained in fine_overlap.
    * It is assumed that the coarse patch contains sufficient data for the
    * stencil width of the refinement operator.
    *
    * @pre dynamic_cast<const CellOverlap *>(&fine_overlap) != 0
    */
   void
   refine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const int dst_component,
      const int src_component,
      const hier::BoxOverlap& fine_overlap,
      const hier::IntVector& ratio) const;

   /**
    * Refine the source component on the coarse patch to the destination
    * component on the fine patch using the cell-centered integer constant
    * interpolation operator.  Interpolation is performed on the intersection
    * of the destination patch and the fine box.   It is assumed that the
    * coarse patch contains sufficient data for the stencil width of the
    * refinement operator.  This differs from the above refine() method
    * only in that it operates on a single fine box instead of a BoxOverlap.
    */
   void
   refine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const int dst_component,
      const int src_component,
      const hier::Box& fine_box,
      const hier::IntVector& ratio) const;
};

}
}
#endif

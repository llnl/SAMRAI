/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Constant averaging operator for node-centered float data on
 *                a  mesh.
 *
 ************************************************************************/

#ifndef included_pdat_NodeFloatInjection
#define included_pdat_NodeFloatInjection

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"

#include <string>

namespace SAMRAI {
namespace pdat {

/**
 * Class NodeFloatInjection implements constant
 * averaging (i.e., injection) for node-centered float patch data defined
 * over a  mesh.  It is derived from the hier::CoarsenOperator base
 * class.  The numerical operations for theaveraging use FORTRAN numerical
 * routines.
 *
 * @see hier::CoarsenOperator
 */

class NodeFloatInjection:
   public hier::CoarsenOperator
{
public:
   /**
    * Uninteresting default constructor.
    */
   NodeFloatInjection();

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~NodeFloatInjection();

   /**
    * The priority of node-centered constant averaging is 0.
    * It will be performed before any user-defined coarsen operations.
    */
   int
   getOperatorPriority() const;

   /**
    * The stencil width of the constant averaging operator is the vector of
    * zeros.  That is, its stencil does not extend outside the fine box.
    */
   hier::IntVector
   getStencilWidth(
      const tbox::Dimension& dim) const;

   /**
    * Coarsen the source component on the fine patch to the destination
    * component on the coarse patch using the node-centered float constant
    * averaging operator.  Coarsening is performed on the intersection of
    * the destination patch and the coarse box.  It is assumed that the
    * fine patch contains sufficient data for the stencil width of the
    * coarsening operator.
    */
   void
   coarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const int dst_component,
      const int src_component,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio) const;
};

}
}
#endif

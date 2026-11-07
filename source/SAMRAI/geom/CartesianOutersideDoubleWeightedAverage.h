/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Weighted averaging operator for outerside double data on
 *                a Cartesian mesh.
 *
 ************************************************************************/

#ifndef included_geom_CartesianOutersideDoubleWeightedAverage
#define included_geom_CartesianOutersideDoubleWeightedAverage

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"

#include <string>
#include <memory>

namespace SAMRAI {
namespace geom {

/**
 * Class CartesianOutersideDoubleWeightedAverage implements conservative
 * side-weighted averaging for outerside double patch data defined over
 * a Cartesian mesh.  It is derived from the hier::CoarsenOperator base class.
 * The numerical operations for theaveraging use FORTRAN numerical routines.
 *
 * @see hier::CoarsenOperator
 */

class CartesianOutersideDoubleWeightedAverage:
   public hier::CoarsenOperator
{
public:
   /**
    * Uninteresting default constructor.
    */
   CartesianOutersideDoubleWeightedAverage();

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~CartesianOutersideDoubleWeightedAverage();

   /**
    * The priority of outerside double weighted averaging is 0.
    * It will be performed before any user-defined coarsen operations.
    */
   int
   getOperatorPriority() const;

   /**
    * The stencil width of the weighted averaging operator is the vector of
    * zeros.  That is, its stencil does not extend outside the fine box.
    */
   hier::IntVector
   getStencilWidth(
      const tbox::Dimension& dim) const;

   /**
    * Coarsen the source component on the fine patch to the destination
    * component on the coarse patch using the outerside double weighted
    * averaging operator.  Coarsening is performed on the intersection of
    * the destination patch and the coarse box.  It is assumed that the
    * fine patch contains sufficient data for the stencil width of the
    * coarsening operator.
    *
    * @pre (fine.getDim() == coarse.getDim()) &&
    *      (fine.getDim() == coarse_box.getDim()) &&
    *      (fine.getDim() == ratio.getDim())
    * @pre fine.getPatchData(src_component) is actually a std::shared_ptr<pdat::OutersideData<double> >
    * @pre coarse.getPatchData(dst_component) is actually a std::shared_ptr<pdat::OutersideData<double> >
    * @pre fine.getPatchData(src_component)->getDepth() == coarse.getPatchData(dst_component)->getDepth()
    * @pre (fine.getDim().getValue() == 1) ||
    *      (fine.getDim().getValue() == 2) || (fine.getDim().getValue() == 3)
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

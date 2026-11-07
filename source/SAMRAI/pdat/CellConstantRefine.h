/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Constant refine operator for cell-centered double data on
 *                a  mesh.
 *
 ************************************************************************/

#ifndef included_pdat_CellConstantRefine
#define included_pdat_CellConstantRefine

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/InvokeOne.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/Utilities.h"

#include <float.h>
#include <math.h>

#include <string>
#include <memory>

extern "C" {

#ifdef __INTEL_COMPILER
#pragma warning (disable:1419)
#endif

// in conrefine1d.f:
void SAMRAI_F77_FUNC(conrefcelldoub1d, CONREFCELLDOUB1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const double *, double *);

void SAMRAI_F77_FUNC(conrefcellflot1d, CONREFCELLFLOT1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const float *, float *);

void SAMRAI_F77_FUNC(conrefcellcplx1d, CONREFCELLCPLX1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);

void SAMRAI_F77_FUNC(conrefcellintg1d, CONREFCELLINTG1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const int *, int *);

// in conrefine2d.f:
void SAMRAI_F77_FUNC(conrefcelldoub2d, CONREFCELLDOUB2D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const double *, double *);

void SAMRAI_F77_FUNC(conrefcellflot2d, CONREFCELLFLOT2D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const float *, float *);


void SAMRAI_F77_FUNC(conrefcellcplx2d, CONREFCELLCPLX2D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);

void SAMRAI_F77_FUNC(conrefcellintg2d, CONREFCELLINTG2D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const int *, int *);

// in conrefine3d.f:
void SAMRAI_F77_FUNC(conrefcelldoub3d, CONREFCELLDOUB3D) (const int&, const int&,
   const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *,
   const double *, double *);

void SAMRAI_F77_FUNC(conrefcellflot3d, CONREFCELLFLOT3D) (const int&, const int&,
   const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *,
   const float *, float *);

void SAMRAI_F77_FUNC(conrefcellcplx3d, CONREFCELLCPLX3D) (const int&, const int&,
   const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);

void SAMRAI_F77_FUNC(conrefcellintg3d, CONREFCELLINTG3D) (const int&, const int&,
   const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *,
   const int *, int *);
} // extern "C"

namespace SAMRAI {
namespace pdat {


template<typename T>
void Call1dFortranCell(const int& ifirstc, const int& ilastc,
   const int& ifirstf, const int& ilastf,
   const int& cilo, const int& cihi,
   const int& filo, const int& fihi,
   const int *ratio,
   const T *arrayc, T *arrayf)
{
   invokeOneOfFour(SAMRAI_F77_FUNC(conrefcelldoub1d, CONREFCELLDOUB1D),
              SAMRAI_F77_FUNC(conrefcellflot1d, CONREFCELLFLOT1D), 
              SAMRAI_F77_FUNC(conrefcellcplx1d, CONREFCELLCPLX1D),
              SAMRAI_F77_FUNC(conrefcellintg1d, CONREFCELLINTG1D),
         ifirstc, ilastc, ifirstf, ilastf,
         cilo, cihi,filo, fihi,
         ratio,
         arrayc,arrayf);

}

template<typename T>
void Call2dFortranCell(const int& ifirstc0, const int& ifirstc1,
   const int& ilastc0, const int& ilastc1,
   const int& ifirstf0, const int& ifirstf1, const int& ilastf0, const int& ilastf1,
   const int& cilo0, const int& cilo1, const int& cihi0, const int& cihi1,
   const int& filo0, const int& filo1, const int& fihi0, const int& fihi1,
   const int *ratio,
   const T *arrayc, T *arrayf)
{
    invokeOneOfFour(SAMRAI_F77_FUNC(conrefcelldoub2d, CONREFCELLDOUB2D),
              SAMRAI_F77_FUNC(conrefcellflot2d, CONREFCELLFLOT2D), 
              SAMRAI_F77_FUNC(conrefcellcplx2d, CONREFCELLCPLX2D),
              SAMRAI_F77_FUNC(conrefcellintg2d, CONREFCELLINTG2D),
         ifirstc0,ifirstc1,
         ilastc0,ilastc1,
         ifirstf0,ifirstf1,
         ilastf0,ilastf1,
         cilo0,cilo1,cihi0,cihi1,
         filo0,filo1,fihi0,fihi1,
         ratio,
         arrayc,arrayf);
}

template<typename T>
void Call3dFortranCell(const int& ifirstc0, const int& ifirstc1, const int& ifirstc2,
   const int& ilastc0, const int& ilastc1, const int& ilastc2,
   const int& ifirstf0, const int& ifirstf1, const int& ifirstf2,
   const int& ilastf0, const int& ilastf1, const int& ilastf2,
   const int& cilo0, const int& cilo1, const int& cilo2,
   const int& cihi0, const int& cihi1, const int& cihi2,
   const int& filo0, const int& filo1, const int& filo2,
   const int& fihi0, const int& fihi1, const int& fihi2,
   const int *ratio,
   const T *arrayc, T *arrayf)
{

    invokeOneOfFour(SAMRAI_F77_FUNC(conrefcelldoub3d, CONREFCELLDOUB3D),
              SAMRAI_F77_FUNC(conrefcellflot3d, CONREFCELLFLOT3D), 
              SAMRAI_F77_FUNC(conrefcellcplx3d, CONREFCELLCPLX3D),
              SAMRAI_F77_FUNC(conrefcellintg3d, CONREFCELLINTG3D),
         ifirstc0,ifirstc1,ifirstc2,
         ilastc0,ilastc1,ilastc2,
         ifirstf0,ifirstf1,ifirstf2,
         ilastf0,ilastf1,ilastf2,
         cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
         filo0,filo1,filo2,fihi0,fihi1,fihi2,
         ratio,
         arrayc,arrayf);
}


/**
 * Class CellConstantRefine implements constant
 * interpolation for cell-centered patch data defined over a
 * mesh.  It is derived from the hier::RefineOperator base class.
 * The numerical operations for interpolation use FORTRAN numerical routines.
 *
 * @see hier::RefineOperator
 */
template<typename T>  // one of double,float,dcomplex,int  
class CellConstantRefine:
   public hier::RefineOperator
{
public:
   /**
    * Uninteresting default constructor.
    */
   CellConstantRefine():
   hier::RefineOperator("CONSTANT_REFINE")
   {
   }
   
   /**
    * Uninteresting destructor.
    */
   virtual ~CellConstantRefine()
   {
   }

   /**
    * The priority of cell-centered constant interpolation is 0.
    * It will be performed before any user-defined interpolation operations.
    */
   int
   getOperatorPriority() const
   {
      return 0;
   }


   /**
    * The stencil width of the constant interpolation operator is the vector
    * of zeros.  That is, its stencil does not extend outside the fine box.
    */
   hier::IntVector
   getStencilWidth(const tbox::Dimension& dim) const
   {
      return hier::IntVector::getZero(dim);
   }

   /**
    * Refine the source component on the coarse patch to the destination
    * component on the fine patch using the cell-centered double constant
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
      const hier::IntVector& ratio) const
   {
      const CellOverlap* t_overlap = CPP_CAST<const CellOverlap *>(&fine_overlap);

      TBOX_ASSERT(t_overlap != 0);

      const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer();
      for (hier::BoxContainer::const_iterator b = boxes.begin();
           b != boxes.end(); ++b) {
         refine(fine,
            coarse,
            dst_component,
            src_component,
            *b,
            ratio);
      }
   }

   /**
    * Refine the source component on the coarse patch to the destination
    * component on the fine patch using the cell-centered double constant
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


} // namespace pdat
} // namespace SAMRAI

#include "SAMRAI/pdat/CellConstantRefine.cpp"
#endif

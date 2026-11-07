/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Constant refine operator for side-centered double data on
 *                a  mesh.
 *
 ************************************************************************/

#ifndef included_pdat_SideConstantRefine
#define included_pdat_SideConstantRefine

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"

#include <string>

extern "C" {

#ifdef __INTEL_COMPILER
#pragma warning (disable:1419)
#endif

// in conrefine1d.f:
void SAMRAI_F77_FUNC(conrefsidedoub1d, CONREFSIDEDOUB1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const double *, double *);
void SAMRAI_F77_FUNC(conrefsideflot1d, CONREFSIDEFLOT1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const float *, float *);
void SAMRAI_F77_FUNC(conrefsidecplx1d, CONREFSIDECPLX1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);
void SAMRAI_F77_FUNC(conrefsideintg1d, CONREFSIDEINTG1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const int *, int *);

// in conrefine2d.f:
// 2d0
void SAMRAI_F77_FUNC(conrefsidedoub2d0, CONREFSIDEDOUB2D0) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const double *, double *);
void SAMRAI_F77_FUNC(conrefsideflot2d0, CONREFSIDEFLOT2D0) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const float *, float *);
void SAMRAI_F77_FUNC(conrefsidecplx2d0, CONREFSIDECPLX2D0) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);
void SAMRAI_F77_FUNC(conrefsideintg2d0, CONREFSIDEINTG2D0) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const int *, int *);

//2d1
void SAMRAI_F77_FUNC(conrefsidedoub2d1, CONREFSIDEDOUB2D1) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const double *, double *);
void SAMRAI_F77_FUNC(conrefsideflot2d1, CONREFSIDEFLOT2D1) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const float *, float *);
void SAMRAI_F77_FUNC(conrefsidecplx2d1, CONREFSIDECPLX2D1) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);
void SAMRAI_F77_FUNC(conrefsideintg2d1, CONREFSIDEINTG2D1) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const int *, int *);
// in conrefine3d.f:
// 3d0
void SAMRAI_F77_FUNC(conrefsidedoub3d0, CONREFSIDEDOUB3D0) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conrefsideflot3d0, CONREFSIDEflot3D0) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conrefsidecplx3d0, CONREFSIDECPLX3D0) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conrefsideintg3d0, CONREFSIDEINTG3D0) (const int&, const int&,
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

//3d1
void SAMRAI_F77_FUNC(conrefsidedoub3d1, CONREFSIDEDOUB3D1) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conrefsideflot3d1, CONREFSIDEFLOT3D1) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conrefsidecplx3d1, CONREFSIDECPLX3D1) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conrefsideintg3d1, CONREFSIDEINTG3D1) (const int&, const int&,
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

//3d2
void SAMRAI_F77_FUNC(conrefsidedoub3d2, CONREFSIDEDOUB3D2) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conrefsideflot3d2, CONREFSIDEFLOT3D2) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conrefsidecplx3d2, CONREFSIDECPLX3D2) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conrefsideintg3d2, CONREFSIDEINTG3D2) (const int&, const int&,
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
void Call1dFortranSide(const int& ifirstc, const int& ilastc,
   const int& ifirstf, const int& ilastf,
   const int& cilo, const int& cihi,
   const int& filo, const int& fihi,
   const int *ratio,
   const T *arrayc, T *arrayf)
{
   invokeOneOfFour(SAMRAI_F77_FUNC(conrefsidedoub1d, CONREFSIDEDOUB1D),
              SAMRAI_F77_FUNC(conrefsideflot1d, CONREFSIDEFLOT1D), 
              SAMRAI_F77_FUNC(conrefsidecplx1d, CONREFSIDECPLX1D),
              SAMRAI_F77_FUNC(conrefsideintg1d, CONREFSIDEINTG1D),
         ifirstc, ilastc, ifirstf, ilastf,
         cilo, cihi,filo, fihi,
         ratio,
         arrayc,arrayf);

}

template<typename T>
void Call2dFortranSide_d0(const int& ifirstc0, const int& ifirstc1,
   const int& ilastc0, const int& ilastc1,
   const int& ifirstf0, const int& ifirstf1, const int& ilastf0, const int& ilastf1,
   const int& cilo0, const int& cilo1, const int& cihi0, const int& cihi1,
   const int& filo0, const int& filo1, const int& fihi0, const int& fihi1,
   const int *ratio,
   const T *arrayc, T *arrayf)
{
    invokeOneOfFour(SAMRAI_F77_FUNC(conrefsidedoub2d0, CONREFSIDEDOUB2D0),
              SAMRAI_F77_FUNC(conrefsideflot2d0, CONREFSIDEFLOT2D0), 
              SAMRAI_F77_FUNC(conrefsidecplx2d0, CONREFSIDECPLX2D0),
              SAMRAI_F77_FUNC(conrefsideintg2d0, CONREFSIDEINTG2D0),
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
void Call2dFortranSide_d1(const int& ifirstc0, const int& ifirstc1,
   const int& ilastc0, const int& ilastc1,
   const int& ifirstf0, const int& ifirstf1, const int& ilastf0, const int& ilastf1,
   const int& cilo0, const int& cilo1, const int& cihi0, const int& cihi1,
   const int& filo0, const int& filo1, const int& fihi0, const int& fihi1,
   const int *ratio,
   const T *arrayc, T *arrayf)
{
    invokeOneOfFour(SAMRAI_F77_FUNC(conrefsidedoub2d1, CONREFSIDEDOUB2D1),
              SAMRAI_F77_FUNC(conrefsideflot2d1, CONREFSIDEFLOT2D1), 
              SAMRAI_F77_FUNC(conrefsidecplx2d1, CONREFSIDECPLX2D1),
              SAMRAI_F77_FUNC(conrefsideintg2d1, CONREFSIDEINTG2D1),
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
void Call3dFortranSide_d0(const int& ifirstc0, const int& ifirstc1, const int& ifirstc2,
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

    invokeOneOfFour(SAMRAI_F77_FUNC(conrefsidedoub3d0, CONREFSIDEDOUB3D0),
              SAMRAI_F77_FUNC(conrefsideflot3d0, CONREFSIDEFLOT3D0), 
              SAMRAI_F77_FUNC(conrefsidecplx3d0, CONREFSIDECPLX3D0),
              SAMRAI_F77_FUNC(conrefsideintg3d0, CONREFSIDEINTG3D0),
         ifirstc0,ifirstc1,ifirstc2,
         ilastc0,ilastc1,ilastc2,
         ifirstf0,ifirstf1,ifirstf2,
         ilastf0,ilastf1,ilastf2,
         cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
         filo0,filo1,filo2,fihi0,fihi1,fihi2,
         ratio,
         arrayc,arrayf);
}

template<typename T>
void Call3dFortranSide_d1(const int& ifirstc0, const int& ifirstc1, const int& ifirstc2,
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

    invokeOneOfFour(SAMRAI_F77_FUNC(conrefsidedoub3d1, CONREFSIDEDOUB3D1),
              SAMRAI_F77_FUNC(conrefsideflot3d1, CONREFSIDEFLOT3D1), 
              SAMRAI_F77_FUNC(conrefsidecplx3d1, CONREFSIDECPLX3D1),
              SAMRAI_F77_FUNC(conrefsideintg3d1, CONREFSIDEINTG3D1),
         ifirstc0,ifirstc1,ifirstc2,
         ilastc0,ilastc1,ilastc2,
         ifirstf0,ifirstf1,ifirstf2,
         ilastf0,ilastf1,ilastf2,
         cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
         filo0,filo1,filo2,fihi0,fihi1,fihi2,
         ratio,
         arrayc,arrayf);
}

template<typename T>
void Call3dFortranSide_d2(const int& ifirstc0, const int& ifirstc1, const int& ifirstc2,
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

    invokeOneOfFour(SAMRAI_F77_FUNC(conrefsidedoub3d2, CONREFSIDEDOUB3D2),
              SAMRAI_F77_FUNC(conrefsideflot3d2, CONREFSIDEFLOT3D2), 
              SAMRAI_F77_FUNC(conrefsidecplx3d2, CONREFSIDECPLX3D2),
              SAMRAI_F77_FUNC(conrefsideintg3d2, CONREFSIDEINTG3D2),
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
 * Class SideConstantRefine implements constant
 * interpolation for side-centered double patch data defined over a
 * mesh.  It is derived from the hier::RefineOperator base class.
 * The numerical operations for interpolation use FORTRAN numerical routines.
 *
 * @see hier::RefineOperator
 */

template<typename T> // one of double,float,dcomplex,int
class SideConstantRefine:
   public hier::RefineOperator
{
public:
   /**
    * Uninteresting default constructor.
    */
   SideConstantRefine():
      hier::RefineOperator("CONSTANT_REFINE")
   {
   }

   /**
    * Uninteresting destructor.
    */
   virtual ~SideConstantRefine()
   {
   }

   /**
    * The priority of side-centered double constant interpolation is 0.
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
    * component on the fine patch using the side-centered double constant
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
#include "SAMRAI/pdat/SideConstantRefine.cpp"
#endif

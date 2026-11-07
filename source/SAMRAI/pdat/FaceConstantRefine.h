/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Constant refine operator for face-centered double data on
 *                a  mesh.
 *
 ************************************************************************/

#ifndef included_pdat_FaceConstantRefine
#define included_pdat_FaceConstantRefine

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/pdat/InvokeOne.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"

#include <string>
#include <memory>

extern "C" {

#ifdef __INTEL_COMPILER
#pragma warning (disable:1419)
#endif

// in conrefine1d.f:
void SAMRAI_F77_FUNC(conreffacedoub1d, CONREFFACEDOUB1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const double *, double *);
void SAMRAI_F77_FUNC(conreffaceflot1d, CONREFFACEFLOT1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const float *, float *);
void SAMRAI_F77_FUNC(conreffacecplx1d, CONREFFACECPLX1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);
void SAMRAI_F77_FUNC(conreffaceintg1d, CONREFFACEINTG1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const int *, int *);

// in conrefine2d.f:
// 2d0
void SAMRAI_F77_FUNC(conreffacedoub2d0, CONREFFACEDOUB2D0) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const double *, double *);
void SAMRAI_F77_FUNC(conreffaceflot2d0, CONREFFACEFLOT2D0) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const float *, float *);
void SAMRAI_F77_FUNC(conreffacecplx2d0, CONREFFACECPLX2D0) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);
void SAMRAI_F77_FUNC(conreffaceintg2d0, CONREFFACEINTG2D0) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const int *, int *);

//2d1
void SAMRAI_F77_FUNC(conreffacedoub2d1, CONREFFACEDOUB2D1) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const double *, double *);
void SAMRAI_F77_FUNC(conreffaceflot2d1, CONREFFACEFLOT2D1) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const float *, float *);
void SAMRAI_F77_FUNC(conreffacecplx2d1, CONREFFACECPLX2D1) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);
void SAMRAI_F77_FUNC(conreffaceintg2d1, CONREFFACEINTG2D1) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const int *, int *);
// in conrefine3d.f:
// 3d0
void SAMRAI_F77_FUNC(conreffacedoub3d0, CONREFFACEDOUB3D0) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conreffaceflot3d0, CONREFFACEflot3D0) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conreffacecplx3d0, CONREFFACECPLX3D0) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conreffaceintg3d0, CONREFFACEINTG3D0) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conreffacedoub3d1, CONREFFACEDOUB3D1) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conreffaceflot3d1, CONREFFACEFLOT3D1) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conreffacecplx3d1, CONREFFACECPLX3D1) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conreffaceintg3d1, CONREFFACEINTG3D1) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conreffacedoub3d2, CONREFFACEDOUB3D2) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conreffaceflot3d2, CONREFFACEFLOT3D2) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conreffacecplx3d2, CONREFFACECPLX3D2) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conreffaceintg3d2, CONREFFACEINTG3D2) (const int&, const int&,
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
void Call1dFortranFace(const int& ifirstc, const int& ilastc,
   const int& ifirstf, const int& ilastf,
   const int& cilo, const int& cihi,
   const int& filo, const int& fihi,
   const int *ratio,
   const T *arrayc, T *arrayf)
{
   invokeOneOfFour(SAMRAI_F77_FUNC(conreffacedoub1d, CONREFFACEDOUB1D),
              SAMRAI_F77_FUNC(conreffaceflot1d, CONREFFACEFLOT1D), 
              SAMRAI_F77_FUNC(conreffacecplx1d, CONREFFACECPLX1D),
              SAMRAI_F77_FUNC(conreffaceintg1d, CONREFFACEINTG1D),
         ifirstc, ilastc, ifirstf, ilastf,
         cilo, cihi,filo, fihi,
         ratio,
         arrayc,arrayf);

}

template<typename T>
void Call2dFortranFace_d0(const int& ifirstc0, const int& ifirstc1,
   const int& ilastc0, const int& ilastc1,
   const int& ifirstf0, const int& ifirstf1, const int& ilastf0, const int& ilastf1,
   const int& cilo0, const int& cilo1, const int& cihi0, const int& cihi1,
   const int& filo0, const int& filo1, const int& fihi0, const int& fihi1,
   const int *ratio,
   const T *arrayc, T *arrayf)
{
    invokeOneOfFour(SAMRAI_F77_FUNC(conreffacedoub2d0, CONREFFACEDOUB2D0),
              SAMRAI_F77_FUNC(conreffaceflot2d0, CONREFFACEFLOT2D0), 
              SAMRAI_F77_FUNC(conreffacecplx2d0, CONREFFACECPLX2D0),
              SAMRAI_F77_FUNC(conreffaceintg2d0, CONREFFACEINTG2D0),
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
void Call2dFortranFace_d1(const int& ifirstc0, const int& ifirstc1,
   const int& ilastc0, const int& ilastc1,
   const int& ifirstf0, const int& ifirstf1, const int& ilastf0, const int& ilastf1,
   const int& cilo0, const int& cilo1, const int& cihi0, const int& cihi1,
   const int& filo0, const int& filo1, const int& fihi0, const int& fihi1,
   const int *ratio,
   const T *arrayc, T *arrayf)
{
    invokeOneOfFour(SAMRAI_F77_FUNC(conreffacedoub2d1, CONREFFACEDOUB2D1),
              SAMRAI_F77_FUNC(conreffaceflot2d1, CONREFFACEFLOT2D1), 
              SAMRAI_F77_FUNC(conreffacecplx2d1, CONREFFACECPLX2D1),
              SAMRAI_F77_FUNC(conreffaceintg2d1, CONREFFACEINTG2D1),
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
void Call3dFortranFace_d0(const int& ifirstc0, const int& ifirstc1, const int& ifirstc2,
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

    invokeOneOfFour(SAMRAI_F77_FUNC(conreffacedoub3d0, CONREFFACEDOUB3D0),
              SAMRAI_F77_FUNC(conreffaceflot3d0, CONREFFACEFLOT3D0), 
              SAMRAI_F77_FUNC(conreffacecplx3d0, CONREFFACECPLX3D0),
              SAMRAI_F77_FUNC(conreffaceintg3d0, CONREFFACEINTG3D0),
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
void Call3dFortranFace_d1(const int& ifirstc0, const int& ifirstc1, const int& ifirstc2,
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

    invokeOneOfFour(SAMRAI_F77_FUNC(conreffacedoub3d1, CONREFFACEDOUB3D1),
              SAMRAI_F77_FUNC(conreffaceflot3d1, CONREFFACEFLOT3D1), 
              SAMRAI_F77_FUNC(conreffacecplx3d1, CONREFFACECPLX3D1),
              SAMRAI_F77_FUNC(conreffaceintg3d1, CONREFFACEINTG3D1),
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
void Call3dFortranFace_d2(const int& ifirstc0, const int& ifirstc1, const int& ifirstc2,
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

    invokeOneOfFour(SAMRAI_F77_FUNC(conreffacedoub3d2, CONREFFACEDOUB3D2),
              SAMRAI_F77_FUNC(conreffaceflot3d2, CONREFFACEFLOT3D2), 
              SAMRAI_F77_FUNC(conreffacecplx3d2, CONREFFACECPLX3D2),
              SAMRAI_F77_FUNC(conreffaceintg3d2, CONREFFACEINTG3D2),
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
 * Class FaceConstantRefine implements constant
 * interpolation for face-centered double patch data defined over a
 * mesh.  It is derived from the hier::RefineOperator base class.
 * The numerical operations for interpolation use FORTRAN numerical routines.
 *
 * @see hier::RefineOperator
 */

template<typename T> // one of double,float,dcomplex,int
class FaceConstantRefine:
   public hier::RefineOperator
{
public:
   /**
    * Uninteresting default constructor.
    */
   FaceConstantRefine():
      hier::RefineOperator("CONSTANT_REFINE")
   {
   }

   /**
    * Uninteresting destructor.
    */
   virtual ~FaceConstantRefine()
   {
   }

   /**
    * The priority of face-centered double constant interpolation is 0.
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
    * component on the fine patch using the face-centered double constant
    * interpolation operator.  Interpolation is performed on the intersection
    * of the destination patch and the boxes contained in fine_overlap.
    * It is assumed that the coarse patch contains sufficient data for the
    * stencil width of the refinement operator.
    *
    * @pre dynamic_cast<const FaceOverlap *>(&fine_overlap) != 0
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
#include "SAMRAI/pdat/FaceConstantRefine.cpp"
#endif

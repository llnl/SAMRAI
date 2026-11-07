/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Constant refine operator for edge-centered double data on
 *                a  mesh.
 *
 ************************************************************************/

#ifndef included_pdat_EdgeConstantRefine
#define included_pdat_EdgeConstantRefine

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
void SAMRAI_F77_FUNC(conrefedgedoub1d, CONREFEDGEDOUB1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const double *, double *);
void SAMRAI_F77_FUNC(conrefedgeflot1d, CONREFEDGEFLOT1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const float *, float *);
void SAMRAI_F77_FUNC(conrefedgecplx1d, CONREFEDGECPLX1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);
void SAMRAI_F77_FUNC(conrefedgeintg1d, CONREFEDGEINTG1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const int *, int *);

// in conrefine2d.f:
// 2d0
void SAMRAI_F77_FUNC(conrefedgedoub2d0, CONREFEDGEDOUB2D0) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const double *, double *);
void SAMRAI_F77_FUNC(conrefedgeflot2d0, CONREFEDGEFLOT2D0) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const float *, float *);
void SAMRAI_F77_FUNC(conrefedgecplx2d0, CONREFEDGECPLX2D0) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);
void SAMRAI_F77_FUNC(conrefedgeintg2d0, CONREFEDGEINTG2D0) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const int *, int *);

//2d1
void SAMRAI_F77_FUNC(conrefedgedoub2d1, CONREFEDGEDOUB2D1) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const double *, double *);
void SAMRAI_F77_FUNC(conrefedgeflot2d1, CONREFEDGEFLOT2D1) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const float *, float *);
void SAMRAI_F77_FUNC(conrefedgecplx2d1, CONREFEDGECPLX2D1) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);
void SAMRAI_F77_FUNC(conrefedgeintg2d1, CONREFEDGEINTG2D1) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const int *, int *);
// in conrefine3d.f:
// 3d0
void SAMRAI_F77_FUNC(conrefedgedoub3d0, CONREFEDGEDOUB3D0) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conrefedgeflot3d0, CONREFEDGEflot3D0) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conrefedgecplx3d0, CONREFEDGECPLX3D0) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conrefedgeintg3d0, CONREFEDGEINTG3D0) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conrefedgedoub3d1, CONREFEDGEDOUB3D1) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conrefedgeflot3d1, CONREFEDGEFLOT3D1) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conrefedgecplx3d1, CONREFEDGECPLX3D1) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conrefedgeintg3d1, CONREFEDGEINTG3D1) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conrefedgedoub3d2, CONREFEDGEDOUB3D2) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conrefedgeflot3d2, CONREFEDGEFLOT3D2) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conrefedgecplx3d2, CONREFEDGECPLX3D2) (const int&, const int&,
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
void SAMRAI_F77_FUNC(conrefedgeintg3d2, CONREFEDGEINTG3D2) (const int&, const int&,
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
void Call1dFortranEdge(const int& ifirstc, const int& ilastc,
   const int& ifirstf, const int& ilastf,
   const int& cilo, const int& cihi,
   const int& filo, const int& fihi,
   const int *ratio,
   const T *arrayc, T *arrayf)
{
   invokeOneOfFour(SAMRAI_F77_FUNC(conrefedgedoub1d, CONREFEDGEDOUB1D),
              SAMRAI_F77_FUNC(conrefedgeflot1d, CONREFEDGEFLOT1D), 
              SAMRAI_F77_FUNC(conrefedgecplx1d, CONREFEDGECPLX1D),
              SAMRAI_F77_FUNC(conrefedgeintg1d, CONREFEDGEINTG1D),
         ifirstc, ilastc, ifirstf, ilastf,
         cilo, cihi,filo, fihi,
         ratio,
         arrayc,arrayf);

}

template<typename T>
void Call2dFortranEdge_d0(const int& ifirstc0, const int& ifirstc1,
   const int& ilastc0, const int& ilastc1,
   const int& ifirstf0, const int& ifirstf1, const int& ilastf0, const int& ilastf1,
   const int& cilo0, const int& cilo1, const int& cihi0, const int& cihi1,
   const int& filo0, const int& filo1, const int& fihi0, const int& fihi1,
   const int *ratio,
   const T *arrayc, T *arrayf)
{
    invokeOneOfFour(SAMRAI_F77_FUNC(conrefedgedoub2d0, CONREFEDGEDOUB2D0),
              SAMRAI_F77_FUNC(conrefedgeflot2d0, CONREFEDGEFLOT2D0), 
              SAMRAI_F77_FUNC(conrefedgecplx2d0, CONREFEDGECPLX2D0),
              SAMRAI_F77_FUNC(conrefedgeintg2d0, CONREFEDGEINTG2D0),
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
void Call2dFortranEdge_d1(const int& ifirstc0, const int& ifirstc1,
   const int& ilastc0, const int& ilastc1,
   const int& ifirstf0, const int& ifirstf1, const int& ilastf0, const int& ilastf1,
   const int& cilo0, const int& cilo1, const int& cihi0, const int& cihi1,
   const int& filo0, const int& filo1, const int& fihi0, const int& fihi1,
   const int *ratio,
   const T *arrayc, T *arrayf)
{
    invokeOneOfFour(SAMRAI_F77_FUNC(conrefedgedoub2d1, CONREFEDGEDOUB2D1),
              SAMRAI_F77_FUNC(conrefedgeflot2d1, CONREFEDGEFLOT2D1), 
              SAMRAI_F77_FUNC(conrefedgecplx2d1, CONREFEDGECPLX2D1),
              SAMRAI_F77_FUNC(conrefedgeintg2d1, CONREFEDGEINTG2D1),
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
void Call3dFortranEdge_d0(const int& ifirstc0, const int& ifirstc1, const int& ifirstc2,
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

    invokeOneOfFour(SAMRAI_F77_FUNC(conrefedgedoub3d0, CONREFEDGEDOUB3D0),
              SAMRAI_F77_FUNC(conrefedgeflot3d0, CONREFEDGEFLOT3D0), 
              SAMRAI_F77_FUNC(conrefedgecplx3d0, CONREFEDGECPLX3D0),
              SAMRAI_F77_FUNC(conrefedgeintg3d0, CONREFEDGEINTG3D0),
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
void Call3dFortranEdge_d1(const int& ifirstc0, const int& ifirstc1, const int& ifirstc2,
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

    invokeOneOfFour(SAMRAI_F77_FUNC(conrefedgedoub3d1, CONREFEDGEDOUB3D1),
              SAMRAI_F77_FUNC(conrefedgeflot3d1, CONREFEDGEFLOT3D1), 
              SAMRAI_F77_FUNC(conrefedgecplx3d1, CONREFEDGECPLX3D1),
              SAMRAI_F77_FUNC(conrefedgeintg3d1, CONREFEDGEINTG3D1),
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
void Call3dFortranEdge_d2(const int& ifirstc0, const int& ifirstc1, const int& ifirstc2,
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

    invokeOneOfFour(SAMRAI_F77_FUNC(conrefedgedoub3d2, CONREFEDGEDOUB3D2),
              SAMRAI_F77_FUNC(conrefedgeflot3d2, CONREFEDGEFLOT3D2), 
              SAMRAI_F77_FUNC(conrefedgecplx3d2, CONREFEDGECPLX3D2),
              SAMRAI_F77_FUNC(conrefedgeintg3d2, CONREFEDGEINTG3D2),
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
 * Class EdgeConstantRefine implements constant
 * interpolation for edge-centered double patch data defined over a
 * mesh.  It is derived from the hier::RefineOperator base class.
 * The numerical operations for interpolation use FORTRAN numerical routines.
 *
 * @see hier::RefineOperator
 */
template<typename T> // one of double,float,dcomplex,int
class EdgeConstantRefine:
   public hier::RefineOperator
{
public:
   /**
    * Uninteresting default constructor.
    */
   EdgeConstantRefine():
      hier::RefineOperator("CONSTANT_REFINE")
   {
   }


   /**
    * Uninteresting destructor.
    */
   virtual ~EdgeConstantRefine()
   {
   }

   /**
    * The priority of edge-centered double constant interpolation is 0.
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
    * component on the fine patch using the edge-centered double constant
    * interpolation operator.  Interpolation is performed on the intersection
    * of the destination patch and the boxes contained in fine_overlap.
    * It is assumed that the coarse patch contains sufficient data for the
    * stencil width of the refinement operator.
    *
    * @pre dynamic_cast<const EdgeOverlap *>(&fine_overlap) != 0
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

} // namespace pdat
} // namespace SAMRAI
#include "SAMRAI/pdat/EdgeConstantRefine.cpp"
#endif

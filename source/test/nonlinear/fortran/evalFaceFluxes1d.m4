c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to evaluate face-centered fluxes in div-grad
c                operator in 1d.
c
define(NDIM,1)dnl
include(PDAT_FORTDIR/pdat_m4arrdim1d.i)dnl

      subroutine evalFaceFluxes1d(
     & lo0, hi0, ghostcells,
     & dew,
     & u,
     & dx,
     & gew )

c  Evaluate face-centered fluxes in div-grad operator.

      implicit none

      integer lo0

      integer hi0

      integer ghostcells

      double precision dx(0:NDIM-1)

      double precision dew(SIDE1d(lo,hi,0))
      double precision u(CELL1d(lo,hi,ghostcells))
      double precision gew(SIDE1d(lo,hi,0))

      integer i

      do i = lo0, hi0+1
         gew(i) = dew(i)*(u(i) - u(i-1))/dx(0)
      end do      

      return
      end

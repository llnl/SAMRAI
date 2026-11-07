c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to evaluate face-centered fluxes in div-grad
c                operator in 2d.
c
define(NDIM,2)dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine evalFaceFluxes2d(
     & lo0, hi0, lo1, hi1, ghostcells,
     & dew, dns,
     & u,
     & dx,
     & gew, gns )

c  Evaluate face-centered fluxes in div-grad operator.

      implicit none

      integer lo0
      integer lo1

      integer hi0
      integer hi1

      integer ghostcells

      double precision dx(0:NDIM-1)

      double precision dew(SIDE2d0(lo,hi,0))
      double precision dns(SIDE2d1(lo,hi,0))
      double precision u(CELL2d(lo,hi,ghostcells))
      double precision gew(SIDE2d0(lo,hi,0))
      double precision gns(SIDE2d1(lo,hi,0))

      integer i
      integer j

      do j = lo1, hi1
         do i = lo0, hi0+1
            gew(i,j) = dew(i,j)*(u(i,j) - u(i-1,j))/dx(0)
         end do
      end do
      
      do j = lo1, hi1+1
         do i = lo0, hi0
            gns(i,j) = dns(i,j)*(u(i,j) - u(i,j-1))/dx(1)
         end do
      end do

      return
      end

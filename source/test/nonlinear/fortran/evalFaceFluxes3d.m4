c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to evaluate face-centered fluxes in div-grad
c                operator in 3d.
c
define(NDIM,3)dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine evalFaceFluxes3d(
     & lo0, hi0, lo1, hi1, lo2, hi2, ghostcells,
     & dew, dns, dtb,
     & u,
     & dx,
     & gew, gns, gtb )

c  Evaluate face-centered fluxes in div-grad operator.

      implicit none

      integer lo0
      integer lo1
      integer lo2

      integer hi0
      integer hi1
      integer hi2

      integer ghostcells

      double precision dx(0:NDIM-1)

      double precision dew(SIDE3d0(lo,hi,0))
      double precision dns(SIDE3d1(lo,hi,0))
      double precision dtb(SIDE3d2(lo,hi,0))
      double precision u(CELL3d(lo,hi,ghostcells))
      double precision gew(SIDE3d0(lo,hi,0))
      double precision gns(SIDE3d1(lo,hi,0))
      double precision gtb(SIDE3d2(lo,hi,0))

      integer i
      integer j
      integer k

      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0+1
               gew(i,j,k) = dew(i,j,k)*(u(i,j,k) - u(i-1,j,k))/dx(0)
            end do
         end do
      end do

      do k = lo2, hi2
         do j = lo1, hi1+1
            do i = lo0, hi0
               gns(i,j,k) = dns(i,j,k)*(u(i,j,k) - u(i,j-1,k))/dx(1)
            end do
         end do
      end do

      do k = lo2, hi2+1
         do j = lo1, hi1
            do i = lo0, hi0
               gtb(i,j,k) = dtb(i,j,k)*(u(i,j,k) - u(i,j,k-1))/dx(2)
            end do
         end do
      end do

      return
      end

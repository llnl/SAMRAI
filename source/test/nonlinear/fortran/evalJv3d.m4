c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine for analytic evaluation of Jacobian-vector
c                products for the modified Bratu problem in 3d.
c
define(NDIM,3)dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine bratujv3d(
     & lo0, hi0, lo1, hi1, lo2, hi2, ghostcells,
     & flux0, flux1, flux2,
     & lexpu,
     & v,
     & dx,
     & dt, 
     & z )

c  Analytic evaluation of Jacobian-vector products for the modified
c  Bratu problem.

      implicit none

      integer lo0
      integer lo1
      integer lo2

      integer hi0
      integer hi1
      integer hi2

      integer ghostcells

      double precision dt

      double precision dx(0:NDIM-1)

      double precision flux0(SIDE3d0(lo,hi,0))
      double precision flux1(SIDE3d1(lo,hi,0))
      double precision flux2(SIDE3d2(lo,hi,0))
      double precision lexpu(CELL3d(lo,hi,0))
      double precision v(CELL3d(lo,hi,ghostcells))
      double precision z(CELL3d(lo,hi,0))

      integer i
      integer j
      integer k

      double precision vol

      double precision one
      parameter      ( one=1.0d0 )

      vol  = dx(0)*dx(1)*dx(2)

      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0
               z(i,j,k) = ((one - dt*lexpu(i,j,k))*v(i,j,k)
     &             - dt*( (flux0(i+1,j,k) - flux0(i,j,k))/dx(0) +
     &                    (flux1(i,j+1,k) - flux1(i,j,k))/dx(1) +
     &                    (flux2(i,j,k+1) - flux2(i,j,k))/dx(2) ))
            end do
         end do
      end do

      return
      end

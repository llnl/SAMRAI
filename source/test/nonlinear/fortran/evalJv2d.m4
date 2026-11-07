c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine for analytic evaluation of Jacobian-vector
c                products for the modified Bratu problem in 2d.
c
define(NDIM,2)dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine bratujv2d(
     & lo0, hi0, lo1, hi1, ghostcells,
     & flux0, flux1,
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

      integer hi0
      integer hi1

      integer ghostcells

      double precision dt

      double precision dx(0:NDIM-1)

      double precision flux0(SIDE2d0(lo,hi,0))
      double precision flux1(SIDE2d1(lo,hi,0))
      double precision lexpu(CELL2d(lo,hi,0))
      double precision v(CELL2d(lo,hi,ghostcells))
      double precision z(CELL2d(lo,hi,0))

      integer i
      integer j

      double precision one
      parameter      ( one=1.0d0 )

      do j = lo1, hi1
         do i = lo0, hi0
            z(i,j) = (one - dt*lexpu(i,j))*v(i,j)
     &          - dt*( (flux0(i+1,j) - flux0(i,j))/dx(0) +
     &                 (flux1(i,j+1) - flux1(i,j))/dx(1) )
         end do
      end do

      return
      end

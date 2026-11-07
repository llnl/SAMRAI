c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine for analytic evaluation of Jacobian-vector
c                products for the modified Bratu problem in 1d.
c
define(NDIM,1)dnl
include(PDAT_FORTDIR/pdat_m4arrdim1d.i)dnl

      subroutine bratujv1d(
     & lo0, hi0, ghostcells,
     & flux0,
     & lexpu,
     & v,
     & dx,
     & dt, 
     & z )

c  Analytic evaluation of Jacobian-vector products for the modified
c  Bratu problem.

      implicit none

      integer lo0

      integer hi0

      integer ghostcells

      double precision dt

      double precision dx(0:NDIM-1)

      double precision flux0(SIDE1d(lo,hi,0))
      double precision lexpu(CELL1d(lo,hi,0))
      double precision v(CELL1d(lo,hi,ghostcells))
      double precision z(CELL1d(lo,hi,0))

      integer i

      double precision vol

      double precision one
      parameter      ( one=1.0d0 )

      vol  = dx(0)

      do i = lo0, hi0
         z(i) = ((one - dt*lexpu(i))*v(i)
     &       - dt*( (flux0(i+1) - flux0(i))/dx(0) ))*vol
      end do

      return
      end

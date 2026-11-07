c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to evaluate rhs of modified bratu problem at u by
c                subtracting off the terms that are due to BE time
c                discretization in 3d.
c
define(NDIM,3)dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine evalF3d(
     & lo0, hi0, lo1, hi1, lo2, hi2,
     & dx,
     & unew, ucur, 
     & f )

c Evaluate rhs of modified bratu problem at u by subtracting off the
c terms that are due to BE time discretization.

      implicit none

      integer lo0
      integer lo1
      integer lo2

      integer hi0
      integer hi1
      integer hi2

      double precision dx(0:NDIM-1)

      double precision f(CELL3d(lo,hi,0))
      double precision unew(CELL3d(lo,hi,0))
      double precision ucur(CELL3d(lo,hi,0))
     
      integer i
      integer j
      integer k

      double precision vol

      vol = dx(0)*dx(1)*dx(2)

      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0
               f(i,j,k) = f(i,j,k) - vol*(unew(i,j,k) - ucur(i,j,k))
            end do
         end do
      end do

      return
      end

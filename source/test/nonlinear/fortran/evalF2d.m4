c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to evaluate rhs of modified bratu problem at u by
c                subtracting off the terms that are due to BE time
c                discretization in 2d.
c
define(NDIM,2)dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine evalF2d(
     & lo0, hi0, lo1, hi1,
     & dx,
     & unew, ucur, 
     & f )

      implicit none

      integer lo0
      integer lo1

      integer hi0
      integer hi1

      double precision dx(0:NDIM-1)

      double precision f(CELL2d(lo,hi,0))
      double precision unew(CELL2d(lo,hi,0))
      double precision ucur(CELL2d(lo,hi,0))
     
      integer i
      integer j

      double precision vol

      vol = dx(0)*dx(1)

      do j = lo1, hi1
         do i = lo0, hi0
            f(i,j) = f(i,j) - vol*(unew(i,j) - ucur(i,j))
         end do
      end do

      return
      end

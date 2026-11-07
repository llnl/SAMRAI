c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to evaluate rhs of modified bratu problem at u by
c                subtracting off the terms that are due to BE time
c                discretization in 1d.
c
define(NDIM,1)dnl
include(PDAT_FORTDIR/pdat_m4arrdim1d.i)dnl

      subroutine evalF1d(
     & lo0, hi0,
     & dx,
     & unew, ucur, 
     & f )

      implicit none

      integer lo0
      integer hi0

      double precision dx(0:NDIM-1)

      double precision f(CELL1d(lo,hi,0))
      double precision unew(CELL1d(lo,hi,0))
      double precision ucur(CELL1d(lo,hi,0))
     
      integer i

      double precision vol

      vol = dx(0)
      do i = lo0, hi0
         f(i) = f(i) - vol*(unew(i) - ucur(i))
      end do

      return
      end

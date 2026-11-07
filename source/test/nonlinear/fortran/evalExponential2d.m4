c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to evaluate exponential term in modified Bratu
c                problem in 2d.
c
define(NDIM,2)dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine evalExponential2d(
     & lo0, hi0, lo1, hi1, 
     & u,
     & lambda,
     & lexpu )

c  Evaluate exponential term in modified Bratu problem.

      implicit none

      integer lo0
      integer lo1

      integer hi0
      integer hi1

      double precision lambda

      double precision u(CELL2d(lo,hi,0))
      double precision lexpu(CELL2d(lo,hi,0))

      integer i
      integer j

      do j = lo1, hi1
         do i = lo0, hi0
            lexpu(i,j) = lambda*exp(u(i,j))
         end do
      end do

      return
      end

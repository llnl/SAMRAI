c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to evaluate exponential term in modified Bratu
c                problem in 1d.
c
define(NDIM,1)dnl
include(PDAT_FORTDIR/pdat_m4arrdim1d.i)dnl

      subroutine evalExponential1d(
     & lo0, hi0, 
     & u,
     & lambda,
     & lexpu )

c  Evaluate exponential term in modified Bratu problem.

      implicit none

      integer lo0

      integer hi0

      double precision lambda

      double precision u(CELL1d(lo,hi,0))
      double precision lexpu(CELL1d(lo,hi,0))

      integer i

      do i = lo0, hi0
         lexpu(i) = lambda*exp(u(i))
      end do

      return
      end

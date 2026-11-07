c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to evaluate source term in modified Bratu problem
c                in 1d.
c
define(NDIM,1)dnl
include(PDAT_FORTDIR/pdat_m4arrdim1d.i)dnl

      subroutine evalSource1d(
     & lo0, hi0, 
     & lambda,
     & xlo, xhi, dx,
     & t,
     & f )

c  Evaluate source term in modified bratu problem.

      implicit none

      integer lo0

      integer hi0

      double precision lambda
      double precision t

      double precision xlo(0:NDIM-1)
      double precision xhi(0:NDIM-1)
      double precision dx(0:NDIM-1)

      double precision f(CELL1d(lo,hi,0))

      integer i

      double precision vol
      double precision xi, xterm

      double precision half,       one,       two
      parameter      ( half=0.5d0, one=1.0d0, two=2.0d0 )

      vol = dx(0)
      xi = xlo(0) + dx(0)*half
      do i = lo0, hi0
         xterm = xi*(one - xi)
         f(i) = xterm + t*two - lambda*exp(t*xterm)
         xi = xi + dx(0)
      end do

      return
      end

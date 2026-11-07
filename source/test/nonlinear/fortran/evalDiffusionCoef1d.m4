c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to evaluate face-centered diffusion coefficients
c                in 1d.
c
define(NDIM,1)dnl
include(PDAT_FORTDIR/pdat_m4arrdim1d.i)dnl

      subroutine evalDiffusionCoef1d(
     & lo0, hi0, 
     & dx, xlo, xhi,
     & dew )

c  Evaluate face-centered diffusion coefficients for modified 
c  Bratu problem.

      implicit none

      integer lo0

      integer hi0

      double precision dx(0:NDIM-1)
      double precision xlo(0:NDIM-1)
      double precision xhi(0:NDIM-1)

      double precision dew(SIDE1d(lo,hi,0))

      integer i

      double precision xi

      double precision half,       one
      parameter      ( half=0.5d0, one=1.0d0 )

      xi = xlo(0)
      do i = lo0, hi0+1
         dew(i) = one
         xi = xi + dx(0)
      end do
      
      return 
      end

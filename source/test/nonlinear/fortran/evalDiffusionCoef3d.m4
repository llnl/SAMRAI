c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to evaluate face-centered diffusion coefficients
c                in 3d.
c
define(NDIM,3)dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine evalDiffusionCoef3d(
     & lo0, hi0, lo1, hi1, lo2, hi2,
     & dx, xlo, xhi,
     & dew, dns, dtb )

c  Evaluate face-centered diffusion coefficients for modified 
c  Bratu problem.

      implicit none

      integer lo0
      integer lo1
      integer lo2

      integer hi0
      integer hi1
      integer hi2

      double precision dx(0:NDIM-1)
      double precision xlo(0:NDIM-1)
      double precision xhi(0:NDIM-1)

      double precision dew(SIDE3d0(lo,hi,0))
      double precision dns(SIDE3d1(lo,hi,0))
      double precision dtb(SIDE3d2(lo,hi,0))

      integer i
      integer j
      integer k

      double precision xi
      double precision yj
      double precision zk

      double precision half,       one
      parameter      ( half=0.5d0, one=1.0d0 )

      zk = xlo(2) + dx(2)*half
      do k = lo2, hi2
         yj = xlo(1) + dx(1)*half
         do j = lo1, hi1
            xi = xlo(0)
            do i = lo0, hi0+1
               dew(i,j,k) = one
               xi = xi + dx(0)
            end do
            yj = yj + dx(1)
         end do
         zk = zk + dx(2)
      end do

      zk = xlo(2) + dx(2)*half
      do k = lo2, hi2
         yj = xlo(1)
         do j = lo1, hi1+1
            xi = xlo(0) + dx(0)*half
            do i = lo0, hi0
               dns(i,j,k) = one
               xi = xi + dx(0)
            end do
            yj = yj + dx(1)
         end do
         zk = zk + dx(2)
      end do

      zk = xlo(2)
      do k = lo2, hi2+1
         yj = xlo(1) + dx(1)*half
         do j = lo1, hi1
            xi = xlo(0) + dx(0)*half
            do i = lo0, hi0
               dtb(i,j,k) = one
               xi = xi + dx(0)
            end do
            yj = yj + dx(1)
         end do
         zk = zk + dx(2)
      end do
      
      return 
      end

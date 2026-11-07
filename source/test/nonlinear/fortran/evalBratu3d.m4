c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to evaluate modified bratu problem at u by
c                assembling fluxes (in gew), sources (in f and lexpu), and
c                previous time step (in v) in 3d.
c
define(NDIM,3)dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine evalbratu3d(
     & lo0, hi0, lo1, hi1, lo2, hi2, ghostcells,
     & gew, gns, gtb,
     & f, lexpu,
     & v,
     & u,
     & dx, dt,
     & r )

c  Evaluate modified bratu problem at u by assembling fluxes
c  (in gew, gns, gtb), sources (in f and lexpu), and previous
c  time step (in v).

      implicit none

      integer lo0
      integer lo1
      integer lo2

      integer hi0
      integer hi1
      integer hi2

      integer ghostcells

      double precision dt

      double precision dx(0:NDIM-1)

      double precision gew(SIDE3d0(lo,hi,0))
      double precision gns(SIDE3d1(lo,hi,0))
      double precision gtb(SIDE3d2(lo,hi,0))
      double precision f(CELL3d(lo,hi,0))
      double precision lexpu(CELL3d(lo,hi,0))
      double precision v(CELL3d(lo,hi,0))
      double precision u(CELL3d(lo,hi,ghostcells))
      double precision r(CELL3d(lo,hi,0))
     
      integer i
      integer j
      integer k

      do k = lo2, hi2
         do j = lo1, hi1
            do i = lo0, hi0
               r(i,j,k) = u(i,j,k) - v(i,j,k)
     &                  - dt*((gew(i+1,j,k) - gew(i,j,k))/dx(0) +
     &                        (gns(i,j+1,k) - gns(i,j,k))/dx(1) +
     &                        (gtb(i,j,k+1) - gtb(i,j,k))/dx(2) +
     &                        lexpu(i,j,k) + f(i,j,k))
            end do
         end do
      end do

      return
      end

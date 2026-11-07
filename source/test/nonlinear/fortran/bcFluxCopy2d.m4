c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to copy face fluxes in the x and y directions.
c
define(NDIM,2)dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine ewfluxcopy2d(
     & lo0, hi0, lo1, hi1,
     & ewflux,
     & sideflux,
     & side )

c  Copy face fluxes in the x-direction.

      implicit none

      integer lo0
      integer lo1

      integer hi0
      integer hi1

      integer side

      double precision ewflux(SIDE2d0(lo,hi,0))
      double precision sideflux(OUTERSIDE2d0(lo,hi,0))

      integer ihi
      integer ilo
      integer j

      if (side .eq. 0) then
         ilo = lo0
         do j = lo1, hi1
            sideflux(j) = ewflux(ilo,j)
         end do
      else if (side .eq. 1) then
         ihi = hi0+1
         do j = lo1, hi1
            sideflux(j) = ewflux(ihi,j)
         end do
      end if

      return
      end

      subroutine nsfluxcopy2d(
     & lo0, hi0, lo1, hi1,
     & nsflux,
     & sideflux,
     & side )

c  Copy face fluxes in the y-direction.

      implicit none

      integer lo0
      integer lo1

      integer hi0
      integer hi1

      integer side

      double precision nsflux(SIDE2d1(lo,hi,0))
      double precision sideflux(OUTERSIDE2d1(lo,hi,0))

      integer i
      integer jhi
      integer jlo

      if (side .eq. 0) then
         jlo = lo1
         do i = lo0, hi0
            sideflux(i) = nsflux(i,jlo)
         end do
      else if (side .eq. 1) then
         jhi = hi1+1
         do i = lo0, hi0
            sideflux(i) = nsflux(i,jhi)
         end do
      end if

      return
      end

c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to copy face fluxes in the x, y, and z directions.
c
define(NDIM,3)dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine ewfluxcopy3d(
     & lo0, hi0, lo1, hi1, lo2, hi2,
     & ewflux,
     & sideflux,
     & side )

c  Copy face fluxes in the x-direction.

      implicit none

      integer lo0
      integer lo1
      integer lo2

      integer hi0
      integer hi1
      integer hi2

      integer side

      double precision ewflux(SIDE3d0(lo,hi,0))
      double precision sideflux(OUTERSIDE3d0(lo,hi,0))

      integer ihi
      integer ilo
      integer j
      integer k

      if (side .eq. 0) then
         ilo = lo0
         do k = lo2, hi2
            do j = lo1, hi1
               sideflux(j,k) = ewflux(ilo,j,k)
            end do
         end do
      else if (side .eq. 1) then
         ihi = hi0+1
         do k = lo2, hi2
            do j = lo1, hi1
               sideflux(j,k) = ewflux(ihi,j,k)
            end do
         end do
      end if

      return
      end

      subroutine nsfluxcopy3d(
     & lo0, hi0, lo1, hi1, lo2, hi2,
     & nsflux,
     & sideflux,
     & side )

c  Copy face fluxes in the y-direction.

      implicit none

      integer lo0
      integer lo1
      integer lo2

      integer hi0
      integer hi1
      integer hi2

      integer side

      double precision nsflux(SIDE3d1(lo,hi,0))
      double precision sideflux(OUTERSIDE3d1(lo,hi,0))

      integer i
      integer jhi
      integer jlo
      integer k

      if (side .eq. 0) then
         jlo = lo1
         do k = lo2, hi2
            do i = lo0, hi0
               sideflux(i,k) = nsflux(i,jlo,k)
            end do
         end do
      else if (side .eq. 1) then
         jhi = hi1+1
         do k = lo2, hi2
            do i = lo0, hi0
               sideflux(i,k) = nsflux(i,jhi,k)
            end do
         end do
      end if

      return
      end

      subroutine tbfluxcopy3d(
     & lo0, hi0, lo1, hi1, lo2, hi2,
     & tbflux,
     & sideflux,
     & side )

c  Copy face fluxes in the z-direction.

      implicit none

      integer lo0
      integer lo1
      integer lo2

      integer hi0
      integer hi1
      integer hi2

      integer side

      double precision tbflux(SIDE3d2(lo,hi,0))
      double precision sideflux(OUTERSIDE3d2(lo,hi,0))

      integer i
      integer j
      integer khi
      integer klo

      if (side .eq. 0) then
         klo = lo2
         do j = lo1, hi1
            do i = lo0, hi0
               sideflux(i,j) = tbflux(i,j,klo)
            end do
         end do
      else if (side .eq. 1) then
         khi = hi2+1
         do j = lo1, hi1
            do i = lo0, hi0
               sideflux(i,j) = tbflux(i,j,khi)
            end do
         end do
      end if 

      return
      end


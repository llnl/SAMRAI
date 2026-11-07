c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to fixup boundary fluxes in the x, y, and z
c                directions.
c
define(NDIM,3)dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine ewbcfluxfix3d(
     & lo0, hi0, lo1, hi1, lo2, hi2,
     & ghostcells,
     & dx,
     & u,
     & ewflux,
     & seglo, seghi, 
     & face )

c  Fixup routine for boundary fluxes in the x-direction.

      implicit none

      integer lo0
      integer lo1
      integer lo2

      integer hi0
      integer hi1
      integer hi2

      integer ghostcells

      integer segHi(0:NDIM-1)
      integer segLo(0:NDIM-1)

      integer face

      double precision dx(0:NDIM-1)
      double precision u(CELL3d(lo,hi,ghostcells))
      double precision ewflux(SIDE3d0(lo,hi,0))

      integer ihi
      integer ilo
      integer j
      integer jhi
      integer jlo
      integer k
      integer khi
      integer klo

      double precision two,        three      
      parameter      ( two= 2.0d0, three=3.0d0 )
      double precision eight,       nine
      parameter      ( eight=8.0d0, nine=9.0d0 )

      jlo = max(lo1, segLo(1))
      jhi = min(hi1, segHi(1))
      klo = max(lo2, segLo(2))
      khi = min(hi2, segHi(2))

      if (face .eq. 0) then
         ilo = lo0
         do k = klo, khi
            do j = jlo, jhi
c  This is for first-order approximations to the flux.
c               ewflux(ilo,j,k) = two*ewflux(ilo,j,k)
c  This is for second-order approximations to the flux.
               ewflux(ilo,j,k) = (-eight*u(lo0-1,j,k) + 
     &                             nine*u(lo0,j,k)    - 
     &                             u(lo0+1,j,k))/(three*dx(0))
            end do
         end do
      else if (face .eq. 1) then
         ihi = hi0+1
         do k = klo, khi
            do j = jlo, jhi
c  This is for first-order approximations to the flux.
c               ewflux(ihi,j,k) = two*ewflux(ihi,j,k)
c  This is for second-order approximations to the flux.
               ewflux(ihi,j,k) = -(-eight*u(hi0+1,j,k) +
     &                              nine*u(hi0,j,k)    -
     &                              u(hi0-1,j,k))/(three*dx(0))
            end do
         end do
      end if

      return
      end

      subroutine nsbcfluxfix3d(
     & lo0, hi0, lo1, hi1, lo2, hi2,
     & ghostcells,
     & dx,
     & u,
     & nsflux,
     & seglo, seghi, 
     & face )

c  Fixup routine for boundary fluxes in the y-direction.

      implicit none

      integer lo0
      integer lo1
      integer lo2

      integer hi0
      integer hi1
      integer hi2

      integer ghostcells

      integer segHi(0:NDIM-1)
      integer segLo(0:NDIM-1)

      integer face

      double precision dx(0:NDIM-1)
      double precision u(CELL3d(lo,hi,ghostcells))
      double precision nsflux(SIDE3d1(lo,hi,0))

      integer i
      integer ihi
      integer ilo
      integer jhi
      integer jlo
      integer k
      integer khi
      integer klo

      double precision two,        three      
      parameter      ( two= 2.0d0, three=3.0d0 )
      double precision eight,       nine
      parameter      ( eight=8.0d0, nine=9.0d0 )

      ilo = max(lo0, segLo(0))
      ihi = min(hi0, segHi(0))
      klo = max(lo2, segLo(2))
      khi = min(hi2, segHi(2))

      if (face .eq. 2) then
         jlo = lo1
         do k = klo, khi
            do i = ilo, ihi
c  This is for first-order approximations to the flux.
c               nsflux(i,jlo,k) = two*nsflux(i,jlo,k)
c  This is for second-order approximations to the flux.
               nsflux(i,jlo,k) = (-eight*u(i,lo1-1,k) +
     &                             nine*u(i,lo1,k)    -
     &                             u(i,lo1+1,k))/(three*dx(1))
            end do
         end do
      else if (face .eq. 3) then
         jhi = hi1+1
         do k = klo, khi
            do i = ilo, ihi
c  This is for first-order approximations to the flux.
c               nsflux(i,jhi,k) = two*nsflux(i,jhi,k)
c  This is for second-order approximations to the flux.
               nsflux(i,jhi,k) = - (-eight*u(i,hi1+1,k) +
     &                               nine*u(i,hi1,k)    -
     &                               u(i,hi1-1,k))/(three*dx(1))
            end do
         end do
      end if

      return
      end

      subroutine tbbcfluxfix3d(
     & lo0, hi0, lo1, hi1, lo2, hi2,
     & ghostcells,
     & dx,
     & u,
     & tbflux,
     & seglo, seghi, 
     & face )

c  Fixup routine for boundary fluxes in the z-direction.

      implicit none

      integer lo0
      integer lo1
      integer lo2

      integer hi0
      integer hi1
      integer hi2

      integer segHi(0:NDIM-1)
      integer segLo(0:NDIM-1)

      integer ghostcells

      integer face

      double precision dx(0:NDIM-1)
      double precision u(CELL3d(lo,hi,ghostcells))
      double precision tbflux(SIDE3d2(lo,hi,0))

      integer i
      integer ihi
      integer ilo
      integer j
      integer jhi
      integer jlo
      integer khi
      integer klo

      double precision two,        three      
      parameter      ( two= 2.0d0, three=3.0d0 )
      double precision eight,       nine
      parameter      ( eight=8.0d0, nine=9.0d0 )

      ilo = max(lo0, segLo(0))
      ihi = min(hi0, segHi(0))
      jlo = max(lo1, segLo(1))
      jhi = min(hi1, segHi(1))

      if (face .eq. 4) then
         klo = lo2
         do j = jlo, jhi
            do i = ilo, ihi
c  This is for first-order approximations to the flux.
c               tbflux(i,j,klo) = two*tbflux(i,j,klo)
c  This is for second-order approximations to the flux.
               tbflux(i,j,klo) = (-eight*u(i,j,lo2-1) +
     &                             nine*u(i,j,lo2)    -
     &                             u(i,j,lo2+1))/(three*dx(2))
            end do
         end do
      else if (face .eq. 5) then
         khi = hi2+1
         do j = jlo, jhi
            do i = ilo, ihi
c  This is for first-order approximations to the flux.
c               tbflux(i,j,khi) = two*tbflux(i,j,khi)
c  This is for second-order approximations to the flux.
               tbflux(i,j,khi) = -(-eight*u(i,j,hi2+1) + 
     &                              nine*u(i,j,hi2)    -
     &                              u(i,j,hi2-1))/(three*dx(2))
            end do
         end do
      end if

      return
      end

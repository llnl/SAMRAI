c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to fixup boundary fluxes in the x and y
c                directions.
c
define(NDIM,2)dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine ewbcfluxfix2d(
     & lo0, hi0, lo1, hi1,
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

      integer hi0
      integer hi1

      integer ghostcells

      integer segHi(0:NDIM-1)
      integer segLo(0:NDIM-1)

      integer face

      double precision dx(0:NDIM-1)
      double precision u(CELL2d(lo,hi,ghostcells))
      double precision ewflux(SIDE2d0(lo,hi,0))

      integer ihi
      integer ilo
      integer j
      integer jhi
      integer jlo

      double precision two,       three
      parameter      ( two=2.0d0, three=3.0d0 )
      double precision eight,       nine
      parameter      ( eight=8.0d0, nine=9.0d0 )

      jlo = max(lo1, segLo(1))
      jhi = min(hi1, segHi(1))

      if (face .eq. 0) then
         ilo = lo0
         do j = jlo, jhi
c  This is for first-order approximations to the flux.
c            ewflux(ilo,j) = two*ewflux(ilo,j)
c  This is for second-order approximations to the flux.
            ewflux(ilo,j) = (-eight*u(lo0-1,j) + 
     &                       nine*u(lo0,j)    - 
     &                       u(lo0+1,j))/(three*dx(0))
         end do
      else if (face .eq. 1) then
         ihi = hi0+1
         do j = jlo, jhi
c  This is for first-order approximations to the flux.
c            ewflux(ihi,j) = two*ewflux(ihi,j)
c  This is for second-order approximations to the flux.
            ewflux(ihi,j) = -(-eight*u(hi0+1,j) +
     &                         nine*u(hi0,j)    -
     &                         u(hi0-1,j))/(three*dx(0))
         end do
      end if

      return
      end

      subroutine nsbcfluxfix2d(
     & lo0, hi0, lo1, hi1,
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

      integer hi0
      integer hi1

      integer ghostcells

      integer segHi(0:NDIM-1)
      integer segLo(0:NDIM-1)

      integer face

      double precision dx(0:NDIM-1)
      double precision u(CELL2d(lo,hi,ghostcells))
      double precision nsflux(SIDE2d1(lo,hi,0))

      integer i
      integer ihi
      integer ilo
      integer jhi
      integer jlo

      double precision two,        three      
      parameter      ( two= 2.0d0, three=3.0d0 )
      double precision eight,       nine
      parameter      ( eight=8.0d0, nine=9.0d0 )

      ilo = max(lo0, segLo(0))
      ihi = min(hi0, segHi(0))

      if (face .eq. 2) then
         jlo = lo1
         do i = ilo, ihi
c  This is for first-order approximations to the flux.
c            nsflux(i,jlo) = two*nsflux(i,jlo)
c  This is for second-order approximations to the flux.
            nsflux(i,jlo) = (-eight*u(i,lo1-1) +
     &                        nine*u(i,lo1)    -
     &                        u(i,lo1+1))/(three*dx(1))
         end do
      else if (face .eq. 3) then
         jhi = hi1+1
         do i = ilo, ihi
c  This is for first-order approximations to the flux.
c            nsflux(i,jhi) = two*nsflux(i,jhi)
c  This is for second-order approximations to the flux.
            nsflux(i,jhi) = - (-eight*u(i,hi1+1) +
     &                          nine*u(i,hi1)    -
     &                          u(i,hi1-1))/(three*dx(1))
         end do
      end if

      return
      end

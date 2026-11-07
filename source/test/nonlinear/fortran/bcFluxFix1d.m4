c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to fixup boundary fluxes in the x-direction.
c
define(NDIM,1)dnl
include(PDAT_FORTDIR/pdat_m4arrdim1d.i)dnl

      subroutine ewbcfluxfix1d(
     & lo0, hi0, 
     & ghostcells,
     & dx,
     & u,
     & ewflux,
     & seglo, seghi, 
     & face )

c  Fixup routine for boundary fluxes in the x-direction.

      implicit none

      integer lo0

      integer hi0

      integer ghostcells

      integer segHi(0:NDIM-1)
      integer segLo(0:NDIM-1)

      integer face

      double precision dx(0:NDIM-1)
      double precision u(CELL1d(lo,hi,ghostcells))
      double precision ewflux(SIDE1d(lo,hi,0))

      integer ihi
      integer ilo

      double precision two,       three
      parameter      ( two=2.0d0, three=3.0d0 )
      double precision eight,       nine
      parameter      ( eight=8.0d0, nine=9.0d0 )

      if (face .eq. 0) then
         ilo = lo0
c  This is for first-order approximations to the flux.
c         ewflux(ilo) = two*ewflux(ilo)
c  This is for second-order approximations to the flux.
          ewflux(ilo) = (-eight*u(lo0-1) + 
     &                    nine*u(lo0)    - 
     &                    u(lo0+1))/(three*dx(0))
      else if (face .eq. 1) then
         ihi = hi0+1
c  This is for first-order approximations to the flux.
c         ewflux(ihi) = two*ewflux(ihi)
c  This is for second-order approximations to the flux.
          ewflux(ihi) = -(-eight*u(hi0+1) +
     &                    nine*u(hi0)    -
     &                    u(hi0-1))/(three*dx(0))
      end if

      return
      end

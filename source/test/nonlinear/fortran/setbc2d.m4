c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to set boundary conditions in modified Bratu
c                problem in 2d.
c
define(NDIM,2)dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine setbc2d(
     & lo0, hi0, lo1, hi1, ghostcells,
     & u,
     & bdrySegLo, bdrySegHi,
     & bdrySide )

c  Set boundary conditions in modified Bratu problem.

      implicit none

      integer lo0
      integer lo1

      integer hi0
      integer hi1

      integer ghostcells

      integer bdrySide

      integer bdrySegHi(0:NDIM-1)
      integer bdrySegLo(0:NDIM-1)

      double precision u(CELL2d(lo,hi,ghostcells))

      integer i
      integer ihi
      integer ilo
      integer j
      integer jhi
      integer jlo

      double precision zero
      parameter      ( zero=0.0d0 )

      if (bdrySide .eq. 0) then

         ilo = bdrySegLo(0)+1
         do j = bdrySegLo(1), bdrySegHi(1)
            u(ilo-1,j) = zero
         end do

      else if (bdrySide .eq. 1) then

         ihi = bdrySegHi(0)-1
         do j = bdrySegLo(1), bdrySegHi(1)
            u(ihi+1,j) = zero
         end do

      else if (bdrySide .eq. 2) then

         jlo = bdrySegLo(1)+1
         do i = bdrySegLo(0), bdrySegHi(0)
             u(i,jlo-1) = zero
         end do

      else if (bdrySide .eq. 3) then

         jhi = bdrySegHi(1)-1
         do i = bdrySegLo(0), bdrySegHi(0)
            u(i,jhi+1) = zero
         end do

      end if

      return 
      end

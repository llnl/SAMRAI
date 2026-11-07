c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines in 1d.
c
define(NDIM,1)dnl
include(PDAT_FORTDIR/pdat_m4arrdim1d.i)dnl

      subroutine compfacdiag1d(
     &  ifirst0,ilast0,
     &  dt,cellvol,
     &  expu,dsrcdu,
     &  diag)
c***********************************************************************
      implicit none
      double precision one
      parameter(one=1.d0)
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0
      double precision dt, cellvol
      double precision 
     &  expu(CELL1d(ifirst,ilast,0)),
     &  dsrcdu(CELL1d(ifirst,ilast,0))
c ouput arrays:
      double precision 
     &  diag(CELL1d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer ic0

      do ic0=ifirst0,ilast0
         diag(ic0) = cellvol
     &       *(one - dt*(expu(ic0) + dsrcdu(ic0)))
      enddo

      return
      end
c
c
c
      subroutine compfacoffdiag1d(
     &  ifirst0,ilast0,
     &  dt, cellvol,
     &  sidediff0,
     &  offdiag0)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0
      double precision dt, cellvol
      double precision
     &  sidediff0(SIDE1d(ifirst,ilast,0))
c ouput arrays:
      double precision
     &  offdiag0(FACE1d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer ie0
      double precision factor

      factor = cellvol*dt
c
      do ie0=ifirst0,ilast0+1
         offdiag0(ie0) = factor*sidediff0(ie0)
      enddo
c
      return
      end

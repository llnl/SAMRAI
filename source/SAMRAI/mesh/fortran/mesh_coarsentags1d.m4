c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for coarsening 1d integer tag values.
c
define(NDIM,1)dnl
include(PDAT_FORTDIR/pdat_m4arrdim1d.i)dnl
c
c***********************************************************************
c Constant averaging for 1d cell-centered int tag data. The operation
c uses the rule that if any value on the fine mesh contained in a cell of 
c the coarse mesh is not equal to zero, the value on the coarse mesh is 
c set to one.
c***********************************************************************
c
      subroutine coarsentags1d(
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      integer zero, one
      parameter (zero=0)
      parameter (one=1)
c
      integer
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0
      integer ratio(0:NDIM-1)
      integer
     &  arrayf(CELL1d(filo,fihi,0)),
     &  arrayc(CELL1d(cilo,cihi,0))
      integer ic0,if0,ir0
c
c***********************************************************************
c

      do ic0=ifirstc0,ilastc0
         arrayc(ic0) = zero
      enddo

      do ir0=0,ratio(0)-1
         do ic0=ifirstc0,ilastc0
            if0=ic0*ratio(0)+ir0
            if (arrayf(if0) .ne. zero) 
     &         arrayc(ic0)=one
         enddo
      enddo
c
      return
      end
c

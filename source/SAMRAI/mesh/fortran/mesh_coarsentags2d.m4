c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for coarsening 2d integer tag values.
c
define(NDIM,2)dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl
c
c***********************************************************************
c Constant averaging for 2d cell-centered double data.   The operation
c uses the rule that if any value on the fine mesh contained in a cell of 
c the coarse mesh is not equal to zero, the value on the coarse mesh is 
c set to one.
c***********************************************************************
c
      subroutine coarsentags2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      integer zero, one
      parameter (zero=0)
      parameter (one=1)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:NDIM-1)
      integer
     &  arrayf(CELL2d(filo,fihi,0)),
     &  arrayc(CELL2d(cilo,cihi,0))
      integer ic0,ic1,if0,if1,ir0,ir1
c
c***********************************************************************
c
      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            arrayc(ic0,ic1)=zero
         enddo
      enddo

      do ir1=0,ratio(1)-1
         do ir0=0,ratio(0)-1
            do ic1=ifirstc1,ilastc1
               if1=ic1*ratio(1)+ir1
               do ic0=ifirstc0,ilastc0
                  if0=ic0*ratio(0)+ir0
                  if (arrayf(if0,if1) .ne. zero)
     &               arrayc(ic0,ic1)=one
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c

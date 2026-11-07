c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for edge operations in 2d.
c
define(NDIM,2)dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl
c
c***********************************************************************
c Set edge values from cell-centered values
c***********************************************************************
c
      subroutine setedges2d(
     &  ifirst0,ifirst1,
     &  ilast0,ilast1,
     &  ngcell0,ngcell1,
     &  ngedge0,ngedge1,
     &  cell,
     &  edge0,
     &  edge1)
c***********************************************************************
      implicit none
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  ngcell0,ngcell1,
     &  ngedge0,ngedge1
      double precision
     &  cell(CELL2dVECG(ifirst,ilast,ngcell)),
     &  edge0(EDGE2d0VECG(ifirst,ilast,ngedge)),
     &  edge1(EDGE2d1VECG(ifirst,ilast,ngedge))
      integer ic0,ic1
c
c***********************************************************************
c
c  edge0(ifirst0:ilast0,ifirst1:ilast1+1)
c  --x---x--
c  |   |   |
c  --x---x--
c  |   |   |
c  --x---x--

      do ic1=ifirst1,ilast1+1
         do ic0=ifirst0,ilast0
            edge0(ic0,ic1) = cell(ic0,ic1-1) + cell(ic0,ic1)
         enddo
      enddo

c  edge1(ifirst0:ilast0+1,ifirst1:ilast1)
c  ---------
c  x   x   x
c  ---------
c  x   x   x
c  ---------
c
      do ic1=ifirst1,ilast1
         do ic0=ifirst0,ilast0+1
            edge1(ic0,ic1) = cell(ic0-1,ic1) + cell(ic0,ic1)
         enddo
      enddo

      return
      end
c

c
c***********************************************************************
c Check edge values
c***********************************************************************
c
      subroutine checkedges2d(
     &  ifirst0,ifirst1,
     &  ilast0,ilast1,
     &  ngedge0,ngedge1,
     &  correct_val,
     &  all_correct,
     &  edge0,
     &  edge1)
c***********************************************************************
      implicit none
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  ngedge0,ngedge1
      double precision correct_val
      integer all_correct 
      double precision
     &  edge0(EDGE2d0VECG(ifirst,ilast,ngedge)),
     &  edge1(EDGE2d1VECG(ifirst,ilast,ngedge))
      integer ic0,ic1
c
c***********************************************************************
c
      do ic1=ifirst1,ilast1+1
         do ic0=ifirst0,ilast0
           if (edge0(ic0,ic1).ne.correct_val) then
              print*, "(i,j): ",ic0,ic1,"   incorrect:",edge0(ic0,ic1),
     &                "     correct:",correct_val
              all_correct = 0
           endif
         enddo
      enddo
c
      do ic1=ifirst1,ilast1
         do ic0=ifirst0,ilast0+1
           if (edge1(ic0,ic1).ne.correct_val) then
              print*, "(i,j): ",ic0,ic1,"   incorrect:",edge1(ic0,ic1),
     &                "     correct:",correct_val
              all_correct = 0
           endif
         enddo
      enddo
c
      return
      end
c

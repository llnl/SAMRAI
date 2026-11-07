c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for edge operations in 3d.
c
define(NDIM,3)dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl
c
c***********************************************************************
c Set edge values from cell-centered values
c***********************************************************************
c
      subroutine setedges3d(
     &  ifirst0,ifirst1,ifirst2,
     &  ilast0,ilast1,ilast2,
     &  ngcell0,ngcell1,ngcell2,
     &  ngedge0,ngedge1,ngedge2,
     &  cell,
     &  edge0,
     &  edge1,
     &  edge2)
c***********************************************************************
      implicit none
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  ngcell0,ngcell1,ngcell2,
     &  ngedge0,ngedge1,ngedge2
      double precision
     &  cell(CELL3dVECG(ifirst,ilast,ngcell)),
     &  edge0(EDGE3d0VECG(ifirst,ilast,ngedge)),
     &  edge1(EDGE3d1VECG(ifirst,ilast,ngedge)),
     &  edge2(EDGE3d2VECG(ifirst,ilast,ngedge))
      integer ic0,ic1,ic2
c
c***********************************************************************
c
c  edge0(ifirst0:ilast0,ifirst1:ilast1+1,ifirst2:ilast2+1)


      do ic2=ifirst2,ilast2+1
         do ic1=ifirst1,ilast1+1
            do ic0=ifirst0,ilast0
               edge0(ic0,ic1,ic2) = 
     &                 cell(ic0,ic1-1,ic2-1) + cell(ic0,ic1,ic2-1) +
     &                 cell(ic0,ic1-1,ic2) + cell(ic0,ic1,ic2) 
            enddo
         enddo
      enddo

c  edge1(ifirst0:ilast0+1,ifirst1:ilast1,ifirst2:ilast2+1)
c
      do ic2=ifirst2,ilast2+1
         do ic1=ifirst1,ilast1
            do ic0=ifirst0,ilast0+1
               edge1(ic0,ic1,ic2) = 
     &                 cell(ic0-1,ic1,ic2-1) + cell(ic0,ic1,ic2-1) +
     &                 cell(ic0-1,ic1,ic2) + cell(ic0,ic1,ic2) 
            enddo
         enddo
      enddo

c  edge2(ifirst0:ilast0+1,ifirst1:ilast1+1,ifirst2:ilast2)
c
      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1+1
            do ic0=ifirst0,ilast0+1
               edge2(ic0,ic1,ic2) = 
     &                 cell(ic0-1,ic1-1,ic2) + cell(ic0,ic1-1,ic2) +
     &                 cell(ic0-1,ic1,ic2) + cell(ic0,ic1,ic2) 
            enddo
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
      subroutine checkedges3d(
     &  ifirst0,ifirst1,ifirst2,
     &  ilast0,ilast1,ilast2,
     &  ngedge0,ngedge1,ngedge2,
     &  correct_val,
     &  all_correct,
     &  edge0,
     &  edge1,
     &  edge2)
c***********************************************************************
      implicit none
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,ifirst2,ilast2,
     &  ngedge0,ngedge1,ngedge2
      double precision correct_val
      integer all_correct 
      double precision
     &  edge0(EDGE3d0VECG(ifirst,ilast,ngedge)),
     &  edge1(EDGE3d1VECG(ifirst,ilast,ngedge)),
     &  edge2(EDGE3d2VECG(ifirst,ilast,ngedge))
      integer ic0,ic1,ic2
c
c***********************************************************************
c
      do ic2=ifirst2,ilast2+1
         do ic1=ifirst1,ilast1+1
            do ic0=ifirst0,ilast0
              if (edge0(ic0,ic1,ic2).ne.correct_val) then  
                 print*, "(i,j,k): ",ic0,ic1,ic2,"   incorrect:",
     &               edge0(ic0,ic1,ic2),"     correct:",correct_val
                 all_correct = 0
              endif    
            enddo
         enddo
      enddo

      do ic2=ifirst2,ilast2+1
         do ic1=ifirst1,ilast1
            do ic0=ifirst0,ilast0+1
              if (edge1(ic0,ic1,ic2).ne.correct_val) then  
                 print*, "(i,j,k): ",ic0,ic1,ic2,"   incorrect:",
     &               edge1(ic0,ic1,ic2),"     correct:",correct_val
                 all_correct = 0
              endif    
            enddo
         enddo
      enddo

      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1+1
            do ic0=ifirst0,ilast0+1
              if (edge2(ic0,ic1,ic2).ne.correct_val) then  
                 print*, "(i,j,k): ",ic0,ic1,ic2,"   incorrect:",
     &               edge2(ic0,ic1,ic2),"     correct:",correct_val
                 all_correct = 0
              endif    
            enddo
         enddo
      enddo
c
      return
      end
c

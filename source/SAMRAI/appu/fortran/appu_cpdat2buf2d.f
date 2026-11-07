c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   Copies data from 2D fortran array to 1D double buffer.
c

      subroutine cpddat2buf2d(
     & gidxlo0, gidxlo1,
     & bidxlo0, bidxlo1, 
     & bidxhi0, bidxhi1,
     & gidxhi0, gidxhi1, 
     & darray, buffer, bufsize)
c     =============================================================
      implicit none
      integer gidxlo0, gidxlo1, 
     &        bidxlo0, bidxlo1,
     &        bidxhi0, bidxhi1, 
     &        gidxhi0, gidxhi1,
     &        bufsize
      double precision darray(gidxlo0:gidxhi0,
     &                        gidxlo1:gidxhi1)
      double precision buffer(0:bufsize-1)
      integer in0,in1,mark 
c     =============================================================

      mark = 0 
      do in1=bidxlo1,bidxhi1
         do in0=bidxlo0,bidxhi0
            buffer(mark) = darray(in0,in1)
            mark = mark + 1
         enddo
      enddo

      return
      end


      subroutine cpfdat2buf2d(
     & gidxlo0, gidxlo1, 
     & bidxlo0, bidxlo1,
     & bidxhi0, bidxhi1,
     & gidxhi0, gidxhi1,
     & farray, buffer, bufsize)
c     =============================================================
      implicit none
      integer gidxlo0, gidxlo1,
     &        bidxlo0, bidxlo1,
     &        bidxhi0, bidxhi1,
     &        gidxhi0, gidxhi1,
     &        bufsize 
      real farray(gidxlo0:gidxhi0,
     &             gidxlo1:gidxhi1)
      double precision buffer(0:bufsize-1)
      integer in0,in1,mark 
c     =============================================================

      mark = 0 
      do in1=bidxlo1,bidxhi1
         do in0=bidxlo0,bidxhi0
            buffer(mark) = dble(farray(in0,in1))
            mark = mark + 1
         enddo
      enddo

      return
      end

      subroutine cpidat2buf2d(
     & gidxlo0, gidxlo1, 
     & bidxlo0, bidxlo1, 
     & bidxhi0, bidxhi1, 
     & gidxhi0, gidxhi1,
     & iarray, buffer, bufsize)
c     =============================================================
      implicit none
      integer gidxlo0, gidxlo1,
     &        bidxlo0, bidxlo1,
     &        bidxhi0, bidxhi1, 
     &        gidxhi0, gidxhi1,
     &        bufsize, 
     &        iarray(gidxlo0:gidxhi0,
     &               gidxlo1:gidxhi1)
      double precision buffer(0:bufsize-1)
      integer in0,in1,mark 
c     =============================================================

      mark = 0 
      do in1=bidxlo1,bidxhi1
         do in0=bidxlo0,bidxhi0
            buffer(mark) = dble(iarray(in0,in1))
            mark = mark + 1
         enddo
      enddo

      return
      end

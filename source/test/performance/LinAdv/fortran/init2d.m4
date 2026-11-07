c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for initialization in 2d.
c
define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine linadvinit2d(data_problem,dx,xlo,xhi,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,
     &  uval,
     &  nintervals,front,
     &  i_uval)
c***********************************************************************
      implicit none
include(FORTDIR/const.i)dnl
include(FORTDIR/probparams.i)dnl
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0,gcw1
      integer data_problem
      integer nintervals
      REAL front(1:nintervals)
      REAL 
     &     dx(0:NDIM-1),xlo(0:NDIM-1),xhi(0:NDIM-1)
      REAL
     &     i_uval(1:nintervals)
      REAL
     &     uval(CELL2dVECG(ifirst,ilast,gcw))
c
c***********************************************************************     
c
      integer ic0,ic1,dir,ifr
      REAL xc(0:NDIM-1)
c
c   dir 0 two linear states (L,R) indp of y,z
c   dir 1 two linear states (L,R) indp of x,z

c     write(6,*) "Inside initplane"
c     write(6,*) "data_problem= ",data_problem
c     write(6,*) "dx= ",dx(0), dx(1)
c     write(6,*) "xlo= ",xlo(0), xlo(1),", xhi = ",xhi(0), xhi(1)
c     write(6,*) "ifirst, ilast= ",ifirst0,ilast0,"and ",ifirst1,ilast1
c     call flush(6)

      dir = 0
      if (data_problem.eq.PIECEWISE_CONSTANT_X) then
         dir = 0
      else if (data_problem.eq.PIECEWISE_CONSTANT_Y) then
         dir = 1
      endif

      if (dir.eq.0) then
         ifr = 1
         do ic0=ifirst0,ilast0
            xc(0) = xlo(0) + dx(0)*(dble(ic0-ifirst0)+half)
            if (xc(dir).gt.front(ifr)) then
              ifr = ifr+1
            endif
            do ic1=ifirst1,ilast1
               uval(ic0,ic1) = i_uval(ifr)
           enddo
         enddo
      else if (dir.eq.1) then
         ifr = 1
         do ic1=ifirst1,ilast1
            xc(1) =xlo(1)+ dx(1)*(dble(ic1-ifirst1)+half)
            if (xc(dir).gt.front(ifr)) then
               ifr = ifr+1
            endif
            do ic0=ifirst0,ilast0
               uval(ic0,ic1) = i_uval(ifr)
           enddo
         enddo
      endif
c
      return
      end

c***********************************************************************
c
c    Initialization routine where we use a spherical profile 
c
c***********************************************************************
      subroutine initsphere2d(data_problem,dx,xlo,xhi,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,
     &  uval,
     &  i_uval,o_uval,
     &  center,radius)
c***********************************************************************
      implicit none
include(FORTDIR/const.i)dnl
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0,gcw1
      integer data_problem
      REAL i_uval,o_uval
      REAL radius,center(0:NDIM-1)
      REAL 
     &     dx(0:NDIM-1),xlo(0:NDIM-1),xhi(0:NDIM-1)
c variables in 2d cell indexed         
      REAL
     &     uval(CELL2dVECG(ifirst,ilast,gcw))
c
c***********************************************************************     
c
      integer ic0,ic1
      REAL xc(0:NDIM-1),x0,x1
      REAL angle

c     write(6,*) "in initsphere"
c     write(6,*) "data_problem= ",data_problem
c     write(6,*) "dx = ",(dx(i),i=0,NDIM-1)
c     write(6,*) "xlo = ",(xlo(i),i=0,NDIM-1)
c     write(6,*) "xhi = ",(xhi(i),i=0,NDIM-1)
c     write(6,*) "ce = ",(ce(i),i=0,NDIM-1)
c     write(6,*) "radius = ",radius
c     write(6,*) "ifirst0,ilast0 = ",ifirst0,ilast0
c     write(6,*) "ifirst1,ilast1 = ",ifirst1,ilast1
c

      do ic1=ifirst1,ilast1
        xc(1) = xlo(1)+dx(1)*(dble(ic1-ifirst1)+half)
        x1 = xc(1)-center(1)
        do ic0=ifirst0,ilast0
           xc(0) = xlo(0)+dx(0)*(dble(ic0-ifirst0)+half)
           x0 = xc(0)-center(0)
           if (x1.eq.zero .and. x0.eq.zero) then
              angle = zero
           else
              angle = atan2(x1,x0)
           endif
           if ((x0**2+x1**2).lt.radius**2) then
              uval(ic0,ic1) = i_uval
           else
              uval(ic0,ic1) = o_uval
           endif
         enddo
      enddo
c
      return
      end

c***********************************************************************
c
c   Sine profile
c
c***********************************************************************
 
      subroutine linadvinitsine2d(data_problem,dx,xlo,
     &  domain_xlo,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,
     &  uval,
     &  nintervals,front,
     &  i_uval,
     &  amplitude)
c***********************************************************************
      implicit none
include(FORTDIR/const.i)dnl
include(FORTDIR/probparams.i)dnl
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0,gcw1
      integer data_problem
      integer nintervals
      REAL
     &     dx(0:NDIM-1),xlo(0:NDIM-1),
     &     domain_xlo(0:NDIM-1)
      REAL front(1:nintervals)
      REAL i_uval(1:nintervals)
      REAL amplitude,period(0:NDIM-1)
c variables in 2d cell indexed
      REAL
     &     uval(CELL2dVECG(ifirst,ilast,gcw))
c
c***********************************************************************
c
      integer ic0,ic1,j,ifr
      REAL xc(0:NDIM-1),xmid(1:10)
      REAL coef(0:NDIM-1),coscoef(1:2)
c
c     write(6,*) "Inside eulerinitsine"
c     write(6,*) "data_problem= ",data_problem
c     write(6,*) "dx= ",dx(0), dx(1)
c     write(6,*) "mc= ",mc(0), mc(1)
c     write(6,*) "xlo= ",xlo(0), xlo(1),", xhi = ",xhi(0), xhi(1)
c     write(6,*) "ifirst, ilast= ",ifirst0,ilast0,ifirst1,ilast1
c     write(6,*) "gamma= ",gamma
c     call flush(6)
 
      if (data_problem.eq.SINE_CONSTANT_Y) then
         write(6,*) "Sorry, Y direction not implemented :-("
         return
      endif
 
      coef(0) = zero
      do j=1,NDIM-1
        coef(j) = two*pi/period(j)
      enddo
 
      do ic1=ifirst1,ilast1
        xc(1) = xlo(1)+dx(1)*(dble(ic1-ifirst1)+half)
        coscoef(1) = amplitude*
     &         cos((xc(1)-domain_xlo(1))*coef(1))
        do j=1,(nintervals-1)
           xmid(j) = front(j) + coscoef(1)
        enddo
        do ic0=ifirst0,ilast0
           xc(0) = xlo(0)+dx(0)*(dble(ic0-ifirst0)+half)
           ifr = 1
           do j=1,(nintervals-1)
              if( xc(0) .gt. xmid(j) ) ifr = j+1
           enddo
           uval(ic0,ic1) = i_uval(ifr)
        enddo
      enddo
 
      return
      end

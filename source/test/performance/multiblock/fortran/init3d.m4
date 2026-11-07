c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for initialization in 3d.
c
define(NDIM,3)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine linadvinit3d(data_problem,dx,xlo,xhi,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  gcw0,gcw1,gcw2,
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
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer gcw0,gcw1,gcw2
      integer data_problem
      integer nintervals
      REAL front(1:nintervals)
      REAL 
     &     dx(0:NDIM-1),xlo(0:NDIM-1),xhi(0:NDIM-1)
      REAL 
     &     i_uval(1:nintervals)
c variables in 2d cell indexed         
      REAL
     &     uval(CELL3dVECG(ifirst,ilast,gcw))
c
c***********************************************************************     
c
      integer ic0,ic1,ic2,dir,ifr
      REAL xc(0:NDIM-1)
c
c   dir 0 two linear states (L,R) indp of y,z
c   dir 1 two linear states (L,R) indp of x,z
c   dir 2 two linear states (L,R) indp of x,y

c     write(6,*) "Inside eulerinit"
c     write(6,*) "data_problem= ",data_problem
c     write(6,*) "dx= ",dx(0), dx(1), dx(2)
c     write(6,*) "xlo= ",xlo(0), xlo(1),xhi(2)
c     write(6,*) "xhi= ",xhi(0), xhi(1),xhi(2)
c     write(6,*) "ifirst= ",ifirst0,ifirst1,ifirst2
c     write(6,*) "ilast= ",ilast0,ilast1,ilast2
c     write(6,*) "gamma= ",gamma
c     call flush(6)

      dir = 0
      if (data_problem.eq.PIECEWISE_CONSTANT_X) then
         dir = 0
      else if (data_problem.eq.PIECEWISE_CONSTANT_Y) then
         dir = 1
      else if (data_problem.eq.PIECEWISE_CONSTANT_Z) then
         dir = 2
      endif

      if (dir.eq.0) then
         ifr = 1
         do ic0=ifirst0,ilast0
            xc(0) = xlo(0)+ dx(0)*(dble(ic0-ifirst0)+half)
            if (xc(dir).gt.front(ifr)) then
               ifr = ifr+1
            endif
            do ic1=ifirst1,ilast1
               do ic2=ifirst2,ilast2
                  uval(ic0,ic1,ic2) = i_uval(ifr)
               enddo
            enddo
         enddo
      else if (dir.eq.1) then
         ifr = 1
         do ic1=ifirst1,ilast1
            xc(1) = xlo(1)+ dx(1)*(dble(ic1-ifirst1)+half)
            if (xc(dir).gt.front(ifr)) then
               ifr = ifr+1
            endif
            do ic2=ifirst2,ilast2
               do ic0=ifirst0,ilast0
                  uval(ic0,ic1,ic2) = i_uval(ifr)
              enddo
           enddo
         enddo
      else if (dir.eq.2) then
         ifr = 1
         do ic2=ifirst2,ilast2
            xc(2) = xlo(2)+ dx(2)*(dble(ic2-ifirst2)+half)
            if (xc(dir).gt.front(ifr)) then
               ifr = ifr+1
            endif
            do ic1=ifirst1,ilast1
               do ic0=ifirst0,ilast0
                  uval(ic0,ic1,ic2) = i_uval(ifr)
              enddo
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
      subroutine initsphere3d(data_problem,dx,xlo,xhi,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  gcw0,gcw1,gcw2,
     &  uval,
     &  l_uval,r_uval,
     &  ce,rad)
c***********************************************************************
      implicit none
include(FORTDIR/const.i)dnl
include(FORTDIR/probparams.i)dnl
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer gcw0,gcw1,gcw2
      integer data_problem
      REAL rad,ce(0:NDIM-1)
c variables in 1d axis indexed
c
      REAL 
     &     dx(0:NDIM-1),xlo(0:NDIM-1),xhi(0:NDIM-1)
c variables in 2d cell indexed         
      REAL
     &     uval(CELL3dVECG(ifirst,ilast,gcw))
      REAL 
     &     l_uval(0:4),
     &     r_uval(0:4)
c
c***********************************************************************     
c
      integer ic0,ic1,ic2,side
      REAL xc(0:NDIM-1),x0,x1,x2
      REAL theta,phi,rad2,rad3
c

c     write(6,*) "    dx = ",dx(0),dx(1),dx(2)
c     write(6,*) "    xlo= ",xlo(0),xlo(1),xlo(2)
c     write(6,*) "    ce = ",ce(0),ce(1),ce(2)," rad ",rad
      do ic2=ifirst2,ilast2
        xc(2) = xlo(2)+dx(2)*(dble(ic2-ifirst2)+half)
        x2 = xc(2) - ce(2)
        do ic1=ifirst1,ilast1
           xc(1) = xlo(1)+dx(1)*(dble(ic1-ifirst1)+half)
           x1 = xc(1) - ce(1)
           do ic0=ifirst0,ilast0
              xc(0) = xlo(0)+dx(0)*(dble(ic0-ifirst0)+half)
              x0 = xc(0)-ce(0)
              if (x1.eq.zero .and. x0.eq.zero) then
                theta = zero
              else
                theta = atan2(x1,x0)
              endif
              rad2 = sqrt(x0**2+x1**2)
              if (rad2.eq.zero .and. x2.eq.zero) then
                phi = zero
              else
                phi  = atan2(rad2,x2)
              endif
              rad3 = sqrt(rad2**2 + x2**2)
              if (rad3.lt.rad) then
                side = 0
              else
                side = 1
              endif
              uval(ic0,ic1,ic2)   = l_uval(side)
           enddo
        enddo
      enddo
      return
      end   

c
c***********************************************************************
c
c    Sine-wave interface
c
c***********************************************************************
 
      subroutine linadvinitsine3d(data_problem,dx,xlo,
     &  domain_xlo,domain_length,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  gcw0,gcw1,gcw2,
     &  uval,
     &  nintervals,front,
     &  i_uval,
     &  amplitude,frequency)
c***********************************************************************
      implicit none
include(FORTDIR/const.i)dnl
include(FORTDIR/probparams.i)dnl
c***********************************************************************
 
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer gcw0,gcw1,gcw2
      integer data_problem
      integer nintervals
      REAL front(1:nintervals)
      REAL
     &     dx(0:NDIM-1),xlo(0:NDIM-1),
     &     domain_xlo(0:NDIM-1),domain_length(0:NDIM-1)
      REAL i_uval(1:nintervals)
      REAL amplitude,frequency(0:NDIM-1)
c variables in 3d cell indexed
      REAL
     &     uval(CELL3dVECG(ifirst,ilast,gcw))
c
c***********************************************************************
c
      integer ic0,ic1,ic2,j,ifr
      REAL xc(0:NDIM-1),xmid(1:10)
      REAL coef(0:NDIM-1),coscoef(1:2)
c
c     write(6,*) "Inside eulerinitrm"
c     write(6,*) "data_problem= ",data_problem
c     write(6,*) "dx= ",dx(0), dx(1), dx(2)
c     write(6,*) "xlo= ",xlo(0), xlo(1), xlo(2)
c     write(6,*) "domain_xlo= ",
c    &   domain_xlo(0), domain_xlo(1), domain_xlo(2)
c     write(6,*) "domain_length= ",domain_length(0),
c    &   domain_length(1), domain_length(2)
c     write(6,*) "ifirst = ",ifirst0,ifirst1,ifirst2
c     write(6,*) "ilast = ",ilast0,ilast1,ilast2
c     write(6,*) "front= ",front(1)
c     write(6,*) "amplitude= ",amplitude
c     write(6,*) "frequency= ",frequency(0),
c    &  frequency(1),frequency(2)
c     call flush(6)
 
      if (data_problem.eq.SINE_CONSTANT_Y .or.
     &    data_problem.eq.SINE_CONSTANT_Z) then
         write(6,*) "Sorry, Y and Z directions not implemented :-("
         return
      endif
 
      coef(0) = zero
      do j=1,NDIM-1
        coef(j) = two*pi*frequency(j)/domain_length(j)
      enddo
c
      do ic2=ifirst2,ilast2
         xc(2) = xlo(2)+dx(2)*(dble(ic2-ifirst2)+half)
         coscoef(2) = amplitude*cos((xc(2)-domain_xlo(2))*coef(2))
         do ic1=ifirst1,ilast1
            xc(1) = xlo(1)+dx(1)*(dble(ic1-ifirst1)+half)
            coscoef(1) =
     &         cos((xc(1)-domain_xlo(1))*coef(1))*coscoef(2)
            do j=1,(nintervals-1)
              xmid(j) = front(j) + coscoef(1)
            enddo
            do ic0=ifirst0,ilast0
               xc(0) = xlo(0)+dx(0)*(dble(ic0-ifirst0)+half)
               ifr = 1
               do j=1,(nintervals-1)
                 if( xc(0) .gt. xmid(j) ) ifr = j+1
               enddo
               uval(ic0,ic1,ic2) = i_uval(ifr)
            enddo
         enddo
      enddo
c
      return
      end

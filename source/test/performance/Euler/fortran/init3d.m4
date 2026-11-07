c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for initialization of 3d euler equations.
c
define(NDIM,3)dnl
define(NEQU,5)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine eulerinit3d(data_problem,dx,xlo,xhi,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  gcw0,gcw1,gcw2,
     &  gamma,
     &  density,velocity,pressure,
     &  nintervals,front,
     &  i_dens,i_vel,i_pres)
c***********************************************************************
      implicit none
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer gcw0,gcw1,gcw2
      integer data_problem 
      integer nintervals
      REAL front(1:nintervals)
      REAL
     &     dx(0:NDIM-1),xlo(0:NDIM-1),xhi(0:NDIM-1),gamma
      REAL i_dens(1:nintervals),
     &     i_vel(0:NDIM-1,1:nintervals),
     &     i_pres(1:nintervals)
c variables in 3d cell indexed         
      REAL
     &     density(CELL3dVECG(ifirst,ilast,gcw)),
     &     velocity(CELL3dVECG(ifirst,ilast,gcw),0:NDIM-1),
     &     pressure(CELL3dVECG(ifirst,ilast,gcw))
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
                  density(ic0,ic1,ic2) = i_dens(ifr)
                  velocity(ic0,ic1,ic2,0) = i_vel(0,ifr)
                  velocity(ic0,ic1,ic2,1) = i_vel(1,ifr)
                  velocity(ic0,ic1,ic2,2) = i_vel(2,ifr)
                  pressure(ic0,ic1,ic2) = i_pres(ifr)
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
                  density(ic0,ic1,ic2) = i_dens(ifr)
                  velocity(ic0,ic1,ic2,0) = i_vel(0,ifr)
                  velocity(ic0,ic1,ic2,1) = i_vel(1,ifr)
                  velocity(ic0,ic1,ic2,2) = i_vel(2,ifr)
                  pressure(ic0,ic1,ic2) = i_pres(ifr)
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
                  density(ic0,ic1,ic2) = i_dens(ifr)
                  velocity(ic0,ic1,ic2,0) = i_vel(0,ifr)
                  velocity(ic0,ic1,ic2,1) = i_vel(1,ifr)
                  velocity(ic0,ic1,ic2,2) = i_vel(2,ifr)
                  pressure(ic0,ic1,ic2) = i_pres(ifr)
              enddo
           enddo
         enddo
      endif
c
      return
      end   
c
c***********************************************************************
c
c    Initialization routine where we use a spherical profile 
c
c***********************************************************************
      subroutine eulerinitsphere3d(data_problem,dx,xlo,xhi,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  gcw0,gcw1,gcw2,
     &  gamma,
     &  density,velocity,pressure,
     &  i_dens,i_vel,i_pres,
     &  o_dens,o_vel,o_pres,
     &  center,radius)
c***********************************************************************
      implicit none
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer gcw0,gcw1,gcw2
      integer data_problem
      REAL i_dens,i_vel(0:NDIM-1),i_pres,
     &     o_dens,o_vel(0:NDIM-1),o_pres
      REAL center(0:NDIM-1),radius
c variables in 1d axis indexed
c
      REAL 
     &     dx(0:NDIM-1),xlo(0:NDIM-1),xhi(0:NDIM-1),gamma
c variables in 2d cell indexed         
      REAL
     &     density(CELL3dVECG(ifirst,ilast,gcw)),
     &     velocity(CELL3dVECG(ifirst,ilast,gcw),0:NDIM-1),
     &     pressure(CELL3dVECG(ifirst,ilast,gcw))
c
c**********************************************************************     
c
      integer ic0,ic1,ic2
      REAL xc(0:NDIM-1)
      REAL angle,phi,rad2,rad3,x0,x1,x2
c
c     write(6,*) "Inside eulerinitsphere"
c     write(6,*) "data_problem= ",data_problem
c     write(6,*) "dx= ",dx(0), dx(1), dx(2)
c     write(6,*) "xlo= ",xlo(0), xlo(1),xhi(2)
c     write(6,*) "xhi= ",xhi(0), xhi(1),xhi(2)
c     write(6,*) "ifirst= ",ifirst0,ifirst1,ifirst2
c     write(6,*) "ilast= ",ilast0,ilast1,ilast2
c     write(6,*) "gamma= ",gamma
c     write(6,*) "front= ",front(0),front(1),front(2),front(3)
c     write(6,*) "     = ",front(4),front(5)
c     call flush(6)

      do ic2=ifirst2,ilast2
        xc(2) = xlo(2)+dx(2)*(dble(ic2-ifirst2)+half)
        x2 = xc(2)-center(2)
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
              rad2 = sqrt(x0**2+x1**2)
              if (rad2.eq.zero .and. x2.eq.zero) then
                phi = zero
              else
                phi  = atan2(rad2,x2)
              endif
              rad3 = sqrt(rad2**2+x2**2)
              if (rad3.lt.radius) then
                 density(ic0,ic1,ic2)   = i_dens
                 velocity(ic0,ic1,ic2,0)= i_vel(0)*sin(phi)*cos(angle)
                 velocity(ic0,ic1,ic2,1)= i_vel(1)*sin(phi)*sin(angle)
                 velocity(ic0,ic1,ic2,2)= i_vel(2)*cos(phi)
                 pressure(ic0,ic1,ic2)  = i_pres
              else
                 density(ic0,ic1,ic2)   = o_dens
                 velocity(ic0,ic1,ic2,0)= o_vel(0)*sin(phi)*cos(angle)
                 velocity(ic0,ic1,ic2,1)= o_vel(1)*sin(phi)*sin(angle)
                 velocity(ic0,ic1,ic2,2)= o_vel(2)*cos(phi)
                 pressure(ic0,ic1,ic2)  = o_pres
              endif
           enddo
        enddo
      enddo
      return
      end   

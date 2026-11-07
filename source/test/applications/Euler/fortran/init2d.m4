c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines to initialize 2d euler equations.
c
define(NDIM,2)dnl
define(NEQU,4)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine eulerinit2d(data_problem,dx,xlo,xhi,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,
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
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0,gcw1
      integer data_problem
      integer nintervals
      REAL front(1:nintervals)
      REAL 
     &     dx(0:NDIM-1),xlo(0:NDIM-1),xhi(0:NDIM-1)
      REAL i_dens(1:nintervals),
     &     i_vel(0:NDIM-1,1:nintervals),
     &     i_pres(1:nintervals) 
c variables in 2d cell indexed         
      REAL
     &     density(CELL2dVECG(ifirst,ilast,gcw)),
     &     velocity(CELL2dVECG(ifirst,ilast,gcw),0:NDIM-1),
     &     pressure(CELL2dVECG(ifirst,ilast,gcw))
c
c***********************************************************************     
c
      integer ic0,ic1,dir,ifr
      REAL xc(0:NDIM-1)
      REAL gamma
c
c   dir 0 two linear states (L,R) indp of y,z
c   dir 1 two linear states (L,R) indp of x,z
    
c     write(6,*) "Inside eulerinit"
c     write(6,*) "data_problem= ",data_problem
c     write(6,*) "PIECEWISE_CONSTANT_X= ",PIECEWISE_CONSTANT_X
c     write(6,*) "STEP = ",STEP 
c     write(6,*) "PIECEWISE_CONSTANT_Y= ",PIECEWISE_CONSTANT_Y
c     write(6,*) "dx= ",dx(0), dx(1)
c     write(6,*) "xlo= ",xlo(0), xlo(1),", xhi = ",xhi(0), xhi(1)
c     write(6,*) "ifirst, ilast= ",ifirst0,ilast0,"and ",ifirst1,ilast1
c     write(6,*) "gamma= ",gamma
c     call flush(6)

      dir = 0
      if (data_problem.eq.PIECEWISE_CONSTANT_X) then
         dir = 0
      else if (data_problem.eq.STEP) then
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
               density(ic0,ic1) = i_dens(ifr)
               velocity(ic0,ic1,0) = i_vel(0,ifr)
               velocity(ic0,ic1,1) = i_vel(1,ifr)
               pressure(ic0,ic1) = i_pres(ifr)
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
               density(ic0,ic1) = i_dens(ifr)
               velocity(ic0,ic1,0) = i_vel(0,ifr)
               velocity(ic0,ic1,1) = i_vel(1,ifr)
               pressure(ic0,ic1) = i_pres(ifr)
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
      subroutine eulerinitsphere2d(data_problem,dx,xlo,xhi,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,
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
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0,gcw1
      integer data_problem
      REAL i_dens,i_vel(0:NDIM-1),i_pres,
     &     o_dens,o_vel(0:NDIM-1),o_pres
      REAL center(0:NDIM-1),radius
      REAL 
     &     dx(0:NDIM-1),xlo(0:NDIM-1),xhi(0:NDIM-1),gamma
c variables in 2d cell indexed         
      REAL
     &     density(CELL2dVECG(ifirst,ilast,gcw)),
     &     velocity(CELL2dVECG(ifirst,ilast,gcw),0:NDIM-1),
     &     pressure(CELL2dVECG(ifirst,ilast,gcw))
c
c***********************************************************************     
c
      integer ic0,ic1
      REAL xc(0:NDIM-1),x0,x1
      REAL angle
c
c     write(6,*) "Inside eulerinitsphere"
c     write(6,*) "data_problem= ",data_problem
c     write(6,*) "dx= ",dx(0), dx(1)
c     write(6,*) "xlo= ",xlo(0), xlo(1),", xhi = ",xhi(0), xhi(1)
c     write(6,*) "ifirst, ilast= ",ifirst0,ilast0,ifirst1,ilast1
c     write(6,*) "gamma= ",gamma
c     write(6,*) "radius= ",radius
c     write(6,*) "center= ",center(0),center(1)
c     call flush(6)

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
              density(ic0,ic1) = i_dens
              velocity(ic0,ic1,0) = i_vel(0)*cos(angle)
              velocity(ic0,ic1,1) = i_vel(1)*sin(angle)
              pressure(ic0,ic1) = i_pres
           else
              density(ic0,ic1) = o_dens
              velocity(ic0,ic1,0) = o_vel(0)*cos(angle)
              velocity(ic0,ic1,1) = o_vel(1)*sin(angle)
              pressure(ic0,ic1) = o_pres
           endif
c          write(6,*) "cell, state = ",ic0,ic1, density(ic0,ic1),
c    &     pressure(ic0,ic1),velocity(ic0,ic1,0),velocity(ic0,ic1,1)
c          call flush(6)
         enddo
      enddo
c
      return
      end   

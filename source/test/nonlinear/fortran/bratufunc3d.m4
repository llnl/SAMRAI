c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines to solve
c                du/dt = div( D(x,t)*div(u) ) + lambda * exp(u) + source(x,t,u)
c                in 3d.
c
define(NDIM,3)dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

c
c The following macro definitions are used in the routines below to
c define certain expressions used in the solution of the problem:
c
c    du/dt = div( D(x,t)*div(u) ) + lambda * exp(u) + source(x,t,u)
c 
c IMPORTANT NOTE: source term definition and source term derivative 
c                 definition must be consistent
c

c
c DIFFUSION COEFFICIENT DEFINITION:
c
c source_term is written as diffcoef(i0,i1,i2) = some function
c of x = (x(0), x(1), x(2)), t = time
c
define(diffusion_coefficient,`
            sidediff$1(i0,i1,i2) = 1.0d0
')dnl

c
c SOURCE TERM DEFINITION:
c
c source term is written as src(i0,i1,i2) = some function
c of x = (x(0), x(1), x(2)), t = time, and u = u(i0,i1,i2)
c
define(source_term,`
            xterm = x(0)*(1.d0 - x(0))
            yterm = x(1)*(1.d0 - x(1))
            zterm = x(2)*(1.d0 - x(2))
            src(i0,i1,i2) = (xterm*yterm*zterm
     &        + time*(2.0d0*(xterm*zterm + xterm*yterm + yterm*zterm))
     &        - lambda*exp(time*xterm*yterm*zterm))
')dnl

c
c SOURCE TERM DERIVATIVE DEFINITION:
c
c derivative of source term w.r.t. solution is written as 
c dsrcdu(ic0,ic1,ic2) = some function of x = (xc(0), xc(1), xc(2)), 
c t = time, and u = u(ic0,ic1,ic2)
c
define(source_term_derivative,`
            dsrcdu(i0,i1,i2) = 0.0d0
')dnl



      subroutine compdiffcoef3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  xlo, dx, time,
     &  sidediff0,sidediff1,sidediff2)
c***********************************************************************
      implicit none
      double precision half
      parameter(half=0.5d0)
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      double precision xlo(0:NDIM-1), dx(0:NDIM-1), time
      double precision
     &  sidediff0(SIDE3d0(ifirst,ilast,0)),
     &  sidediff1(SIDE3d1(ifirst,ilast,0)),
     &  sidediff2(SIDE3d2(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0,i1,i2
      double precision x(0:NDIM-1)
c
      do i2=ifirst2,ilast2
         x(2) = xlo(2)+dx(2)*(dble(i2-ifirst2)+half) 
         do i1=ifirst1,ilast1
            x(1) = xlo(1)+dx(1)*(dble(i1-ifirst1)+half)
            do i0=ifirst0,ilast0+1
               x(0) = xlo(0)+dx(0)*dble(i0-ifirst0)
diffusion_coefficient(0)dnl
            enddo
         enddo
      enddo
c
      do i2=ifirst2,ilast2
         x(2) = xlo(2)+dx(2)*(dble(i2-ifirst2)+half) 
         do i1=ifirst1,ilast1+1
            x(1) = xlo(1)+dx(1)*dble(i1-ifirst1)
            do i0=ifirst0,ilast0
               x(0) = xlo(0)+dx(0)*(dble(i0-ifirst0)+half)
diffusion_coefficient(1)dnl
            enddo
         enddo
      enddo
c
      do i2=ifirst2,ilast2+1
         x(2) = xlo(2)+dx(2)*dble(i2-ifirst2)
         do i1=ifirst1,ilast1
            x(1) = xlo(1)+dx(1)*(dble(i1-ifirst1)+half)
            do i0=ifirst0,ilast0
               x(0) = xlo(0)+dx(0)*(dble(i0-ifirst0)+half)
diffusion_coefficient(2)dnl
            enddo
         enddo
      enddo
c
      return
      end
c
c
c
      subroutine compexpu3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  ngcu,
     &  lambda,
     &  u,
     &  expu)
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,ngcu
      double precision lambda
      double precision u(CELL3d(ifirst,ilast,ngcu))
c ouput arrays:
      double precision expu(CELL3d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0,i1,i2

      do i2=ifirst2,ilast2
         do i1=ifirst1,ilast1
            do i0=ifirst0,ilast0
               expu(i0,i1,i2) = lambda*exp(u(i0,i1,i2))
            enddo
         enddo
      enddo

      return
      end
c
c
c
      subroutine compsrc3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  ngcu,
     &  xlo, dx, time,
     &  lambda,
     &  u,
     &  src)
c***********************************************************************
      implicit none
      double precision half
      parameter(half=0.5d0)
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,ngcu
      double precision xlo(0:NDIM-1), dx(0:NDIM-1), time, lambda
      double precision u(CELL3d(ifirst,ilast,ngcu))
c ouput arrays:
      double precision src(CELL3d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0,i1,i2
      double precision x(0:NDIM-1)
      double precision xterm,yterm,zterm

      do i2=ifirst2,ilast2
         x(2) = xlo(2)+dx(2)*(dble(i2-ifirst2)+half)
         do i1=ifirst1,ilast1
            x(1) = xlo(1)+dx(1)*(dble(i1-ifirst1)+half)
            do i0=ifirst0,ilast0
               x(0) = xlo(0)+dx(0)*(dble(i0-ifirst0)+half)
source_term()dnl
            enddo
         enddo
      enddo

      return
      end
c
c
c
      subroutine compsideflux3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  ngcu,
     &  dx,
     &  sidediff0,sidediff1,sidediff2,
     &  u,
     &  flux0,flux1,flux2)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,ngcu
      double precision dx(0:NDIM-1)
      double precision 
     &  sidediff0(SIDE3d0(ifirst,ilast,0)),
     &  sidediff1(SIDE3d1(ifirst,ilast,0)),
     &  sidediff2(SIDE3d2(ifirst,ilast,0))
      double precision u(CELL3d(ifirst,ilast,ngcu))
c ouput arrays:
      double precision
     &  flux0(SIDE3d0(ifirst,ilast,0)),
     &  flux1(SIDE3d1(ifirst,ilast,0)),
     &  flux2(SIDE3d2(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0,i1,i2
c
      do i2=ifirst2,ilast2
         do i1=ifirst1,ilast1
            do i0=ifirst0,ilast0+1
               flux0(i0,i1,i2) = sidediff0(i0,i1,i2)
     &               * (u(i0,i1,i2)-u(i0-1,i1,i2))/dx(0)
            enddo
         enddo
      enddo
c
      do i2=ifirst2,ilast2
         do i1=ifirst1,ilast1+1
            do i0=ifirst0,ilast0
               flux1(i0,i1,i2) = sidediff1(i0,i1,i2)
     &               * (u(i0,i1,i2)-u(i0,i1-1,i2))/dx(1)
            enddo
         enddo
      enddo
c
      do i2=ifirst2,ilast2+1
         do i1=ifirst1,ilast1
            do i0=ifirst0,ilast0
               flux2(i0,i1,i2) = sidediff2(i0,i1,i2)
     &               * (u(i0,i1,i2)-u(i0,i1,i2-1))/dx(2)
            enddo
         enddo
      enddo
c
      return
      end
c
c
c
      subroutine fluxbdryfix3d(
     &   ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &   ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,
     &   iside,
     &   btype,
     &   bstate,
     &   flux0,flux1,flux2)
c***********************************************************************
      implicit none
      double precision two
      parameter(two=2.0d0)
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      integer iside
      integer btype
      double precision bstate
c inout arrays:
      double precision
     &  flux0(SIDE3d0(ifirst,ilast,0)),
     &  flux1(SIDE3d1(ifirst,ilast,0)),
     &  flux2(SIDE3d2(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0,i1,i2
c
      if (iside.eq.0) then
c***********************************************************************
c x lower boundary
c***********************************************************************
         i0 = ifirst0
         if (btype.eq.0) then
c
c           dirichlet boundary
c
            do i2=max(ifirst2,ibeg2),min(ilast2,iend2)
               do i1=max(ifirst1,ibeg1),min(ilast1,iend1)
                  flux0(i0,i1,i2) = two*flux0(i0,i1,i2)
               enddo
            enddo
         else
c
c           neumann boundary
c
            do i2=max(ifirst2,ibeg2),min(ilast2,iend2)
               do i1=max(ifirst1,ibeg1),min(ilast1,iend1)
                  flux0(i0,i1,i2) = bstate
               enddo
            enddo
         endif
c
      else if (iside.eq.1) then
c***********************************************************************
c x upper boundary
c***********************************************************************
         i0 = ilast0 + 1
         if (btype.eq.0) then
c
c           dirichlet boundary
c
            do i2=max(ifirst2,ibeg2),min(ilast2,iend2)
               do i1=max(ifirst1,ibeg1),min(ilast1,iend1)
                  flux0(i0,i1,i2) = two*flux0(i0,i1,i2)
               enddo
            enddo
         else
c
c           neumann boundary
c
            do i2=max(ifirst2,ibeg2),min(ilast2,iend2)
               do i1=max(ifirst1,ibeg1),min(ilast1,iend1)
                  flux0(i0,i1,i2) = bstate
               enddo
            enddo
         endif
c
      else if (iside.eq.2) then
c***********************************************************************
c y lower boundary
c***********************************************************************
         i1 = ifirst1
         if (btype.eq.0) then
c
c           dirichlet boundary
c
            do i2=max(ifirst2,ibeg2),min(ilast2,iend2)
               do i0=max(ifirst0,ibeg0),min(ilast0,iend0)
                  flux1(i0,i1,i2) = two*flux1(i0,i1,i2)
               enddo
            enddo
         else
c
c           neumann boundary
c
            do i2=max(ifirst2,ibeg2),min(ilast2,iend2)
               do i0=max(ifirst0,ibeg0),min(ilast0,iend0)
                  flux1(i0,i1,i2) = bstate
               enddo
            enddo
         endif
c
      else if (iside.eq.3) then
c***********************************************************************
c y upper boundary
c***********************************************************************
         i1 = ilast1 + 1
         if (btype.eq.0) then
c
c           dirichlet boundary
c
            do i2=max(ifirst2,ibeg2),min(ilast2,iend2)
               do i0=max(ifirst0,ibeg0),min(ilast0,iend0)
                  flux1(i0,i1,i2) = two*flux1(i0,i1,i2)
               enddo
            enddo
         else
c
c           neumann boundary
c
            do i2=max(ifirst2,ibeg2),min(ilast2,iend2)
               do i0=max(ifirst0,ibeg0),min(ilast0,iend0)
                  flux1(i0,i1,i2) = bstate
               enddo
            enddo
         endif
c
      else if (iside.eq.4) then
c***********************************************************************
c z lower boundary
c***********************************************************************
         i2 = ifirst2
         if (btype.eq.0) then
c
c           dirichlet boundary
c
            do i1=max(ifirst1,ibeg1),min(ilast1,iend1)
               do i0=max(ifirst0,ibeg0),min(ilast0,iend0)
                  flux2(i0,i1,i2) = two*flux2(i0,i1,i2)
               enddo
            enddo
         else
c
c           neumann boundary
c
            do i1=max(ifirst1,ibeg1),min(ilast1,iend1)
               do i0=max(ifirst0,ibeg0),min(ilast0,iend0)
                  flux2(i0,i1,i2) = bstate
               enddo
            enddo
         endif
c
      else if (iside.eq.5) then
c***********************************************************************
c z upper boundary
c***********************************************************************
         i2 = ilast2 + 1
         if (btype.eq.0) then
c
c           dirichlet boundary
c
            do i1=max(ifirst1,ibeg1),min(ilast1,iend1)
               do i0=max(ifirst0,ibeg0),min(ilast0,iend0)
                  flux2(i0,i1,i2) = two*flux2(i0,i1,i2)
               enddo
            enddo
         else
c
c           neumann boundary
c
            do i1=max(ifirst1,ibeg1),min(ilast1,iend1)
               do i0=max(ifirst0,ibeg0),min(ilast0,iend0)
                  flux2(i0,i1,i2) = bstate
               enddo
            enddo
         endif
c
      endif
c
      return
      end
c
c
c
c
      subroutine fluxcopy03d(
     &   ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &   iside,
     &   flux, outerflux)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,iside 
      double precision
     &  flux(SIDE3d0(ifirst,ilast,0))
c ouput arrays:
      double precision
     &  outerflux(OUTERSIDE3d0(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0,i1,i2

      if (iside.eq.0) then
        i0 = ifirst0
      else
        i0 = ilast0+1
      endif

      do i2=ifirst2,ilast2
         do i1=ifirst1,ilast1
            outerflux(i1,i2) = flux(i0,i1,i2)
         enddo
      enddo
c
      return
      end
c
c
c
      subroutine fluxcopy13d(
     &   ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &   iside,
     &   flux, outerflux)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,iside 
      double precision
     &  flux(SIDE3d1(ifirst,ilast,0))
c ouput arrays:
      double precision
     &  outerflux(OUTERSIDE3d1(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0,i1,i2

      if (iside.eq.0) then
        i1 = ifirst1
      else
        i1 = ilast1+1
      endif

      do i2=ifirst2,ilast2
         do i0=ifirst0,ilast0
            outerflux(i0,i2) = flux(i0,i1,i2)
         enddo
      enddo
c
      return
      end
c
c
c
      subroutine fluxcopy23d(
     &   ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &   iside,
     &   flux, outerflux)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,iside
      double precision
     &  flux(SIDE3d2(ifirst,ilast,0))
c ouput arrays:
      double precision
     &  outerflux(OUTERSIDE3d2(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0,i1,i2

      if (iside.eq.0) then
        i2 = ifirst2
      else
        i2 = ilast2+1
      endif

      do i1=ifirst1,ilast1
         do i0=ifirst0,ilast0
            outerflux(i0,i1) = flux(i0,i1,i2)
         enddo
      enddo
c
      return
      end
c
c
c
      subroutine compresidual3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  ngcu,
     &  dx, dt,
     &  u_cur,
     &  u,
     &  expu,
     &  src,
     &  flux0, flux1, flux2,
     &  resid)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,ngcu
      double precision dx(0:NDIM-1), dt
      double precision
     &  u_cur(CELL3d(ifirst,ilast,0)),
     &  u(CELL3d(ifirst,ilast,ngcu)),
     &  expu(CELL3d(ifirst,ilast,0)),
     &  src(CELL3d(ifirst,ilast,0)),
     &  flux0(SIDE3d0(ifirst,ilast,0)),
     &  flux1(SIDE3d1(ifirst,ilast,0)),
     &  flux2(SIDE3d2(ifirst,ilast,0))
c ouput arrays:
      double precision
     &  resid(CELL3d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0,i1,i2
      double precision cellvol

      cellvol = dx(0)*dx(1)*dx(2)
      do i2=ifirst2,ilast2
         do i1=ifirst1,ilast1
            do i0=ifirst0,ilast0
               resid(i0,i1,i2) = (u(i0,i1,i2) - u_cur(i0,i1,i2)
     &            - dt*(  (flux0(i0+1,i1,i2) - flux0(i0,i1,i2))/dx(0)
     &                  + (flux1(i0,i1+1,i2) - flux1(i0,i1,i2))/dx(1) 
     &                  + (flux2(i0,i1,i2+1) - flux2(i0,i1,i2))/dx(2) 
     &                  + expu(i0,i1,i2) + src(i0,i1,i2)))*cellvol
            enddo
         enddo
      enddo
c
      return
      end
c
c
c
      subroutine compsrcderv3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  ngcu,
     &  xlo, dx, time,
     &  lambda,
     &  u,
     &  dsrcdu)
c***********************************************************************
      implicit none
      double precision half
      parameter(half=0.5d0)
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,ngcu
      double precision xlo(0:NDIM-1), dx(0:NDIM-1), time, lambda
      double precision u(CELL3d(ifirst,ilast,ngcu))
c ouput arrays:
      double precision dsrcdu(CELL3d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0,i1,i2
      double precision x(0:NDIM-1)

      do i2=ifirst2,ilast2
         x(2) = xlo(2)+dx(2)*(dble(i2-ifirst2)+half)
         do i1=ifirst1,ilast1
            x(1) = xlo(1)+dx(1)*(dble(i1-ifirst1)+half)
            do i0=ifirst0,ilast0
               x(0) = xlo(0)+dx(0)*(dble(i0-ifirst0)+half)
source_term_derivative()dnl
            enddo
         enddo
      enddo

      return
      end

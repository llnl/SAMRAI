c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines to solve
c                du/dt = div( D(x,t)*div(u) ) + lambda * exp(u) + source(x,t,u)
c                in 2d.
c
define(NDIM,2)dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl

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
c source_term is written as diffcoef(i0,i1) = some function
c of x = (x(0), x(1)), t = time
c
define(diffusion_coefficient,`
            sidediff$1(i0,i1) = 1.0d0
')dnl

c
c SOURCE TERM DEFINITION:
c
c source term is written as src(i0,i1) = some function
c of x = (x(0), x(1)), t = time, and u = u(i0,i1)
c
define(source_term,`
            xterm = x(0)*(1.d0 - x(0))
            yterm = x(1)*(1.d0 - x(1))
            src(i0,i1) = (xterm*yterm
     &        + time*(2.0d0*(xterm + yterm))
     &        - lambda*exp(time*xterm*yterm))
')dnl

c
c SOURCE TERM DERIVATIVE DEFINITION:
c
c derivative of source term w.r.t. solution is written as 
c dsrcdu(ic0,ic1) = some function of x = (xc(0), xc(1)), 
c t = time, and u = u(ic0,ic1)
c
define(source_term_derivative,`
            dsrcdu(i0,i1) = 0.0d0
')dnl



      subroutine compdiffcoef2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  xlo, dx, time,
     &  sidediff0,sidediff1)
c***********************************************************************
      implicit none
      double precision half
      parameter(half=0.5d0)
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      double precision xlo(0:NDIM-1), dx(0:NDIM-1), time
      double precision
     &  sidediff0(SIDE2d0(ifirst,ilast,0)),
     &  sidediff1(SIDE2d1(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0,i1
      double precision x(0:NDIM-1)
c
      do i1=ifirst1,ilast1
         x(1) = xlo(1)+dx(1)*(dble(i1-ifirst1)+half)
         do i0=ifirst0,ilast0+1
            x(0) = xlo(0)+dx(0)*dble(i0-ifirst0)
diffusion_coefficient(0)dnl
         enddo
      enddo
c
      do i1=ifirst1,ilast1+1
         x(1) = xlo(1)+dx(1)*dble(i1-ifirst1)
         do i0=ifirst0,ilast0
            x(0) = xlo(0)+dx(0)*(dble(i0-ifirst0)+half)
diffusion_coefficient(1)dnl
         enddo
      enddo
c
      return
      end
c
c
c
      subroutine compexpu2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  ngcu,
     &  lambda,
     &  u,
     &  expu)
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ngcu
      double precision lambda
      double precision u(CELL2d(ifirst,ilast,ngcu))
c ouput arrays:
      double precision expu(CELL2d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0,i1

      do i1=ifirst1,ilast1
         do i0=ifirst0,ilast0
            expu(i0,i1) = lambda*exp(u(i0,i1))
         enddo
      enddo

      return
      end
c
c
c
      subroutine compsrc2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
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
      integer ifirst0,ilast0,ifirst1,ilast1,ngcu
      double precision xlo(0:NDIM-1), dx(0:NDIM-1), time, lambda
      double precision u(CELL2d(ifirst,ilast,ngcu))
c ouput arrays:
      double precision src(CELL2d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0,i1
      double precision x(0:NDIM-1)
      double precision xterm,yterm

      do i1=ifirst1,ilast1
         x(1) = xlo(1)+dx(1)*(dble(i1-ifirst1)+half)
         do i0=ifirst0,ilast0
            x(0) = xlo(0)+dx(0)*(dble(i0-ifirst0)+half)
source_term()dnl
         enddo
      enddo

      return
      end
c
c
c
      subroutine compsideflux2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  ngcu,
     &  dx,
     &  sidediff0,sidediff1,
     &  u,
     &  flux0,flux1)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ngcu
      double precision dx(0:NDIM-1)
      double precision 
     &  sidediff0(SIDE2d0(ifirst,ilast,0)),
     &  sidediff1(SIDE2d1(ifirst,ilast,0))
      double precision u(CELL2d(ifirst,ilast,ngcu))
c ouput arrays:
      double precision
     &  flux0(SIDE2d0(ifirst,ilast,0)),
     &  flux1(SIDE2d1(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0,i1
c
      do i1=ifirst1,ilast1
         do i0=ifirst0,ilast0+1
            flux0(i0,i1) = sidediff0(i0,i1)
     &            * (u(i0,i1)-u(i0-1,i1))/dx(0)
         enddo
      enddo
c
      do i1=ifirst1,ilast1+1
         do i0=ifirst0,ilast0
            flux1(i0,i1) = sidediff1(i0,i1)
     &            * (u(i0,i1)-u(i0,i1-1))/dx(1)
         enddo
      enddo
c
      return
      end
c
c
c
      subroutine fluxbdryfix2d(
     &   ifirst0,ilast0,ifirst1,ilast1,
     &   ibeg0,iend0,ibeg1,iend1,
     &   iside,
     &   btype,
     &   bstate,
     &   flux0,flux1)
c***********************************************************************
      implicit none
      double precision two
      parameter(two=2.0d0)
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer ibeg0,iend0,ibeg1,iend1
      integer iside
      integer btype
      double precision bstate
c inout arrays:
      double precision
     &  flux0(SIDE2d0(ifirst,ilast,0)),
     &  flux1(SIDE2d1(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0,i1
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
            do i1=max(ifirst1,ibeg1),min(ilast1,iend1)
               flux0(i0,i1) = two*flux0(i0,i1)
            enddo
         else
c
c           neumann boundary
c
            do i1=max(ifirst1,ibeg1),min(ilast1,iend1)
               flux0(i0,i1) = bstate
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
            do i1=max(ifirst1,ibeg1),min(ilast1,iend1)
               flux0(i0,i1) = two*flux0(i0,i1)
            enddo
         else
c
c           neumann boundary
c
            do i1=max(ifirst1,ibeg1),min(ilast1,iend1)
               flux0(i0,i1) = bstate
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
            do i0=max(ifirst0,ibeg0),min(ilast0,iend0)
               flux1(i0,i1) = two*flux1(i0,i1)
            enddo
         else
c
c           neumann boundary
c
            do i0=max(ifirst0,ibeg0),min(ilast0,iend0)
               flux1(i0,i1) = bstate
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
            do i0=max(ifirst0,ibeg0),min(ilast0,iend0)
               flux1(i0,i1) = two*flux1(i0,i1)
            enddo
         else
c
c           neumann boundary
c
            do i0=max(ifirst0,ibeg0),min(ilast0,iend0)
               flux1(i0,i1) = bstate
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
      subroutine fluxcopy02d(
     &   ifirst0,ilast0,ifirst1,ilast1,
     &   iside,
     &   flux, outerflux)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,iside 
      double precision
     &  flux(SIDE2d0(ifirst,ilast,0))
c ouput arrays:
      double precision
     &  outerflux(OUTERSIDE2d0(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0,i1

      if (iside.eq.0) then
        i0 = ifirst0
      else
        i0 = ilast0+1
      endif

      do i1=ifirst1,ilast1
         outerflux(i1) = flux(i0,i1)
      enddo
c
      return
      end
c
c
c
      subroutine fluxcopy12d(
     &   ifirst0,ilast0,ifirst1,ilast1,
     &   iside,
     &   flux, outerflux)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,iside 
      double precision
     &  flux(SIDE2d1(ifirst,ilast,0))
c ouput arrays:
      double precision
     &  outerflux(OUTERSIDE2d1(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0,i1

      if (iside.eq.0) then
        i1 = ifirst1
      else
        i1 = ilast1+1
      endif

      do i0=ifirst0,ilast0
         outerflux(i0) = flux(i0,i1)
      enddo
c
      return
      end
c
c
c
      subroutine compresidual2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  ngcu,
     &  dx, dt,
     &  u_cur,
     &  u,
     &  expu,
     &  src,
     &  flux0, flux1,
     &  resid)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ngcu
      double precision dx(0:NDIM-1), dt
      double precision
     &  u_cur(CELL2d(ifirst,ilast,0)),
     &  u(CELL2d(ifirst,ilast,ngcu)),
     &  expu(CELL2d(ifirst,ilast,0)),
     &  src(CELL2d(ifirst,ilast,0)),
     &  flux0(SIDE2d0(ifirst,ilast,0)),
     &  flux1(SIDE2d1(ifirst,ilast,0))
c ouput arrays:
      double precision
     &  resid(CELL2d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0,i1
      double precision cellvol

      cellvol = dx(0)*dx(1)
      do i1=ifirst1,ilast1
         do i0=ifirst0,ilast0
            resid(i0,i1) = (u(i0,i1) - u_cur(i0,i1)
     &         - dt*(  (flux0(i0+1,i1) - flux0(i0,i1))/dx(0)
     &               + (flux1(i0,i1+1) - flux1(i0,i1))/dx(1) 
     &               + expu(i0,i1) + src(i0,i1)))*cellvol
         enddo
      enddo
c
      return
      end
c
c
c
      subroutine compsrcderv2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
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
      integer ifirst0,ilast0,ifirst1,ilast1,ngcu
      double precision xlo(0:NDIM-1), dx(0:NDIM-1), time, lambda
      double precision u(CELL2d(ifirst,ilast,ngcu))
c ouput arrays:
      double precision dsrcdu(CELL2d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0,i1
      double precision x(0:NDIM-1)

      do i1=ifirst1,ilast1
         x(1) = xlo(1)+dx(1)*(dble(i1-ifirst1)+half)
         do i0=ifirst0,ilast0
            x(0) = xlo(0)+dx(0)*(dble(i0-ifirst0)+half)
source_term_derivative()dnl
         enddo
      enddo

      return
      end

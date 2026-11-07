c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines to solve
c                du/dt = div( D(x,t)*div(u) ) + lambda * exp(u) + source(x,t,u)
c                in 1d.
c
define(NDIM,1)dnl
include(PDAT_FORTDIR/pdat_m4arrdim1d.i)dnl

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
c source_term is written as diffcoef(i0) = some function
c of x = (x(0)), t = time
c
define(diffusion_coefficient,`
            sidediff$1(i0) = 1.0d0
')dnl

c
c SOURCE TERM DEFINITION:
c
c source term is written as src(i0) = some function
c of x = (x(0)), t = time, and u = u(i0)
c
define(source_term,`
            xterm = x(0)*(1.d0 - x(0))
            src(i0) = (xterm
     &        + time*(2.0d0*xterm)
     &        - lambda*exp(time*xterm))
')dnl

c
c SOURCE TERM DERIVATIVE DEFINITION:
c
c derivative of source term w.r.t. solution is written as 
c dsrcdu(ic0) = some function of x = (xc(0)), 
c t = time, and u = u(ic0)
c
define(source_term_derivative,`
            dsrcdu(i0) = 0.0d0
')dnl



      subroutine compdiffcoef1d(
     &  ifirst0,ilast0,
     &  xlo, dx, time,
     &  sidediff0)
c***********************************************************************
      implicit none
      double precision half
      parameter(half=0.5d0)
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0
      double precision xlo(0:NDIM-1), dx(0:NDIM-1), time
      double precision
     &  sidediff0(SIDE1d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0
      double precision x(0:NDIM-1)
c
      do i0=ifirst0,ilast0+1
         x(0) = xlo(0)+dx(0)*dble(i0-ifirst0)
diffusion_coefficient(0)dnl
      enddo

c
      return
      end
c
c
c
      subroutine compexpu1d(
     &  ifirst0,ilast0,
     &  ngcu,
     &  lambda,
     &  u,
     &  expu)
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ngcu
      double precision lambda
      double precision u(CELL1d(ifirst,ilast,ngcu))
c ouput arrays:
      double precision expu(CELL1d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0

      do i0=ifirst0,ilast0
         expu(i0) = lambda*exp(u(i0))
      enddo

      return
      end
c
c
c
      subroutine compsrc1d(
     &  ifirst0,ilast0,
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
      integer ifirst0,ilast0,ngcu
      double precision xlo(0:NDIM-1), dx(0:NDIM-1), time, lambda
      double precision u(CELL1d(ifirst,ilast,ngcu))
c ouput arrays:
      double precision src(CELL1d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0
      double precision x(0:NDIM-1)
      double precision xterm

      do i0=ifirst0,ilast0
         x(0) = xlo(0)+dx(0)*(dble(i0-ifirst0)+half)
source_term()dnl
      enddo

      return
      end
c
c
c
      subroutine compsideflux1d(
     &  ifirst0,ilast0,
     &  ngcu,
     &  dx,
     &  sidediff0,
     &  u,
     &  flux0)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ngcu
      double precision dx(0:NDIM-1)
      double precision 
     &  sidediff0(SIDE1d(ifirst,ilast,0))
      double precision u(CELL1d(ifirst,ilast,ngcu))
c ouput arrays:
      double precision
     &  flux0(SIDE1d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0
c
      do i0=ifirst0,ilast0+1
         flux0(i0) = sidediff0(i0)
     &         * (u(i0)-u(i0-1))/dx(0)
      enddo
c
      return
      end
c
c
c
      subroutine fluxbdryfix1d(
     &   ifirst0,ilast0,
     &   ibeg0,iend0,
     &   iside,
     &   btype,
     &   bstate,
     &   flux0)
c***********************************************************************
      implicit none
      double precision two
      parameter(two=2.0d0)
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0
      integer ibeg0,iend0
      integer iside
      integer btype
      double precision bstate
c inout arrays:
      double precision
     &  flux0(SIDE1d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0
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
            flux0(i0) = two*flux0(i0)
         else
c
c           neumann boundary
c
            flux0(i0) = bstate
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
            flux0(i0) = two*flux0(i0)
         else
c
c           neumann boundary
c
            flux0(i0) = bstate
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
      subroutine fluxcopy01d(
     &   ifirst0,ilast0,
     &   iside,
     &   flux, outerflux)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,iside 
      double precision
     &  flux(SIDE1d(ifirst,ilast,0))
c ouput arrays:
      double precision
     &  outerflux(OUTERSIDE1d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0

      if (iside.eq.0) then
        i0 = ifirst0
      else
        i0 = ilast0+1
      endif

      outerflux(1) = flux(i0)
c
      return
      end
c
c
c
      subroutine compresidual1d(
     &  ifirst0,ilast0,
     &  ngcu,
     &  dx, dt,
     &  u_cur,
     &  u,
     &  expu,
     &  src,
     &  flux0, 
     &  resid)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ngcu
      double precision dx(0:NDIM-1), dt
      double precision
     &  u_cur(CELL1d(ifirst,ilast,0)),
     &  u(CELL1d(ifirst,ilast,ngcu)),
     &  expu(CELL1d(ifirst,ilast,0)),
     &  src(CELL1d(ifirst,ilast,0)),
     &  flux0(SIDE1d(ifirst,ilast,0))
c ouput arrays:
      double precision
     &  resid(CELL1d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0
      double precision cellvol

      cellvol = dx(0)
      do i0=ifirst0,ilast0
         resid(i0) = (u(i0) - u_cur(i0)
     &      - dt*(  (flux0(i0+1) - flux0(i0))/dx(0) 
     &            + expu(i0) + src(i0)))*cellvol
      enddo
c
      return
      end
c
c
c
      subroutine compsrcderv1d(
     &  ifirst0,ilast0,
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
      integer ifirst0,ilast0,ngcu
      double precision xlo(0:NDIM-1), dx(0:NDIM-1), time, lambda
      double precision u(CELL1d(ifirst,ilast,ngcu))
c ouput arrays:
      double precision dsrcdu(CELL1d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer i0
      double precision x(0:NDIM-1)

      do i0=ifirst0,ilast0
         x(0) = xlo(0)+dx(0)*(dble(i0-ifirst0)+half)
source_term_derivative()dnl
      enddo

      return
      end

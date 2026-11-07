c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   2D F77 routines.
c
define(NDIM,2)dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine comprhs2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  ngc0,ngc1,
     &  dx,
     &  y, 
     &  diff0, 
     &  diff1,
     &  rhs)
c***********************************************************************
      implicit none
      double precision one
      parameter(one=1.d0)
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer ngc0,ngc1
      double precision dx(0:NDIM-1)
      double precision 
     &  y(CELL2dVECG(ifirst,ilast,ngc)),
     &  diff0(SIDE2d0(ifirst,ilast,0)),
     &  diff1(SIDE2d1(ifirst,ilast,0))
c output arrays:
      double precision 
     &  rhs(CELL2d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer ic0,ic1
      double precision dgrade0p, dgrade0m, dgrade1p, dgrade1m

c
c  Computes RHS for 1 eqn radiation diffusion application
c 
c      dE/dt = grad dot ( D(E)grad(E) ) = RHS
c
      do ic1=ifirst1,ilast1
         do ic0=ifirst0,ilast0

c        compute  D(E)grad(E) in X and Y 
            
            dgrade0p = diff0(ic0+1,ic1)/dx(0) *
     &                 (y(ic0+1,ic1) - y(ic0,ic1))
            dgrade0m = diff0(ic0,ic1)/dx(0) *
     &                 (y(ic0,ic1) - y(ic0-1,ic1))
            dgrade1p = diff1(ic0,ic1+1)/dx(1) *
     &                 (y(ic0,ic1+1) - y(ic0,ic1))
            dgrade1m = diff1(ic0,ic1)/dx(1) *
     &                 (y(ic0,ic1) - y(ic0,ic1-1))

c        compute  RHS 

            rhs(ic0,ic1) = (dgrade0p - dgrade0m)/dx(0) +
     &                     (dgrade1p - dgrade1m)/dx(1) 

         enddo
      enddo
c
      return
      end
c
c

c
c************************************************************************
c  Subroutine setneufluxvalues sets the outerface flag and neumann flux
c  arrays.  The flag simply holds the boundary type (0-dirichlet, 
c  1-neumann) and the neumann flux array holds the value of the 
c  neumann flux.
c************************************************************************
        subroutine setneufluxvalues2d(
     &    ifirst0,ilast0,ifirst1,ilast1,
     &    bdry_type,
     &    bdry_val,      
     &    flagx0,flagx1,
     &    flagy0,flagy1,
     &    neufluxx0,neufluxx1,
     &    neufluxy0,neufluxy1)
 
c************************************************************************
        implicit none
c input arrays:
      integer
     &    ifirst0,ilast0,ifirst1,ilast1
      integer bdry_type(0:2*NDIM-1)
      double precision bdry_val(0:2*NDIM-1)
c output arrays:
      integer
     &    flagx0(OUTERFACE2d0(ifirst,ilast,0)),
     &    flagx1(OUTERFACE2d0(ifirst,ilast,0)),
     &    flagy0(OUTERFACE2d1(ifirst,ilast,0)),
     &    flagy1(OUTERFACE2d1(ifirst,ilast,0))
      double precision
     &    neufluxx0(OUTERFACE2d0(ifirst,ilast,0)),
     &    neufluxx1(OUTERFACE2d0(ifirst,ilast,0)),
     &    neufluxy0(OUTERFACE2d1(ifirst,ilast,0)),
     &    neufluxy1(OUTERFACE2d1(ifirst,ilast,0))

c
c************************************************************************
c
      integer ic0,ic1 

c  X lower & upper
      do ic1 = ifirst1, ilast1
         flagx0(ic1)    = bdry_type(0)
         flagx1(ic1)    = bdry_type(1)
         neufluxx0(ic1) = bdry_val(0)
         neufluxx1(ic1) = bdry_val(1)
      enddo

c  Y lower & upper
      do ic0 = ifirst0, ilast0
         flagy0(ic0)    = bdry_type(2)
         flagy1(ic0)    = bdry_type(3)
         neufluxy0(ic0) = bdry_val(2)
         neufluxy1(ic0) = bdry_val(3)
      enddo


      return
      end
c

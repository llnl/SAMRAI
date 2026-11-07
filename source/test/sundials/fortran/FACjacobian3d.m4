c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   3D F77 routines.
c
define(NDIM,3)dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine comprhs3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  ngc0,ngc1,ngc2,
     &  dx,
     &  y,
     &  diff0,
     &  diff1,
     &  diff2,
     &  rhs)
c***********************************************************************
      implicit none
      double precision one
      parameter(one=1.d0)
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ngc0,ngc1,ngc2
      double precision dt, dx(0:NDIM-1)
      double precision
     &  y(CELL3dVECG(ifirst,ilast,ngc)),
     &  diff0(SIDE3d0(ifirst,ilast,0)),
     &  diff1(SIDE3d1(ifirst,ilast,0)),
     &  diff2(SIDE3d2(ifirst,ilast,0))
c output arrays:
      double precision
     &  rhs(CELL3d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer ic0,ic1,ic2
      double precision dgrade0p, dgrade0m, 
     &                 dgrade1p, dgrade1m,
     &                 dgrade2p, dgrade2m

c
c  Computes RHS for 1 eqn radiation diffusion application
c
c      dE/dt = grad dot ( D(E)grad(E) ) = RHS
c
      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            do ic0=ifirst0,ilast0

c        compute  D(E)grad(E) in X, Y, and Z

            dgrade0p = diff0(ic0+1,ic1,ic2)/dx(0) *
     &                 (y(ic0+1,ic1,ic2) - y(ic0,ic1,ic2))
            dgrade0m = diff0(ic0,ic1,ic2)/dx(0) *
     &                 (y(ic0,ic1,ic2) - y(ic0-1,ic1,ic2))
            dgrade1p = diff1(ic0,ic1+1,ic2)/dx(1) *
     &                 (y(ic0,ic1+1,ic2) - y(ic0,ic1,ic2))
            dgrade1m = diff1(ic0,ic1,ic2)/dx(1) *
     &                 (y(ic0,ic1,ic2) - y(ic0,ic1-1,ic2))
            dgrade2p = diff2(ic0,ic1,ic2+1)/dx(2) *
     &                 (y(ic0,ic1,ic2+1) - y(ic0,ic1,ic2))
            dgrade2m = diff1(ic0,ic1,ic2)/dx(1) *
     &                 (y(ic0,ic1,ic2) - y(ic0,ic1,ic2-1))

c        compute  RHS

            rhs(ic0,ic1,ic2) = (dgrade0p - dgrade0m)/dx(0) +
     &                         (dgrade1p - dgrade1m)/dx(1) +
     &                         (dgrade2p - dgrade2m)/dx(2)

            enddo
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
        subroutine setneufluxvalues3d(
     &    ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &    bdry_type,
     &    bdry_val,      
     &    flagx0,flagx1,
     &    flagy0,flagy1,
     &    flagz0,flagz1,
     &    neufluxx0,neufluxx1,
     &    neufluxy0,neufluxy1,
     &    neufluxz0,neufluxz1)
 
c************************************************************************
        implicit none
c input arrays:
      integer
     &    ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer bdry_type(0:2*NDIM-1)
      double precision bdry_val(0:2*NDIM-1)
c output arrays:
      integer
     &    flagx0(OUTERFACE3d0(ifirst,ilast,0)),
     &    flagx1(OUTERFACE3d0(ifirst,ilast,0)),
     &    flagy0(OUTERFACE3d1(ifirst,ilast,0)),
     &    flagy1(OUTERFACE3d1(ifirst,ilast,0)),
     &    flagz0(OUTERFACE3d2(ifirst,ilast,0)),
     &    flagz1(OUTERFACE3d2(ifirst,ilast,0))
      double precision
     &    neufluxx0(OUTERFACE3d0(ifirst,ilast,0)),
     &    neufluxx1(OUTERFACE3d0(ifirst,ilast,0)),
     &    neufluxy0(OUTERFACE3d1(ifirst,ilast,0)),
     &    neufluxy1(OUTERFACE3d1(ifirst,ilast,0)),
     &    neufluxz0(OUTERFACE3d2(ifirst,ilast,0)),
     &    neufluxz1(OUTERFACE3d2(ifirst,ilast,0))

c
c************************************************************************
c
      integer ic0,ic1,ic2

c  X lower & upper faces
      do ic1 = ifirst1, ilast1
         do ic2 = ifirst2, ilast2
            flagx0(ic1,ic2) = bdry_type(0)
            flagx1(ic1,ic2) = bdry_type(1)
            neufluxx0(ic1,ic2) = bdry_val(0)
            neufluxx1(ic1,ic2) = bdry_val(1)
         enddo
      enddo

c  Y lower & upper faces
      do ic2 = ifirst2, ilast2
         do ic0 = ifirst0, ilast0
            flagy0(ic2,ic0) = bdry_type(2)
            flagy1(ic2,ic0) = bdry_type(3)
            neufluxy0(ic2,ic0) = bdry_val(2)
            neufluxy1(ic2,ic0) = bdry_val(3)
         enddo
      enddo

c  Z lower & upper faces
      do ic0 = ifirst0, ilast0
         do ic1 = ifirst1, ilast1
            flagz0(ic0,ic1) = bdry_type(4)
            flagz1(ic0,ic1) = bdry_type(5)
            neufluxz0(ic0,ic1) = bdry_val(4)
            neufluxz1(ic0,ic1) = bdry_val(5)
         enddo
      enddo
c
      return
      end
c

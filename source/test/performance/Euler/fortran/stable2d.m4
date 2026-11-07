c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine for computation of stable dt for 2d euler
c                equations.
c
define(NDIM,2)dnl
define(NEQU,4)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine stabledt2d(dx,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  ngc0,ngc1,
     &  gamma,density,velocity,pressure,stabdt)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(FORTDIR/const.i)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL stabdt,dx(0:NDIM-1)
      integer ifirst0,ilast0,ifirst1,ilast1,ngc0,ngc1
c
      REAL  
     &  gamma,
     &  density(CELL2dVECG(ifirst,ilast,ngc)),
     &  velocity(CELL2dVECG(ifirst,ilast,ngc),0:NDIM-1),
     &  pressure(CELL2dVECG(ifirst,ilast,ngc))
c    
      integer ic0,ic1
      integer ighoslft(0:NDIM-1),ighosrgt(0:NDIM-1)

      REAL maxspeed(0:NDIM-1),lambda
c
      ighoslft(0) = ifirst0 - ngc0
      ighoslft(1) = ifirst1 - ngc1
      ighosrgt(0) = ilast0 + ngc0
      ighosrgt(1) = ilast1 + ngc1

      maxspeed(0)=zero
      maxspeed(1)=zero

      do  ic1=ighoslft(1),ighosrgt(1)
         do  ic0=ighoslft(0),ighosrgt(0)
            lambda = 
     &       sqrt(max(zero,gamma*pressure(ic0,ic1)/density(ic0,ic1)))
            maxspeed(0) = max(maxspeed(0),
     &           abs(velocity(ic0,ic1,0))+lambda)
            maxspeed(1) = max(maxspeed(1),
     &           abs(velocity(ic0,ic1,1))+lambda)
         enddo
      enddo
      stabdt = min((dx(1)/maxspeed(1)),(dx(0)/maxspeed(0)))
      return
      end 

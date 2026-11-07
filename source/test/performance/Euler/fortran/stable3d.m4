c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine for computation of stable dt for 3d euler
c                equations.
c
define(NDIM,3)dnl
define(NEQU,5)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine stabledt3d(dx,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  ngc0,ngc1,ngc2,
     &  gamma,density,velocity,pressure,stabdt)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(FORTDIR/const.i)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL stabdt,dx(0:NDIM-1)
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  ngc0,ngc1,ngc2
c
      REAL  
     &  density(CELL3dVECG(ifirst,ilast,ngc)),
     &  velocity(CELL3dVECG(ifirst,ilast,ngc),0:NDIM-1),
     &  pressure(CELL3dVECG(ifirst,ilast,ngc))
c    
      integer ic0,ic1,ic2
      integer ighoslft(0:NDIM-1),ighosrgt(0:NDIM-1)

      REAL maxspeed(0:NDIM-1),gamma,lambda
c
      ighoslft(0) = ifirst0 - ngc0
      ighoslft(1) = ifirst1 - ngc1
      ighoslft(2) = ifirst2 - ngc2
      ighosrgt(0) = ilast0 + ngc0
      ighosrgt(1) = ilast1 + ngc1
      ighosrgt(2) = ilast2 + ngc2

      maxspeed(0)=zero
      maxspeed(1)=zero
      maxspeed(2)=zero

      do  ic2=ighoslft(2),ighosrgt(2)
         do  ic1=ighoslft(1),ighosrgt(1)
            do  ic0=ighoslft(0),ighosrgt(0)
               lambda = sqrt(max(zero,
     &              gamma*pressure(ic0,ic1,ic2)/density(ic0,ic1,ic2)))
               maxspeed(0) = max(maxspeed(0),
     &              abs(velocity(ic0,ic1,ic2,0))+lambda)
               maxspeed(1) = max(maxspeed(1),
     &              abs(velocity(ic0,ic1,ic2,1))+lambda)
               maxspeed(2) = max(maxspeed(2),
     &              abs(velocity(ic0,ic1,ic2,2))+lambda)
            enddo
         enddo
      enddo
      stabdt = min((dx(1)/maxspeed(1)),(dx(0)/maxspeed(0)))
      stabdt = min((dx(2)/maxspeed(2)),stabdt)
c     write(6,*) " dx(0),maxspeed(0)= ",dx(0),maxspeed(0)
c      write(6,*) " dx(1),maxspeed(1)= ",dx(1),maxspeed(1)
c      write(6,*) " dx(2),maxspeed(2)= ",dx(2),maxspeed(2)
c      write(6,*) "        stabdt= ",stabdt
      return
      end 

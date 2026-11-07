c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for flux compuation for 2d euler equations.
c
define(NDIM,2)dnl
define(NEQU,4)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl
include(FORTDIR/m4flux2d.i)dnl

      subroutine fluxcorrec(dt,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  dx,
     &  gamma,density,velocity,pressure,
     &  flux0,flux1,
     &  trlft0,trlft1,
     &  trrgt0,trrgt1)
c***********************************************************************
      implicit none 
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      REAL dt 
c variables in 1d axis indexed
c
      REAL 
     &     dx(0:NDIM-1)
c variables in 2d cell indexed         
      REAL
     &     density(CELL2d(ifirst,ilast,CELLG)),
     &     velocity(CELL2d(ifirst,ilast,CELLG),0:NDIM-1),
     &     pressure(CELL2d(ifirst,ilast,CELLG)),
     &     flux0(FACE2d0(ifirst,ilast,FLUXG),NEQU),
     &     flux1(FACE2d1(ifirst,ilast,FLUXG),NEQU),
     &     trlft0(FACE2d0(ifirst,ilast,FACEG),NEQU),
     &     trrgt0(FACE2d0(ifirst,ilast,FACEG),NEQU),
     &     trlft1(FACE2d1(ifirst,ilast,FACEG),NEQU),
     &     trrgt1(FACE2d1(ifirst,ilast,FACEG),NEQU)
c
c***********************************************************************     
c
      integer ic0,ic1,k
      REAL trnsvers(NEQU)
c     REAL ttvlft(NEQU),ttvrgt(NEQU)     
      REAL ttv(NEQU)
      REAL v2norm,rho,vel1,vel0,gamma,
     &     gam_min_one
     
c     write(6,*) "In fluxcorrec()"
      if (FLUXG.lt.1) then
         write(6,*) "flux ghosts < 1!"
c         call flush(6)
         stop
      endif
c     ******************************************************************
c     * complete tracing at cell edges
c     ******************************************************************
      gam_min_one = gamma-one
c   correct the 1-direction with 0-fluxes
c     write(6,*) " correct the 1-direction with 0-fluxes"
      do ic1=ifirst1-(FLUXG),ilast1+(FLUXG)
        do ic0=ifirst0-(FLUXG-1),ilast0+(FLUXG-1)
          rho  = density(ic0,ic1)
          vel0 = velocity(ic0,ic1,0)
          vel1 = velocity(ic0,ic1,1)
          v2norm= vel1**2 +vel0**2
          do k=1,NEQU
            trnsvers(k)=
     &        (flux0(ic0+1,ic1,k)-flux0(ic0,ic1,k))*half/dx(0)
          enddo

          ttv(1)= trnsvers(1) 
          ttv(2)= (trnsvers(2) - vel0*trnsvers(1))/rho
          ttv(3)= (trnsvers(3) - vel1*trnsvers(1))/rho
          ttv(NEQU)= (v2norm*half*trnsvers(1) -
     &             vel0*trnsvers(2) - vel1*trnsvers(3) +
     &             trnsvers(NEQU))*gam_min_one
          do k=1,NEQU
            trrgt1(ic1  ,ic0,k)= trrgt1(ic1  ,ic0,k) - ttv(k)
            trlft1(ic1+1,ic0,k)= trlft1(ic1+1,ic0,k) - ttv(k)
          enddo
          trrgt1(ic1  ,ic0,NEQU) = max(smallr,trrgt1(ic1  ,ic0,NEQU))
          trlft1(ic1+1,ic0,NEQU) = max(smallr,trlft1(ic1+1,ic0,NEQU))
        enddo
      enddo
c     call flush(6)     

c   correct the 0-direction with 1-fluxes
c     write(6,*) " correct the 0-direction with 1-fluxes"
      do ic0=ifirst0-(FLUXG),ilast0+(FLUXG)
        do ic1=ifirst1-(FLUXG-1),ilast1+(FLUXG-1)
          rho  = density(ic0,ic1)
          vel0 = velocity(ic0,ic1,0)
          vel1 = velocity(ic0,ic1,1)
          v2norm= vel1**2 +vel0**2
          do k=1,NEQU
            trnsvers(k)=
     &        (flux1(ic1+1,ic0,k)-flux1(ic1,ic0,k))*half/dx(1)
          enddo
c         write(6,*) "   flux1(",ic1+1,ic0,")= ",flux1(ic1+1,ic0,1),
c    &      flux1(ic1+1,ic0,2),flux1(ic1+1,ic0,3),flux1(ic1+1,ic0,4)
c         write(6,*) "   flux1(",ic1,ic0,")= ",flux1(ic1,ic0,1),
c    &      flux1(ic1,ic0,2),flux1(ic1,ic0,3),flux1(ic1,ic0,4)
c      
          ttv(1)= trnsvers(1)
          ttv(2)= (trnsvers(2) - vel0*trnsvers(1))/rho
          ttv(3)= (trnsvers(3) - vel1*trnsvers(1))/rho
          ttv(NEQU)= (v2norm*half*trnsvers(1) -
     &             vel0*trnsvers(2) - vel1*trnsvers(3) +
     &             trnsvers(NEQU))*gam_min_one
c         write(6,*) " old trrgt0(",ic0,ic1,"=",trrgt0(ic0,ic1,1),
c    &               trrgt0(ic0,ic1,2),trrgt0(ic0,ic1,3),
c    &               trrgt0(ic0,ic1,4) 
c         write(6,*) " old trlft0(",ic0+1,ic1,"=",trlft0(ic0+1,ic1,1),
c    &               trlft0(ic0+1,ic1,2),trlft0(ic0+1,ic1,3),
c    &               trlft0(ic0+1,ic1,4) 
          do k=1,NEQU
            trrgt0(ic0  ,ic1,k)= trrgt0(ic0  ,ic1,k) - ttv(k)
            trlft0(ic0+1,ic1,k)= trlft0(ic0+1,ic1,k) - ttv(k)
          enddo
          trrgt0(ic0  ,ic1,NEQU) = max(smallr,trrgt0(ic0  ,ic1,NEQU))
          trlft0(ic0+1,ic1,NEQU) = max(smallr,trlft0(ic0+1,ic1,NEQU))
c         write(6,*) " new trrgt0(",ic0,ic1,"=",trrgt0(ic0,ic1,1),
c    &               trrgt0(ic0,ic1,2),trrgt0(ic0,ic1,3),
c    &               trrgt0(ic0,ic1,4) 
c         write(6,*) " new trlft0(",ic0+1,ic1,"=",trlft0(ic0+1,ic1,1),
c    &               trlft0(ic0+1,ic1,2),trlft0(ic0+1,ic1,3),
c    &               trlft0(ic0+1,ic1,4) 
        enddo
      enddo
c      call flush(6)     
      return
      end   
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
      subroutine fluxcalculation2d(dt,extra_cell,visco,dx,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gamma,rpchoice,
     &  density,velocity,pressure,
     &  flux0,flux1,
     &  trlft0,trlft1,trrgt0,trrgt1)
     
c***********************************************************************
      implicit none
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer extra_cell,visco,rpchoice
      REAL dt,dx(0:NDIM-1)
c variables in 2d cell indexed         
      REAL
     &     density(CELL2d(ifirst,ilast,CELLG)),
     &     velocity(CELL2d(ifirst,ilast,CELLG),NDIM),
     &     pressure(CELL2d(ifirst,ilast,CELLG))
c variables in 2d side indexed         
      REAL
     &     flux0(FACE2d0(ifirst,ilast,FLUXG),NEQU),
     &     flux1(FACE2d1(ifirst,ilast,FLUXG),NEQU), 
     &     trlft0(FACE2d0(ifirst,ilast,FACEG),NEQU),
     &     trrgt0(FACE2d0(ifirst,ilast,FACEG),NEQU),
     &     trlft1(FACE2d1(ifirst,ilast,FACEG),NEQU),
     &     trrgt1(FACE2d1(ifirst,ilast,FACEG),NEQU)
c
c***********************************************************************     
c
      integer ic0,ic1,ie0,ie1,j
      REAL   stateL(NEQU),stateR(NEQU),
     &       riemst(NEQU),gamma
      REAL mom0,mom1,Hent,v2norm,vel_tan
      REAL mom0L,mom1L,enerL,mom0R,mom1R,enerR
      REAL maxeig, vcoef,vcorr(NEQU),gam_min_one
c
c variables for hllc scheme
      REAL aLsq,aRsq,keL,keR,flux(NEQU),diff(NEQU)
      REAL mfL,mfR,star(NEQU),sL,sM,sR
      REAL w,omw,hat(NEQU+1),denom
c
      if (FLUXG.lt.1) then
         write(6,*) "flux ghosts < 1!" 
c         call flush(6)
         stop
      endif
c
c***********************************************************************
c solve riemann problems for conservative flux
c  arguments: ( axis for RP, other axis, extra cells-direction)
c***********************************************************************
c      
c     write(6,*) "In fluxcalculation2d(",extra_cell,")"
c     write(6,*) "ifirst0,ilast0,ifirst1,ilast1,extra_cell",
c    &       ifirst0,ilast0,ifirst1,ilast1,extra_cell
      gam_min_one = gamma-one
riemann_solve(0,1,(extra_cell+(FLUXG-1)))dnl
c
riemann_solve(1,0,(extra_cell+(FLUXG-1)))dnl

      if (visco.eq.1) then
      write(6,*) "doing artificial viscosity"
c
artificial_viscosity1(0,1)dnl
c
artificial_viscosity1(1,0)dnl
c
      endif
c      call flush(6)     
      return
      end 
c***********************************************************************
c***********************************************************************
c***********************************************************************
  
      subroutine consdiff2d(ifirst0,ilast0,ifirst1,ilast1,dx,
     &  flux0,flux1,
     &  gamma,density,velocity,pressure)
c***********************************************************************
      implicit none
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
      integer ifirst0, ilast0,ifirst1, ilast1
      REAL dx(0:NDIM-1)
      REAL
     &     flux0(FACE2d0(ifirst,ilast,FLUXG),NEQU),
     &     flux1(FACE2d1(ifirst,ilast,FLUXG),NEQU),
     &     gamma,
     &     density(CELL2d(ifirst,ilast,CELLG)),
     &     velocity(CELL2d(ifirst,ilast,CELLG),NDIM),
     &     pressure(CELL2d(ifirst,ilast,CELLG))
c
      integer ic0,ic1,k
      REAL temp,v2norm,mom(NDIM),energy,
     &     gam_min_one
      
c***********************************************************************
c update conserved variables to full time
c note the permutation of indices in 2nd coordinate direction
c***********************************************************************
      gam_min_one = gamma-one 

      do ic1=ifirst1,ilast1
        do ic0=ifirst0,ilast0
          mom(1) = density(ic0,ic1)*velocity(ic0,ic1,1)
          mom(2) = density(ic0,ic1)*velocity(ic0,ic1,2)
          v2norm = (velocity(ic0,ic1,1)**2+velocity(ic0,ic1,2)**2)
          energy = pressure(ic0,ic1)/gam_min_one +
     &                        half*density(ic0,ic1)*v2norm
 
          density(ic0,ic1) = density(ic0,ic1)
     &          -(flux0(ic0+1,ic1,1)-flux0(ic0,ic1,1))/dx(0)
     &          -(flux1(ic1+1,ic0,1)-flux1(ic1,ic0,1))/dx(1)
          density(ic0,ic1) = max(smallr,density(ic0,ic1))
          do k=1,NDIM
            mom(k) = mom(k)
     &        -(flux0(ic0+1,ic1,k+1)-flux0(ic0,ic1,k+1))/dx(0)
     &        -(flux1(ic1+1,ic0,k+1)-flux1(ic1,ic0,k+1))/dx(1)
            velocity(ic0,ic1,k) = mom(k)/density(ic0,ic1)
          enddo
          energy = energy
     &       -(flux0(ic0+1,ic1,NEQU)-flux0(ic0,ic1,NEQU))/dx(0)
     &       -(flux1(ic1+1,ic0,NEQU)-flux1(ic1,ic0,NEQU))/dx(1)
c      
          v2norm = (velocity(ic0,ic1,1)**2+velocity(ic0,ic1,2)**2)
          temp = energy - half*density(ic0,ic1)*v2norm
          pressure(ic0,ic1) = gam_min_one*temp
          pressure(ic0,ic1) = max(smallr,pressure(ic0,ic1))
        enddo
      enddo
c
      return
      end
c***********************************************************************
include(FORTDIR/gas1d_approxrp2d.i)dnl
include(FORTDIR/gas1d_exactrp2d.i)dnl

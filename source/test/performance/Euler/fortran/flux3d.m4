c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for flux compuation for 3d euler equations.
c
define(NDIM,3)dnl
define(NEQU,5)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl
include(FORTDIR/m4flux3d.i)dnl

      subroutine fluxcorrec2d(dt,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  dx,gamma,idir,
     &  density,velocity,pressure,
     &  flux0,flux1,flux2,
     &  tracelft0,tracelft1,tracelft2,
     &  tracergt0,tracergt1,tracergt2,
     &  ttracelft0,ttracelft1,ttracelft2,
     &  ttracergt0,ttracergt1,ttracergt2)
c***********************************************************************
      implicit none
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      REAL dt 
      integer idir
c variables in 1d axis indexed
c
      REAL 
     &     dx(0:NDIM-1)
c variables in 2d cell indexed         
      REAL
     &     gamma,
     &     density(CELL3d(ifirst,ilast,CELLG)),
     &     velocity(CELL3d(ifirst,ilast,CELLG),0:NDIM-1),
     &     pressure(CELL3d(ifirst,ilast,CELLG)),
c
     &     flux0(FACE3d0(ifirst,ilast,FLUXG),NEQU), 
     &     flux1(FACE3d1(ifirst,ilast,FLUXG),NEQU), 
     &     flux2(FACE3d2(ifirst,ilast,FLUXG),NEQU), 
c
     &     tracelft0(FACE3d0(ifirst,ilast,FACEG),NEQU),
     &     tracelft1(FACE3d1(ifirst,ilast,FACEG),NEQU),
     &     tracelft2(FACE3d2(ifirst,ilast,FACEG),NEQU),
     &     tracergt0(FACE3d0(ifirst,ilast,FACEG),NEQU),
     &     tracergt1(FACE3d1(ifirst,ilast,FACEG),NEQU),
     &     tracergt2(FACE3d2(ifirst,ilast,FACEG),NEQU),
c
     &     ttracelft0(FACE3d0(ifirst,ilast,FACEG),NEQU),
     &     ttracelft1(FACE3d1(ifirst,ilast,FACEG),NEQU),
     &     ttracelft2(FACE3d2(ifirst,ilast,FACEG),NEQU),
     &     ttracergt0(FACE3d0(ifirst,ilast,FACEG),NEQU),
     &     ttracergt1(FACE3d1(ifirst,ilast,FACEG),NEQU),
     &     ttracergt2(FACE3d2(ifirst,ilast,FACEG),NEQU)
c
c***********************************************************************     
c
      integer ic0,ic1,ic2,k
      REAL trnsvers(NEQU)
c     REAL ttvlft(NEQU),ttvrgt(NEQU)     
      REAL ttv(NEQU)
      REAL v2norm,rho,vel1,vel2,vel0,gam_min_one
c     write(6,*) "In fluxcorrec2d()"
      if (FLUXG.lt.1) then
         write(6,*) "flux ghosts < 1!"
         stop
      endif
c    
      gam_min_one = gamma - one
c
c  "Forward" computation of transverse flux terms
c
      if (idir.eq.1) then
c
correc_flux2d(0,`ic1,ic2',1,`ic2,ic0',2)dnl
c
correc_flux2d(1,`ic2,ic0',0,`ic1,ic2',2)dnl
c
correc_flux2d(2,`ic0,ic1',0,`ic1,ic2',1)dnl
c
c  "Backward" computation of transverse flux terms
c
      elseif (idir.eq.-1) then
c
correc_flux2d(0,`ic1,ic2',2,`ic0,ic1',1)dnl
c
correc_flux2d(1,`ic2,ic0',2,`ic0,ic1',0)dnl
c
correc_flux2d(2,`ic0,ic1',1,`ic2,ic0',0)dnl
c
      endif
c
      return
      end   
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
      subroutine fluxcorrec3d(dt,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  dx,
     &  gamma,density,velocity,pressure,
     &  fluxa0,fluxa1,fluxa2,
     &  fluxb0,fluxb1,fluxb2,
     &  tracelft0,tracelft1,tracelft2,
     &  tracergt0,tracergt1,tracergt2)
c***********************************************************************
      implicit none
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      REAL dt 
c variables in 1d axis indexed
c
      REAL 
     &     dx(0:NDIM-1)
c variables in 2d cell indexed         
      REAL
     &     gamma,
     &     density(CELL3d(ifirst,ilast,CELLG)),
     &     velocity(CELL3d(ifirst,ilast,CELLG),0:NDIM-1),
     &     pressure(CELL3d(ifirst,ilast,CELLG)),
     &     fluxa0(FACE3d0(ifirst,ilast,FLUXG),NEQU),
     &     fluxa1(FACE3d1(ifirst,ilast,FLUXG),NEQU), 
     &     fluxa2(FACE3d2(ifirst,ilast,FLUXG),NEQU), 
     &     fluxb0(FACE3d0(ifirst,ilast,FLUXG),NEQU),
     &     fluxb1(FACE3d1(ifirst,ilast,FLUXG),NEQU), 
     &     fluxb2(FACE3d2(ifirst,ilast,FLUXG),NEQU), 
     &     tracelft0(FACE3d0(ifirst,ilast,FACEG),NEQU),
     &     tracergt0(FACE3d0(ifirst,ilast,FACEG),NEQU),
     &     tracelft1(FACE3d1(ifirst,ilast,FACEG),NEQU),
     &     tracergt1(FACE3d1(ifirst,ilast,FACEG),NEQU),
     &     tracelft2(FACE3d2(ifirst,ilast,FACEG),NEQU),
     &     tracergt2(FACE3d2(ifirst,ilast,FACEG),NEQU)
c
c***********************************************************************     
c
      integer ic0,ic1,ic2,k
      REAL trnsvers(NEQU)
c     REAL ttvlft(NEQU),ttvrgt(NEQU)     
      REAL ttv(NEQU)
      REAL v2norm,rho,vel1,vel2,vel0,gam_min_one 
c
      if (FLUXG.lt.1) then
         write(6,*) "flux ghosts < 1!"
         stop
      endif
     
c     ******************************************************************
c     * complete tracing at cell edges
c     ******************************************************************
      gam_min_one = gamma - one
c         
correc_flux3d(2,0,1,a0,a1,`ic1,ic2',`ic2,ic0')dnl
c         
correc_flux3d(1,2,0,a2,b0,`ic0,ic1',`ic1,ic2')dnl
c         
correc_flux3d(0,1,2,b1,b2,`ic2,ic0',`ic0,ic1')dnl
c   
      return
      end   
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
      subroutine fluxcalculation3d(dt,xcell0,xcell1,visco,dx,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  gamma,rpchoice,
     &  density,velocity,pressure,
     &  flux0,flux1,flux2,
     &  trlft0,trlft1,trlft2,
     &  trrgt0,trrgt1,trrgt2)
     
c***********************************************************************
      implicit none
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer xcell0,xcell1,visco,rpchoice
      REAL dt,dx(0:NDIM-1),gamma 
c variables in 2d cell indexed         
      REAL
     &     density(CELL3d(ifirst,ilast,CELLG)),
     &     velocity(CELL3d(ifirst,ilast,CELLG),NDIM),
     &     pressure(CELL3d(ifirst,ilast,CELLG))
c variables in 2d side indexed         
      REAL
     &     flux0(FACE3d0(ifirst,ilast,FLUXG),NEQU),
     &     flux1(FACE3d1(ifirst,ilast,FLUXG),NEQU), 
     &     flux2(FACE3d2(ifirst,ilast,FLUXG),NEQU), 
     &     trlft0(FACE3d0(ifirst,ilast,FACEG),NEQU),
     &     trrgt0(FACE3d0(ifirst,ilast,FACEG),NEQU),
     &     trlft1(FACE3d1(ifirst,ilast,FACEG),NEQU),
     &     trrgt1(FACE3d1(ifirst,ilast,FACEG),NEQU),
     &     trlft2(FACE3d2(ifirst,ilast,FACEG),NEQU),
     &     trrgt2(FACE3d2(ifirst,ilast,FACEG),NEQU)
c
c***********************************************************************     
c
      integer ic0,ic1,ic2,ie0,ie1,ie2,j
      REAL   stateL(NEQU),stateR(NEQU),
     &       riemst(NEQU)
      REAL mom0,mom1,mom2,Hent,v2norm,vel(0:NDIM-1)
      REAL mom0L,mom1L,mom2L,enerL,mom0R,mom1R,mom2R,enerR
      REAL maxeig, vcoef,vcorr(NEQU),gam_min_one
c
c variables for hllc scheme
      REAL aLsq,aRsq,keL,keR,flux(NEQU),diff(NEQU)
      REAL mfL,mfR,star(NEQU),sL,sM,sR
      REAL w,omw,hat(NEQU+1),denom
c
      if (FLUXG.lt.1) then
         write(6,*) "flux ghosts < 1!"
         stop
      endif
c
c***********************************************************************
c solve riemann problems for conservative flux
c  arguments: ( axis for RP, other axis, extra cells-direction)
c***********************************************************************
c      
      gam_min_one = gamma-one

c     write(6,*) "in fluxcalculation3d..."
c     write(6,*) "ifirst = ",ifirst0,ifirst1,ifirst2
c     write(6,*) "ilast = ",ilast0,ilast1,ilast2
c     write(6,*) "xcell0,xcell1 = ",xcell0,xcell1
c     write(6,*) "visco = ",visco
c     write(6,*) "gamma = ",gamma
c     write(6,*) "rpchoice = ",rpchoice
c     write(6,*) "gam_min_one = ",gam_min_one
c     call flush(6)
riemann_solve(0,1,2,`ic1,ic2',(xcell0+FLUXG-1),(xcell1+FLUXG-1))dnl

c
riemann_solve(1,0,2,`ic2,ic0',(xcell0+FLUXG-1),(xcell1+FLUXG-1))dnl

c
riemann_solve(2,0,1,`ic0,ic1',(xcell0+FLUXG-1),(xcell1+FLUXG-1))dnl

      if (visco.eq.1) then
c     write(6,*) "doing artificial viscosity"
c
artificial_viscosity1(0,1,2)dnl
c
artificial_viscosity1(1,2,0)dnl
c
artificial_viscosity1(2,0,1)dnl
c
      endif
      return
      end 
c***********************************************************************
c***********************************************************************
c***********************************************************************
  
      subroutine consdiff3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  dx,
     &  flux0,flux1,flux2,
     &  gamma,density,velocity,pressure)
c***********************************************************************
      implicit none
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
      integer ifirst0, ilast0,ifirst1, ilast1,ifirst2,ilast2
      REAL dx(0:NDIM-1)
      REAL
     &     flux0(FACE3d0(ifirst,ilast,FLUXG),NEQU),
     &     flux1(FACE3d1(ifirst,ilast,FLUXG),NEQU),
     &     flux2(FACE3d2(ifirst,ilast,FLUXG),NEQU),
     &     gamma,
     &     density(CELL3d(ifirst,ilast,CELLG)),
     &     velocity(CELL3d(ifirst,ilast,CELLG),NDIM),
     &     pressure(CELL3d(ifirst,ilast,CELLG))
c
      integer ic0,ic1,ic2,k
      REAL temp,v2norm,mom(NDIM),energy,
     &     gam_min_one
      
c***********************************************************************
c update conserved to full time
c note the permutation of indices in 2nd, 3rd coordinate directions
c***********************************************************************
      gam_min_one = gamma - one

      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
           do ic0=ifirst0,ilast0
             mom(1) = density(ic0,ic1,ic2)*velocity(ic0,ic1,ic2,1)
             mom(2) = density(ic0,ic1,ic2)*velocity(ic0,ic1,ic2,2)
             mom(3) = density(ic0,ic1,ic2)*velocity(ic0,ic1,ic2,3)
             v2norm = (velocity(ic0,ic1,ic2,1)**2+
     &                 velocity(ic0,ic1,ic2,2)**2+
     &                 velocity(ic0,ic1,ic2,3)**2)
             energy = pressure(ic0,ic1,ic2)/gam_min_one +
     &                        half*density(ic0,ic1,ic2)*v2norm
 
             density(ic0,ic1,ic2) = density(ic0,ic1,ic2)
     &          -(flux0(ic0+1,ic1,ic2,1)-flux0(ic0,ic1,ic2,1))/dx(0)
     &          -(flux1(ic1+1,ic2,ic0,1)-flux1(ic1,ic2,ic0,1))/dx(1)
     &          -(flux2(ic2+1,ic0,ic1,1)-flux2(ic2,ic0,ic1,1))/dx(2)
             density(ic0,ic1,ic2) = max(smallr,density(ic0,ic1,ic2))

             do k=1,3
               mom(k) = mom(k)
     &          -(flux0(ic0+1,ic1,ic2,k+1)-flux0(ic0,ic1,ic2,k+1))/dx(0)
     &          -(flux1(ic1+1,ic2,ic0,k+1)-flux1(ic1,ic2,ic0,k+1))/dx(1)
     &          -(flux2(ic2+1,ic0,ic1,k+1)-flux2(ic2,ic0,ic1,k+1))/dx(2)
               velocity(ic0,ic1,ic2,k) = mom(k)/density(ic0,ic1,ic2)
             enddo
             energy = energy
     &        -(flux0(ic0+1,ic1,ic2,NEQU)-flux0(ic0,ic1,ic2,NEQU))/dx(0)
     &        -(flux1(ic1+1,ic2,ic0,NEQU)-flux1(ic1,ic2,ic0,NEQU))/dx(1)
     &        -(flux2(ic2+1,ic0,ic1,NEQU)-flux2(ic2,ic0,ic1,NEQU))/dx(2)
c      
             v2norm = (velocity(ic0,ic1,ic2,1)**2+
     &            velocity(ic0,ic1,ic2,2)**2+velocity(ic0,ic1,ic2,3)**2)
             temp = energy - half*density(ic0,ic1,ic2)*v2norm
             pressure(ic0,ic1,ic2) = gam_min_one*temp
             pressure(ic0,ic1,ic2) = max(smallr,pressure(ic0,ic1,ic2))
           enddo
         enddo
      enddo
c
      return
      end
c***********************************************************************
include(FORTDIR/gas1d_approxrp3d.i)dnl
include(FORTDIR/gas1d_exactrp3d.i)dnl

c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to compute corner fluxes for 3d euler equations.
c
define(NDIM,3)dnl
define(NEQU,5)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl
include(FORTDIR/m4fluxcorner3d.i)dnl

      subroutine onethirdstate(dt,dx,idir,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  gamma,density,velocity,pressure,
     &  flux0,flux1,flux2,
     &  st3)
c***********************************************************************
      implicit none
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
c input arrays:
      integer idir,ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      REAL dt
c variables in 1d axis indexed
c
      REAL
     &     dx(0:NDIM-1), gamma
c variables in 3d cell indexed
      REAL
     &     density(CELL3d(ifirst,ilast,CELLG)),
     &     velocity(CELL3d(ifirst,ilast,CELLG),0:NDIM-1),
     &     pressure(CELL3d(ifirst,ilast,CELLG)),
     &     flux0(FACE3d0(ifirst,ilast,FLUXG),NEQU),
     &     flux1(FACE3d1(ifirst,ilast,FLUXG),NEQU),
     &     flux2(FACE3d2(ifirst,ilast,FLUXG),NEQU),
     &     st3(CELL3d(ifirst,ilast,CELLG),NEQU)
c
c***********************************************************************
c
      integer ic0,ic1,ic2,k
      REAL trnsvers(NEQU)
      REAL ttv(NEQU)
      REAL v2norm,rho,vel1,vel2,vel0,pres,
     &     gam_min_one
c
      if (FLUXG.lt.1) then
         write(6,*) "flux ghosts < 1!"
         stop
      endif

c     ******************************************************************
c     * complete tracing at cell edges
c     ******************************************************************
      gam_min_one = gamma - one
!$OMP PARALLEL DEFAULT(none)
!$OMPc SHARED(dt,dx,idir,gamma,gam_min_one,
!$OMPc        ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
!$OMPc        density,velocity,pressure,flux0,flux1,flux2,st3,
!$OMPc        CELLG,FLUXG,FACEG)
!$OMPc PRIVATE(ic0,ic1,ic2,k,
!$OMPc         trnsvers,ttv,rho,pres,vel0,vel1,vel2,v2norm)
c
      if (idir.eq.0) then
c
!$OMP DO SCHEDULE(DYNAMIC)
st_third(0,1,2,`ic1,ic2')dnl
!$OMP END DO
c
      elseif (idir.eq.1) then
c
!$OMP DO SCHEDULE(DYNAMIC)
st_third(1,2,0,`ic2,ic0')dnl
!$OMP END DO
c
      elseif (idir.eq.2) then
c
!$OMP DO SCHEDULE(DYNAMIC)
st_third(2,0,1,`ic0,ic1')dnl
!$OMP END DO
c
      endif
!$OMP END PARALLEL
      return
      end
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
      subroutine fluxthird(dt,dx,idir,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  gamma,rpchoice,
     &  density,velocity,pressure,
     &  st3,
     &  flux0,flux1,flux2 )

c***********************************************************************
      implicit none
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
c***********************************************************************
c input arrays:
      integer idir,ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      REAL dt,dx(0:NDIM-1),gamma
      integer rpchoice
      REAL
     &     density(CELL3d(ifirst,ilast,CELLG)),
     &     velocity(CELL3d(ifirst,ilast,CELLG),0:NDIM-1),
     &     pressure(CELL3d(ifirst,ilast,CELLG))
c variables in 2d side indexed
      REAL
     &     flux0(FACE3d0(ifirst,ilast,FLUXG),NEQU),
     &     flux1(FACE3d1(ifirst,ilast,FLUXG),NEQU),
     &     flux2(FACE3d2(ifirst,ilast,FLUXG),NEQU),
     &     st3(CELL3d(ifirst,ilast,CELLG),NEQU)
c
c***********************************************************************
c
      integer ic0,ic1,ic2
      REAL   stateL(4),stateR(4),
     &       riemst(4)
      REAL mom2,mom1,mom0,Hent,v2norm,vel(0:3-1),
     &     gam_min_one
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
      gam_min_one = gamma - one
!$OMP PARALLEL DEFAULT(none)
!$OMPc SHARED(dt,dx,idir,gamma,gam_min_one,rpchoice,
!$OMPc        ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
!$OMPc        density,velocity,pressure,flux0,flux1,flux2,st3,
!$OMPc        CELLG,FLUXG,FACEG,
!$OMPc        APPROX_RIEM_SOLVE,EXACT_RIEM_SOLVE,
!$OMPc        HLLC_RIEM_SOLVE)
!$OMPc PRIVATE(ic0,ic1,ic2,
!$OMPc         stateL,stateR,riemst,mom0,mom1,mom2,Hent,
!$OMPc         vel,v2norm)
c
      if (idir.eq.0) then
c
!$OMP DO SCHEDULE(DYNAMIC)
f_third(1,0,2,`ic2,ic0',`ic0,ic1-1,ic2')dnl
!$OMP END DO
c
!$OMP DO SCHEDULE(DYNAMIC)
f_third(2,0,1,`ic0,ic1',`ic0,ic1,ic2-1')dnl
!$OMP END DO
c
      elseif (idir.eq.1) then
c
!$OMP DO SCHEDULE(DYNAMIC)
f_third(0,1,2,`ic1,ic2',`ic0-1,ic1,ic2')dnl
!$OMP END DO
c
!$OMP DO SCHEDULE(DYNAMIC)
f_third(2,1,0,`ic0,ic1',`ic0,ic1,ic2-1')dnl
!$OMP END DO
c
      elseif (idir.eq.2) then
c
!$OMP DO SCHEDULE(DYNAMIC)
f_third(1,2,0,`ic2,ic0',`ic0,ic1-1,ic2')dnl
!$OMP END DO
c
!$OMP DO SCHEDULE(DYNAMIC)
f_third(0,2,1,`ic1,ic2',`ic0-1,ic1,ic2')dnl
!$OMP END DO
c
      endif
c
!$OMP END PARALLEL
      return
      end
c***********************************************************************
c***********************************************************************
c***********************************************************************
      subroutine fluxcorrecjt(dt,dx,idir,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  gamma,density,velocity,pressure,
     &  flux0,flux1,flux2,
     &  tracelft0,tracelft1,tracelft2,
     &  tracergt0,tracergt1,tracergt2)
c***********************************************************************
      implicit none
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
c***********************************************************************
c input arrays:
      integer idir,ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      REAL dt
c variables in 1d axis indexed
c
      REAL
     &     dx(0:3-1)
c variables in 2d cell indexed
      REAL
     &     density(CELL3d(ifirst,ilast,CELLG)),
     &     velocity(CELL3d(ifirst,ilast,CELLG),0:NDIM-1),
     &     pressure(CELL3d(ifirst,ilast,CELLG)),
     &     flux0(FACE3d0(ifirst,ilast,FLUXG),NEQU),
     &     flux1(FACE3d1(ifirst,ilast,FLUXG),NEQU),
     &     flux2(FACE3d2(ifirst,ilast,FLUXG),NEQU),
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
      REAL ttv(NEQU)
      REAL v2norm,rho,vel1,vel2,vel0,gamma,
     &     gam_min_one
c
      if (FLUXG.lt.1) then
         write(6,*) "flux ghosts < 1!"
         stop
      endif

c     ******************************************************************
c     * complete tracing at cell edges
c     ******************************************************************
      gam_min_one = gamma - one
!$OMP PARALLEL DEFAULT(none)
!$OMPc SHARED(dt,dx,idir,gamma,gam_min_one,
!$OMPc        ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
!$OMPc        density,velocity,pressure,flux0,flux1,flux2,
!$OMPc        CELLG,FLUXG,FACEG,
!$OMPc        tracelft0,tracelft1,tracelft2,
!$OMPc        tracergt0,tracergt1,tracergt2)
!$OMPc PRIVATE(ic0,ic1,ic2,k,
!$OMPc         vel0,vel1,vel2,rho,
!$OMPc         trnsvers,ttv,v2norm)
c
      if (idir.eq.0) then
c
correc_fluxjt(0,2,1,`ic2,ic0',`ic0,ic1')dnl
c
      elseif (idir.eq.1) then
c
correc_fluxjt(1,0,2,`ic0,ic1',`ic1,ic2')dnl
c
      elseif (idir.eq.2) then
c
correc_fluxjt(2,1,0,`ic1,ic2',`ic2,ic0')dnl
c
      endif
c
!$OMP END PARALLEL
      return
      end
c
c***********************************************************************

c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for corner fluxes in 2d.
c
define(NDIM,3)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl
include(FORTDIR/m4fluxcorner3d.i)dnl

      subroutine onethirdstate3d(dt,dx,idir,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  advecspeed,uval,
     &  flux0,flux1,flux2,
     &  st3)
c***********************************************************************
      implicit none
include(FORTDIR/const.i)dnl
include(FORTDIR/probparams.i)dnl
c***********************************************************************     
c input arrays:
      integer idir,ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      REAL dt 
c variables in 1d axis indexed
c
      REAL dx(0:NDIM-1)
c variables in 2d cell indexed         
      REAL
     &     advecspeed(0:NDIM-1),
     &     uval(CELL3d(ifirst,ilast,CELLG)),
     &     flux0(FACE3d0(ifirst,ilast,FLUXG)),
     &     flux1(FACE3d1(ifirst,ilast,FLUXG)),
     &     flux2(FACE3d2(ifirst,ilast,FLUXG)),
     &     st3(CELL3d(ifirst,ilast,CELLG))
c
c***********************************************************************     
c
      integer ic0,ic1,ic2
      REAL trnsvers
     
c     ******************************************************************
c     * complete tracing at cell edges
c     ******************************************************************
c         
      if (idir.eq.0) then
c
st_third(0,1,2,`ic1,ic2')dnl
c
      elseif (idir.eq.1) then
c
st_third(1,2,0,`ic2,ic0')dnl
c
      elseif (idir.eq.2) then
c
st_third(2,0,1,`ic0,ic1')dnl
c
      endif
      return
      end   
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
      subroutine fluxthird3d(dt,dx,idir,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  advecspeed,
     &  uval,
     &  st3,
     &  flux0,flux1,flux2)
     
c***********************************************************************
      implicit none
include(FORTDIR/const.i)dnl
include(FORTDIR/probparams.i)dnl
c***********************************************************************
c***********************************************************************
c input arrays:
      integer idir,ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      REAL dt 
      REAL 
     &     dx(0:NDIM-1)
c variables in 2d cell indexed         
      REAL
     &     advecspeed(0:NDIM-1),
     &     uval(CELL3d(ifirst,ilast,CELLG))
c variables in 2d side indexed         
      REAL
     &     flux0(FACE3d0(ifirst,ilast,FLUXG)),
     &     flux1(FACE3d1(ifirst,ilast,FLUXG)),
     &     flux2(FACE3d2(ifirst,ilast,FLUXG)),
     &     st3(CELL3d(ifirst,ilast,CELLG)) 
c
c***********************************************************************     
c
      integer ic0,ic1,ic2
      REAL   riemst
c
c***********************************************************************
c solve riemann problems for conservative flux
c  arguments: ( axis for RP, other axis, extra cells-direction)
c***********************************************************************
c      
      if (idir.eq.0) then
c
f_third(1,0,2,`ic2,ic0',`ic0,ic1-1,ic2')dnl
c
f_third(2,0,1,`ic0,ic1',`ic0,ic1,ic2-1')dnl
c
      elseif (idir.eq.1) then
c
f_third(0,1,2,`ic1,ic2',`ic0-1,ic1,ic2')dnl
c
f_third(2,1,0,`ic0,ic1',`ic0,ic1,ic2-1')dnl
c
      elseif (idir.eq.2) then
c
f_third(1,2,0,`ic2,ic0',`ic0,ic1-1,ic2')dnl
c
f_third(0,2,1,`ic1,ic2',`ic0-1,ic1,ic2')dnl
c
      endif
c
c      call flush(6)     
      return
      end 
c***********************************************************************
c***********************************************************************
c***********************************************************************
      subroutine fluxcorrecjt3d(dt,dx,idir,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  advecspeed,uval,
     &  flux0,flux1,flux2,
     &  tracelft0,tracelft1,tracelft2,
     &  tracergt0,tracergt1,tracergt2)
c***********************************************************************
      implicit none
include(FORTDIR/const.i)dnl
include(FORTDIR/probparams.i)dnl
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
     &     advecspeed(0:NDIM-1),
     &     uval(CELL3d(ifirst,ilast,CELLG)),
     &     flux0(FACE3d0(ifirst,ilast,FLUXG)),
     &     flux1(FACE3d1(ifirst,ilast,FLUXG)),
     &     flux2(FACE3d2(ifirst,ilast,FLUXG)),
     &     tracelft0(FACE3d0(ifirst,ilast,FACEG)),
     &     tracergt0(FACE3d0(ifirst,ilast,FACEG)),
     &     tracelft1(FACE3d1(ifirst,ilast,FACEG)),
     &     tracergt1(FACE3d1(ifirst,ilast,FACEG)),
     &     tracelft2(FACE3d2(ifirst,ilast,FACEG)),
     &     tracergt2(FACE3d2(ifirst,ilast,FACEG))
c
c***********************************************************************     
c
      integer ic0,ic1,ic2
      REAL trnsvers
c     REAL ttvlft,ttvrgt     
     
c     ******************************************************************
c     * complete tracing at cell edges
c     ******************************************************************
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
c      call flush(6)     
      return
      end   
c
c***********************************************************************
c***********************************************************************
c***********************************************************************

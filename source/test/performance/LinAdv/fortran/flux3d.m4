c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines to evaluate fluxes in 3d.
c
define(NDIM,3)dnl
define(NEQU,1)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl
include(FORTDIR/m4flux3d.i)dnl

      subroutine fluxcorrec3d2d(dt,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  dx,advecspeed,idir,
     &  uval,
     &  flux0,flux1,flux2,
     &  tracelft0,tracelft1,tracelft2,
     &  tracergt0,tracergt1,tracergt2,
     &  ttracelft0,ttracelft1,ttracelft2,
     &  ttracergt0,ttracergt1,ttracergt2)
c***********************************************************************
      implicit none
include(FORTDIR/const.i)dnl
include(FORTDIR/probparams.i)dnl
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
     &     advecspeed(0:NDIM-1),
     &     uval(CELL3d(ifirst,ilast,CELLG)),
c
     &     flux0(FACE3d0(ifirst,ilast,FLUXG)),
     &     flux1(FACE3d1(ifirst,ilast,FLUXG)), 
     &     flux2(FACE3d2(ifirst,ilast,FLUXG)), 
c
     &     tracelft0(FACE3d0(ifirst,ilast,FACEG)),
     &     tracelft1(FACE3d1(ifirst,ilast,FACEG)),
     &     tracelft2(FACE3d2(ifirst,ilast,FACEG)),
     &     tracergt0(FACE3d0(ifirst,ilast,FACEG)),
     &     tracergt1(FACE3d1(ifirst,ilast,FACEG)),
     &     tracergt2(FACE3d2(ifirst,ilast,FACEG)),
c
     &     ttracelft0(FACE3d0(ifirst,ilast,FACEG)),
     &     ttracelft1(FACE3d1(ifirst,ilast,FACEG)),
     &     ttracelft2(FACE3d2(ifirst,ilast,FACEG)),
     &     ttracergt0(FACE3d0(ifirst,ilast,FACEG)),
     &     ttracergt1(FACE3d1(ifirst,ilast,FACEG)),
     &     ttracergt2(FACE3d2(ifirst,ilast,FACEG)) 
c
c***********************************************************************     
c
      integer ic0,ic1,ic2
      REAL trnsvers
     
c     ******************************************************************
c     * complete tracing at cell edges
c     ******************************************************************
c         
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
      subroutine fluxcorrec3d3d(dt,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  dx,
     &  advecspeed,uval,
     &  fluxa0,fluxa1,fluxa2,
     &  fluxb0,fluxb1,fluxb2,
     &  tracelft0,tracelft1,tracelft2,
     &  tracergt0,tracergt1,tracergt2)
c***********************************************************************
      implicit none
include(FORTDIR/const.i)dnl
include(FORTDIR/probparams.i)dnl
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
     &     advecspeed(0:NDIM-1),
     &     uval(CELL3d(ifirst,ilast,CELLG)),
     &     fluxa0(FACE3d0(ifirst,ilast,FLUXG)),
     &     fluxa1(FACE3d1(ifirst,ilast,FLUXG)), 
     &     fluxa2(FACE3d2(ifirst,ilast,FLUXG)), 
     &     fluxb0(FACE3d0(ifirst,ilast,FLUXG)),
     &     fluxb1(FACE3d1(ifirst,ilast,FLUXG)), 
     &     fluxb2(FACE3d2(ifirst,ilast,FLUXG)), 
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
     
c     ******************************************************************
c     * complete tracing at cell edges
c     ******************************************************************
c         
correc_flux3d(2,0,1,a0,a1,`ic1,ic2',`ic2,ic0')dnl
c         
correc_flux3d(1,2,0,a2,b0,`ic0,ic1',`ic1,ic2')dnl
c         
correc_flux3d(0,1,2,b1,b2,`ic2,ic0',`ic0,ic1')dnl
c   
c      call flush(6)     
      return
      end   
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
      subroutine fluxcalculation3d(dt,xcell0,xcell1,visco,dx,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  advecspeed,
     &  uval,
     &  flux0,flux1,flux2,
     &  trlft0,trlft1,trlft2,
     &  trrgt0,trrgt1,trrgt2)
     
c***********************************************************************
      implicit none
include(FORTDIR/const.i)dnl
include(FORTDIR/probparams.i)dnl
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer xcell0,xcell1,visco
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
     &     trlft0(FACE3d0(ifirst,ilast,FACEG)),
     &     trrgt0(FACE3d0(ifirst,ilast,FACEG)),
     &     trlft1(FACE3d1(ifirst,ilast,FACEG)),
     &     trrgt1(FACE3d1(ifirst,ilast,FACEG)),
     &     trlft2(FACE3d2(ifirst,ilast,FACEG)),
     &     trrgt2(FACE3d2(ifirst,ilast,FACEG))
c
c***********************************************************************     
c
      integer ic0,ic1,ic2,ie0,ie1,ie2
      REAL   riemst
c
c***********************************************************************
c solve riemann problems for conservative flux
c  arguments: ( axis for RP, other axis, extra cells-direction)
c***********************************************************************
c      

riemann_solve(0,1,2,`ic1,ic2',xcell0,xcell1)dnl

c
riemann_solve(1,0,2,`ic2,ic0',xcell0,xcell1)dnl

c
riemann_solve(2,0,1,`ic0,ic1',xcell0,xcell1)dnl

      if (visco.eq.1) then
      write(6,*) "doing artificial viscosity"
c
crtificial_viscosity1(0,1,2)dnl
c
crtificial_viscosity1(1,2,0)dnl
c
crtificial_viscosity1(2,0,1)dnl
c
      endif
c      call flush(6)     
      return
      end 
c***********************************************************************
c***********************************************************************
c***********************************************************************
  
      subroutine consdiff3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  dx,
     &  flux0,flux1,flux2,
     &  advecspeed,uval)
c***********************************************************************
      implicit none
include(FORTDIR/const.i)dnl
include(FORTDIR/probparams.i)dnl
c***********************************************************************
      integer ifirst0, ilast0,ifirst1, ilast1,ifirst2,ilast2
      REAL dx(0:NDIM-1)
      REAL
     &     flux0(FACE3d0(ifirst,ilast,FLUXG)),
     &     flux1(FACE3d1(ifirst,ilast,FLUXG)),
     &     flux2(FACE3d2(ifirst,ilast,FLUXG)),
     &     advecspeed(0:NDIM-1),
     &     uval(CELL3d(ifirst,ilast,CELLG))
c
      integer ic0,ic1,ic2
      
c***********************************************************************
c update velocity to full time
c note the reversal of indices in 2nd coordinate direction
c***********************************************************************
c***********************************************************************
      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
           do ic0=ifirst0,ilast0
             uval(ic0,ic1,ic2) = uval(ic0,ic1,ic2)
     &          -(flux0(ic0+1,ic1,ic2)-flux0(ic0,ic1,ic2))/dx(0)
     &          -(flux1(ic1+1,ic2,ic0)-flux1(ic1,ic2,ic0))/dx(1)
     &          -(flux2(ic2+1,ic0,ic1)-flux2(ic2,ic0,ic1))/dx(2)
           enddo
         enddo
      enddo
      return
      end
c***********************************************************************

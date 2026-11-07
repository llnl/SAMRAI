c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for flux computation in 2d.
c
define(NDIM,2)dnl
define(NEQU,1)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl
include(FORTDIR/m4flux2d.i)dnl

      subroutine fluxcorrec(dt,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  dx,
     &  advecspeed,
     &  flux0,flux1,
     &  trlft0,trlft1,
     &  trrgt0,trrgt1)
c***********************************************************************
      implicit none
include(FORTDIR/const.i)dnl
include(FORTDIR/probparams.i)dnl
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
     &     advecspeed(0:NDIM-1),
     &     flux0(FACE2d0(ifirst,ilast,FLUXG)),
     &     flux1(FACE2d1(ifirst,ilast,FLUXG)), 
     &     trlft0(FACE2d0(ifirst,ilast,FACEG)),
     &     trrgt0(FACE2d0(ifirst,ilast,FACEG)),
     &     trlft1(FACE2d1(ifirst,ilast,FACEG)),
     &     trrgt1(FACE2d1(ifirst,ilast,FACEG))
c
c***********************************************************************     
c
      integer ic0,ic1
      REAL trnsvers
     
c     write(6,*) "In fluxcorrec()"
c     ******************************************************************
c     * complete tracing at cell edges
c     ******************************************************************
c   correct the 1-direction with 0-fluxes
c     write(6,*) " correct the 1-direction with 0-fluxes"
      do ic1=ifirst1-1,ilast1+1
        do ic0=ifirst0-1,ilast0+1
          trnsvers= (flux0(ic0+1,ic1)-flux0(ic0,ic1))*0.5/dx(0)

          trrgt1(ic1  ,ic0)= trrgt1(ic1  ,ic0) - trnsvers
          trlft1(ic1+1,ic0)= trlft1(ic1+1,ic0) - trnsvers
        enddo
      enddo
c     call flush(6)     

c   correct the 0-direction with 1-fluxes
c     write(6,*) " correct the 0-direction with 1-fluxes"
      do ic0=ifirst0-1,ilast0+1
        do ic1=ifirst1-1,ilast1+1
          trnsvers= (flux1(ic1+1,ic0)-flux1(ic1,ic0))*0.5/dx(1)
          trrgt0(ic0  ,ic1)= trrgt0(ic0  ,ic1) - trnsvers
          trlft0(ic0+1,ic1)= trlft0(ic0+1,ic1) - trnsvers
        enddo
      enddo
c     call flush(6)     
      return
      end   
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
      subroutine fluxcalculation2d(dt,extra_cell,visco,dx,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  advecspeed,
     &  flux0,flux1,
     &  trlft0,trlft1,trrgt0,trrgt1)
c***********************************************************************
      implicit none
include(FORTDIR/const.i)dnl
include(FORTDIR/probparams.i)dnl
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,extra_cell,visco
      REAL dt 
      REAL 
     &     dx(0:NDIM-1)
c variables in 2d cell indexed         
      REAL
     &     advecspeed(0:NDIM-1),
c variables in 2d side indexed         
     &     flux0(FACE2d0(ifirst,ilast,FLUXG)),
     &     flux1(FACE2d1(ifirst,ilast,FLUXG)), 
     &     trlft0(FACE2d0(ifirst,ilast,FACEG)),
     &     trrgt0(FACE2d0(ifirst,ilast,FACEG)),
     &     trlft1(FACE2d1(ifirst,ilast,FACEG)),
     &     trrgt1(FACE2d1(ifirst,ilast,FACEG))
c
c***********************************************************************     
c
      integer ic0,ic1,ie0,ie1
      REAL   riemst
c
c***********************************************************************
c solve riemann problems for conservative flux
c  arguments: ( axis for RP, other axis, extra cells-direction)
c***********************************************************************
c      
c     write(6,*) "In fluxcalculation2d(",extra_cell,")"
c     write(6,*) "ifirst0,ilast0,ifirst1,ilast1,extra_cell",
c    &       ifirst0,ilast0,ifirst1,ilast1,extra_cell

riemann_solve(0,1,extra_cell)dnl
c
riemann_solve(1,0,extra_cell)dnl

      if (visco.eq.1) then
      write(6,*) "doing artificial viscosity"
c
crtificial_viscosity1(0,1)dnl
c
crtificial_viscosity1(1,0)dnl
c
      endif
c     call flush(6)     
      return
      end 
c***********************************************************************
c***********************************************************************
c***********************************************************************
  
      subroutine consdiff2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  dx,
     &  flux0,flux1,
     &  advecspeed,src,uval)
c***********************************************************************
      implicit none
include(FORTDIR/const.i)dnl
include(FORTDIR/probparams.i)dnl
c***********************************************************************
      integer ifirst0, ilast0,ifirst1, ilast1
      REAL dx(0:NDIM-1)
      REAL
     &     flux0(FACE2d0(ifirst,ilast,FLUXG)),
     &     flux1(FACE2d1(ifirst,ilast,FLUXG)),
     &     advecspeed(0:NDIM-1),src,
     &     uval(CELL2d(ifirst,ilast,CELLG))
c
      integer ic0,ic1
      
c***********************************************************************
c update velocity to full time
c note the reversal of indices in 2nd coordinate direction
c***********************************************************************
c***********************************************************************
c     write(6,*) "at top of consdiff2d"
c         call flush(6)
c     write(6,*) "flux0"
c     do  ic1=ifirst1,ilast1
c       do  ic0=ifirst0,ilast0+1
c         write(6,*) "ic0,flux0= ",ic0,ic1,
c    &                  flux0(ic0,ic1,1),flux0(ic0,ic1,2),
c    &                  flux0(ic0,ic1,3),flux0(ic0,ic1,4)
c         call flush(6)
c       enddo
c     enddo
c     write(6,*) "flux1"
c     do  ic1=ifirst1,ilast1
c       do  ic0=ifirst0,ilast0+1
c         write(6,*) "ic0,flux1= ",ic0,ic1,
c    &                  flux1(ic0,ic1,1),flux1(ic0,ic1,2),
c    &                  flux1(ic0,ic1,3),flux1(ic0,ic1,4)
c         call flush(6)
c       enddo
c     enddo
c     do  ic1=ifirst1,ilast1
c       do  ic0=ifirst0,ilast0
c         write(6,*) "ic0,ic1,all = ",ic0,ic1,
c    &                        density(ic0,ic1),
c    &                        velocity(ic0,ic1,1),velocity(ic0,ic1,2),
c    &                        pressure(ic0,ic1)
c         call flush(6)
c       enddo
c     enddo
c***********************************************************************
      do ic1=ifirst1,ilast1
        do ic0=ifirst0,ilast0
          uval(ic0,ic1) = uval(ic0,ic1) +src
     &          -(flux0(ic0+1,ic1)-flux0(ic0,ic1))/dx(0)
     &          -(flux1(ic1+1,ic0)-flux1(ic1,ic0))/dx(1)
        enddo
      enddo
c***********************************************************************
c     write(6,*) "in consdiff2d"
c     do  ic1=ifirst1,ilast1+1
c       do  ic0=ifirst0,ilast0+1
c         write(6,*) "ic0,ic1,flux0= ",ic0,ic1,
c    &                        flux0(ic0,ic1,1),flux0(ic0,ic1,2),
c    &                        flux0(ic0,ic1,3),flux0(ic0,ic1,4) 
c         call flush(6)
c       enddo
c     enddo
c     do  ic1=ifirst1,ilast1+1
c       do  ic0=ifirst0,ilast0+1
c         write(6,*) "ic0,ic1,flux1= ",ic0,ic1,
c    &                        flux1(ic0,ic1,1),flux1(ic0,ic1,2),
c    &                        flux1(ic0,ic1,3),flux1(ic0,ic1,4) 
c         call flush(6)
c       enddo
c     enddo
c***********************************************************************
c***********************************************************************
c     do  ic1=ifirst1,ilast1
c       do  ic0=ifirst0,ilast0
c         write(6,*) "ic0,ic1,all = ",ic0,ic1,
c    &                        density(ic0,ic1),
c    &                        velocity(ic0,ic1,1),velocity(ic0,ic1,2),
c    &                        momentum(ic0,ic1,1),momentum(ic0,ic1,2),
c    &                        pressure(ic0,ic1),
c    &                        energy(ic0,ic1)
c         call flush(6)
c       enddo
c     enddo
      return
      end

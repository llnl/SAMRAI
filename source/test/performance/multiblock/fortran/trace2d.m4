c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine for trace in 2d.
c
define(NDIM,2)dnl
define(NEQU,4)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl
include(FORTDIR/m4trace2d.i)dnl

      subroutine inittraceflux2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  uval,
     &  tracelft0,tracelft1,
     &  tracergt0,tracergt1,
     &  fluxriem0,fluxriem1)
c***********************************************************************
      implicit none
include(FORTDIR/const.i)dnl
include(FORTDIR/probparams.i)dnl
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      REAL
     &     uval(CELL2d(ifirst,ilast,CELLG)),
     &     fluxriem0(FACE2d0(ifirst,ilast,FLUXG)),
     &     tracelft0(FACE2d0(ifirst,ilast,FACEG)),
     &     tracergt0(FACE2d0(ifirst,ilast,FACEG)),
     &     fluxriem1(FACE2d1(ifirst,ilast,FLUXG)),
     &     tracelft1(FACE2d1(ifirst,ilast,FACEG)),
     &     tracergt1(FACE2d1(ifirst,ilast,FACEG))
c***********************************************************************
c initialize left and right states at cell edges
c (first-order upwind)
c***********************************************************************
c
      integer ic0,ic1,ie0,ie1
c*********************************************************************** 

trace_init(0,1,2,`ie0-1,ic1',`ie0,ic1')dnl
 
trace_init(1,0,2,`ic0,ie1-1',`ic0,ie1')dnl
 
c
c     we initialize the flux to be zero

      do ic1=ifirst1-FLUXG,ilast1+FLUXG
        do ie0=ifirst0-FLUXG,ilast0+FLUXG+1
            fluxriem0(ie0,ic1) = zero
        enddo
      enddo
c
      do ic0=ifirst0-FLUXG,ilast0+FLUXG
        do ie1=ifirst1-FLUXG,ilast1+FLUXG+1
            fluxriem1(ie1,ic0) = zero
        enddo
      enddo
c
c     call flush(6)
      return
      end 
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
      subroutine chartracing2d0(dt,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  mc,
     &  dx,
     &  advecspeed,igdnv,
     &  uval,
     &  tracelft,tracergt,
     &  ttcelslp,ttedgslp, 
     &  ttraclft,ttracrgt)
c***********************************************************************
      implicit none
include(FORTDIR/probparams.i)dnl
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer mc,igdnv
      REAL dt 
c variables in 1d axis indexed
      REAL 
     &     dx,advecspeed
c variables in 2d axis indexed         
      REAL
     &     uval(CELL2d(ifirst,ilast,CELLG)),
     &     tracelft(FACE2d0(ifirst,ilast,FACEG)),
     &     tracergt(FACE2d0(ifirst,ilast,FACEG))
c  side variables ifirst0 to ifirst0+mc plus ghost cells
      REAL
     &  ttedgslp(ifirst0-FACEG:ifirst0+mc+FACEG),
     &  ttraclft(ifirst0-FACEG:ifirst0+mc+FACEG),
     &  ttracrgt(ifirst0-FACEG:ifirst0+mc+FACEG)
c  cell variables ifirst0 to ifirst0+mc-1 plus ghost cells
      REAL 
     &  ttcelslp(ifirst0-CELLG:ifirst0+mc-1+CELLG)
c***********************************************************************
c
      integer ic0,ic1,idir
c*********************************************************************** 
c***********************************************************************
c trace higher order states at cell edges
c***********************************************************************
      idir = 0
trace_call(0,1)dnl
c***********************************************************************
      return
      end
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
      subroutine chartracing2d1(dt,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  mc,
     &  dx,
     &  advecspeed,igdnv,
     &  uval,
     &  tracelft,tracergt,
     &  ttcelslp,ttedgslp,  
     &  ttraclft,ttracrgt)
c***********************************************************************
      implicit none
include(FORTDIR/probparams.i)dnl
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer mc,igdnv
      REAL dt
c variables in 1d axis indexed
      REAL 
     &     dx,advecspeed
c variables in 2d axis indexed         
      REAL
     &     uval(CELL2d(ifirst,ilast,CELLG)),
     &     tracelft(FACE2d1(ifirst,ilast,FACEG)),
     &     tracergt(FACE2d1(ifirst,ilast,FACEG))
c  side variables ifirst1 to ifirst1+mc plus ghost cells
      REAL
     &  ttedgslp(ifirst1-FACEG:ifirst1+mc+FACEG),
     &  ttraclft(ifirst1-FACEG:ifirst1+mc+FACEG),
     &  ttracrgt(ifirst1-FACEG:ifirst1+mc+FACEG)
c  cell variables ifirst1 to ifirst1+mc-1 plus ghost cells
      REAL 
     &  ttcelslp(ifirst1-CELLG:ifirst1+mc-1+CELLG)
c***********************************************************************
c
      integer ic0,ic1,idir
c*********************************************************************** 
c***********************************************************************
c trace higher order states at cell edges
c***********************************************************************
      idir = 1
trace_call(1,0)dnl
c***********************************************************************
      return
      end
c
c***********************************************************************
c***********************************************************************
c***********************************************************************

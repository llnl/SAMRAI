c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for trace for 2d euler equations.
c
define(NDIM,2)dnl
define(NEQU,4)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl
include(FORTDIR/m4trace2d.i)dnl

      subroutine computesound2d(ifirst0,ilast0,ifirst1,ilast1,
     &  gamma,density, velocity,pressure,sound)
c***********************************************************************
      implicit none
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
c variables indexed as 2dimensional
      REAL
     &     gamma,
     &     density(CELL2d(ifirst,ilast,CELLG)),
     &     velocity(CELL2d(ifirst,ilast,CELLG),0:NDIM-1),
     &     pressure(CELL2d(ifirst,ilast,CELLG)),
     &     sound(CELL2d(ifirst,ilast,CELLG))
c
      integer ic0,ic1
c
c***********************************************************************
c
c     write(6,*) "in computesound"
c     do  ic1=ifirst1-4,ilast1+4
c       do  ic0=ifirst0-4,ilast0+4
c         write(6,*) "ic01,d_p_v=",ic0,ic1,
c    &         density(ic0,ic1),pressure(ic0,ic1),
c    &         velocity(ic0,ic1,0),velocity(ic0,ic1,1) 
c         call flush(6)
c       enddo
c     enddo
c     write(6,*)

      do  ic0=ifirst0-CELLG,ilast0+CELLG
        do  ic1=ifirst1-CELLG,ilast1+CELLG
c         write(6,*) "        density,pressure = ",ic0,ic1,
c    &                        density(ic0,ic1),pressure(ic0,ic1)
c         write(6,*) "        velocity,        = ",ic0,ic1,
c    &                        velocity(ic0,ic1,0),velocity(ic0,ic1,1)
c         call flush(6)
            sound(ic0,ic1) = 
     &        sqrt(max(smallr,gamma*pressure(ic0,ic1)/density(ic0,ic1)))
c           call flush(6)
c
        enddo
      enddo

c     write(6,*) "leaving computesound"
c      call flush(6)
      return
      end 
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
c
      subroutine inittraceflux2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  density,velocity,pressure,
     &  tracelft0,tracelft1,
     &  tracergt0,tracergt1,
     &  fluxriem0,fluxriem1)
c***********************************************************************
      implicit none
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      REAL
     &       density(CELL2d(ifirst,ilast,CELLG)),
     &      velocity(CELL2d(ifirst,ilast,CELLG),NDIM),
     &      pressure(CELL2d(ifirst,ilast,CELLG)),
     &     fluxriem0(FACE2d0(ifirst,ilast,FLUXG),NEQU),
     &     tracelft0(FACE2d0(ifirst,ilast,FACEG),NEQU),
     &     tracergt0(FACE2d0(ifirst,ilast,FACEG),NEQU),
     &     fluxriem1(FACE2d1(ifirst,ilast,FLUXG),NEQU),
     &     tracelft1(FACE2d1(ifirst,ilast,FACEG),NEQU),
     &     tracergt1(FACE2d1(ifirst,ilast,FACEG),NEQU)
      integer ic0,ic1,k,ie0,ie1
c***********************************************************************
c initialize left and right states at cell edges
c (first-order upwind)
c***********************************************************************
      if (FLUXG.lt.1) then
         write(6,*) "flux ghosts < 1!" 
         stop
      endif
c*********************************************************************** 
c     write(6,*) "in init_trace_flux"
c     call flush(6)
c     write(6,*) "ifirst0,ilast0,ifirst1,ilast1",
c    &     ifirst0,ilast0,ifirst1,ilast1
c     write(6,*) "cell ghosts ", CELLG
c     write(6,*) "face ghosts ", FACEG
c     call flush(6)

c     write(6,*) " "
c     write(6,*) " In trace_init0"
trace_init(0,1,2,`ie0-1,ic1',`ie0,ic1')dnl
 
c     write(6,*) " "
c     write(6,*) " In trace_init1"
trace_init(1,0,2,`ic0,ie1-1',`ic0,ie1')dnl
c     write(6,*) " "
 
c
c     we initialize the flux to be 0

      do ic1=ifirst1-FLUXG,ilast1+FLUXG
        do ie0=ifirst0-FLUXG,ilast0+FLUXG+1
          do k=1,NEQU  
            fluxriem0(ie0,ic1,k) = zero
          enddo  
        enddo
      enddo
c
      do ic0=ifirst0-FLUXG,ilast0+FLUXG
        do ie1=ifirst1-FLUXG,ilast1+FLUXG+1
          do k=1,NEQU  
            fluxriem1(ie1,ic0,k) = zero
          enddo  
        enddo
      enddo
c
c      call flush(6)
      return
      end 
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
      subroutine chartracing2d0(dt,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  mc, dx,
     &  gamma,igdnv,
     &  sound,
     &  tracelft,tracergt,
     &  ttcelslp, ttedgslp, 
     &  ttsound,
     &  ttraclft, ttracrgt)
c***********************************************************************
      implicit none
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer mc,igdnv 
      REAL dt 
c variables in 1d axis indexed
      REAL 
     &     dx, gamma
      REAL 
     &     sound(CELL2d(ifirst,ilast,CELLG))
      REAL 
     &     tracelft(FACE2d0(ifirst,ilast,FACEG),NEQU),
     &     tracergt(FACE2d0(ifirst,ilast,FACEG),NEQU)
c  side variables ifirst0 to ifirst0+mc plus ghost cells
      REAL
     &  ttedgslp(ifirst0-FACEG:ifirst0+mc+FACEG,NEQU),
     &  ttraclft(ifirst0-FACEG:ifirst0+mc+FACEG,NEQU),
     &  ttracrgt(ifirst0-FACEG:ifirst0+mc+FACEG,NEQU)
c  cell variables ifirst0 to ifirst0+mc-1 plus ghost cells
      REAL 
     &  ttcelslp(ifirst0-CELLG:ifirst0+mc-1+CELLG,NEQU),
     &  ttsound(ifirst0-CELLG:ifirst0+mc-1+CELLG)
c***********************************************************************
c
      integer ic0,ic1,k,idir
c*********************************************************************** 
c***********************************************************************
c trace higher order states at cell edges
c***********************************************************************
c     write(6,*) "traced right chartracing2d0"
c     do  ic1=ifirst1-FACEG,ilast1+FACEG
c       do  ic0=ifirst0-FACEG,ilast0+FACEG+1
c         write(6,*) "ic,state=",ic0,ic1,
c    &         tracergt(ic0,ic1,1),tracergt(ic0,ic1,2),
c    &         tracergt(ic0,ic1,3),tracergt(ic0,ic1,4) 
c         call flush(6)
c       enddo
c     enddo
c     write(6,*)
c
c     call flush(6)
c
      idir = 0
trace_call(0,1)dnl
c
c***********************************************************************
c     write(6,*) "leaving chartracing2d"
c      call flush(6)
      return
      end
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
      subroutine chartracing2d1(dt,
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  mc, dx,
     &  gamma,igdnv,
     &  sound,
     &  tracelft,tracergt,
     &  ttcelslp, ttedgslp,  
     &  ttsound,
     &  ttraclft, ttracrgt)
c***********************************************************************
      implicit none
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer mc,igdnv
      REAL dt 
c variables in 1d axis indexed
      REAL 
     &     dx, gamma
      REAL 
     &     sound(CELL2d(ifirst,ilast,CELLG))
      REAL 
     &     tracelft(FACE2d1(ifirst,ilast,FACEG),NEQU),
     &     tracergt(FACE2d1(ifirst,ilast,FACEG),NEQU)
c  side variables ifirst1 to ifirst1+mc plus ghost cells
      REAL
     &  ttedgslp(ifirst1-FACEG:ifirst1+mc+FACEG,NEQU),
     &  ttraclft(ifirst1-FACEG:ifirst1+mc+FACEG,NEQU),
     &  ttracrgt(ifirst1-FACEG:ifirst1+mc+FACEG,NEQU)
c  cell variables ifirst1 to ifirst1+mc-1 plus ghost cells
      REAL 
     &  ttcelslp(ifirst1-CELLG:ifirst1+mc-1+CELLG,NEQU),
     &  ttsound(ifirst1-CELLG:ifirst1+mc-1+CELLG)
c***********************************************************************
c
      integer ic0,ic1,k,idir
c*********************************************************************** 
c***********************************************************************
c trace higher order states at cell edges
c***********************************************************************
c     write(6,*) "Entering chartracing2d1"
c     call flush(6)
c
      idir = 1
trace_call(1,0)dnl
c
c***********************************************************************
c     write(6,*) "leaving chartracing2d"
c      call flush(6)
      return
      end

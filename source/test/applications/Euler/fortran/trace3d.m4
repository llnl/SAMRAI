c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for trace for 3d euler equations.
c
define(NDIM,3)dnl
define(NEQU,5)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl
include(FORTDIR/m4trace3d.i)dnl

      subroutine computesound3d(ifirst0,ilast0,ifirst1,ilast1,
     &  ifirst2,ilast2,gamma,density,velocity,pressure,sound)
c***********************************************************************
      implicit none
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
c variables indexed as 3dimensional
      REAL
     &      gamma,
     &      density(CELL3d(ifirst,ilast,CELLG)),
     &      velocity(CELL3d(ifirst,ilast,CELLG),0:NDIM-1),
     &      pressure(CELL3d(ifirst,ilast,CELLG)),
     &      sound(CELL3d(ifirst,ilast,CELLG))
c
      integer ic0,ic1,ic2
c
c***********************************************************************
c
!$OMP PARALLEL DEFAULT(none)
!$OMPc SHARED(density,velocity,pressure,sound,
!$OMPc        ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
!$OMPc        gamma,
!$OMPc        CELLG,FLUXG,FACEG)
!$OMPc PRIVATE(ic0,ic1,ic2)

!$OMP DO SCHEDULE(DYNAMIC)
      do  ic0=ifirst0-CELLG,ilast0+CELLG
        do  ic1=ifirst1-CELLG,ilast1+CELLG
          do  ic2=ifirst2-CELLG,ilast2+CELLG
            sound(ic0,ic1,ic2) = sqrt(max(smallr,
     &                                    gamma*pressure(ic0,ic1,ic2)
     &                                    /density(ic0,ic1,ic2)))
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      return
      end 
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
c
      subroutine inittraceflux3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  density,velocity,pressure,
     &  tracelft0,tracelft1,tracelft2,
     &  tracergt0,tracergt1,tracergt2,
     &  fluxriem0,fluxriem1,fluxriem2)
c***********************************************************************
      implicit none
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      REAL
     &       density(CELL3d(ifirst,ilast,CELLG)),
     &      velocity(CELL3d(ifirst,ilast,CELLG),NDIM),
     &      pressure(CELL3d(ifirst,ilast,CELLG)),
     &     fluxriem0(FACE3d0(ifirst,ilast,FLUXG),NEQU),
     &     tracelft0(FACE3d0(ifirst,ilast,FACEG),NEQU),
     &     tracergt0(FACE3d0(ifirst,ilast,FACEG),NEQU),
     &     fluxriem1(FACE3d1(ifirst,ilast,FLUXG),NEQU),
     &     tracelft1(FACE3d1(ifirst,ilast,FACEG),NEQU),
     &     tracergt1(FACE3d1(ifirst,ilast,FACEG),NEQU),
     &     fluxriem2(FACE3d2(ifirst,ilast,FLUXG),NEQU),
     &     tracelft2(FACE3d2(ifirst,ilast,FACEG),NEQU),
     &     tracergt2(FACE3d2(ifirst,ilast,FACEG),NEQU)
      integer ic0,ic1,k,ie0,ie1,ic2,ie2
c***********************************************************************
c initialize left and right states at cell edges
c (first-order upwind)
c***********************************************************************
      if (FLUXG.lt.1) then
         write(6,*) "flux ghosts < 1!"
         stop
      endif
c
c***********************************************************************

!$OMP PARALLEL DEFAULT(none)
!$OMPc SHARED(density,velocity,pressure,
!$OMPc        tracelft0,tracelft1,tracelft2,
!$OMPc        tracergt0,tracergt1,tracergt2,
!$OMPc        fluxriem0,fluxriem1,fluxriem2,
!$OMPc        ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
!$OMPc        CELLG,FLUXG,FACEG)
!$OMPc PRIVATE(ic0,ic1,ic2,ie0,ie1,ie2)

!$OMP DO SCHEDULE(DYNAMIC)
trace_init(0,1,2,`ie0-1,ic1,ic2',`ie0,ic1,ic2')dnl
!$OMP END DO

!$OMP DO SCHEDULE(DYNAMIC)
trace_init(1,2,0,`ic0,ie1-1,ic2',`ic0,ie1,ic2')dnl
!$OMP END DO

!$OMP DO SCHEDULE(DYNAMIC)
trace_init(2,0,1,`ic0,ic1,ie2-1',`ic0,ic1,ie2')dnl
!$OMP END DO

c
c     we initialize the flux to be 0

!$OMP DO SCHEDULE(DYNAMIC)
      do ic2=ifirst2-FLUXG,ilast2+FLUXG
         do ic1=ifirst1-FLUXG,ilast1+FLUXG
           do ie0=ifirst0-FLUXG,ilast0+FLUXG+1
             do k=1,NEQU  
               fluxriem0(ie0,ic1,ic2,k) = zero
             enddo  
           enddo
         enddo
      enddo
!$OMP END DO
c
!$OMP DO SCHEDULE(DYNAMIC)
      do ic2=ifirst2-FLUXG,ilast2+FLUXG
         do ic0=ifirst0-FLUXG,ilast0+FLUXG
            do ie1=ifirst1-FLUXG,ilast1+FLUXG+1
               do k=1,NEQU  
                 fluxriem1(ie1,ic2,ic0,k) = zero
               enddo  
            enddo
         enddo
      enddo
!$OMP END DO
c
!$OMP DO SCHEDULE(DYNAMIC)
      do ic1=ifirst1-FLUXG,ilast1+FLUXG
         do ic0=ifirst0-FLUXG,ilast0+FLUXG
            do ie2=ifirst2-FLUXG,ilast2+FLUXG+1
               do k=1,NEQU  
                 fluxriem2(ie2,ic0,ic1,k) = zero
               enddo  
            enddo
         enddo
      enddo
!$OMP END DO
c
!$OMP END PARALLEL
      return
      end 
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
c
      subroutine chartracing3d0(dt,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  mc,dx,
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
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer mc,igdnv
      REAL dt 
c variables in 1d axis indexed
      REAL 
     &     dx, gamma
      REAL 
     &     sound(CELL3d(ifirst,ilast,CELLG))
      REAL 
     &     tracelft(FACE3d0(ifirst,ilast,FACEG),NEQU),
     &     tracergt(FACE3d0(ifirst,ilast,FACEG),NEQU)
c  side variables ifirst0 to ifirst0+mc plus ghost cells
      REAL
     & ttedgslp(ifirst0-FACEG:ifirst0+mc+FACEG,NEQU),
     & ttraclft(ifirst0-FACEG:ifirst0+mc+FACEG,NEQU),
     & ttracrgt(ifirst0-FACEG:ifirst0+mc+FACEG,NEQU)
c  cell variables ifirst0 to ifirst0+mc-1 plus ghost cells
      REAL 
     &  ttcelslp(ifirst0-CELLG:ifirst0+mc-1+CELLG,NEQU),
     &  ttsound(ifirst0-CELLG:ifirst0+mc-1+CELLG)
c***********************************************************************
      integer ic0,ic1,ic2,k,idir

c***********************************************************************
c trace higher order states at cell edges
c***********************************************************************
      idir = 0
!$OMP PARALLEL DEFAULT(none)
!$OMPc SHARED(dt,mc,dx,gamma,igdnv,sound,
!$OMPc        ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
!$OMPc        tracelft,tracergt,
!$OMPc        CELLG,FLUXG,FACEG,
!$OMPc        idir)
!$OMPc PRIVATE(ic0,ic1,ic2,
!$OMPc         ttsound,ttcelslp,ttedgslp,ttraclft,ttracrgt,
!$OMPc         k)
!$OMP DO SCHEDULE(DYNAMIC)
trace_call(0,1,2)dnl
!$OMP END DO
c***********************************************************************
!$OMP END PARALLEL
      return
      end
c***********************************************************************
c***********************************************************************
c***********************************************************************
      subroutine chartracing3d1(dt,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  mc,dx,
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
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer mc,igdnv
      REAL dt 
c variables in 1d axis indexed
      REAL 
     &     dx, gamma
      REAL 
     &     sound(CELL3d(ifirst,ilast,CELLG))
      REAL 
     &     tracelft(FACE3d1(ifirst,ilast,FACEG),NEQU),
     &     tracergt(FACE3d1(ifirst,ilast,FACEG),NEQU)
      REAL
     & ttedgslp(ifirst1-FACEG:ifirst1+mc+FACEG,NEQU),
     & ttraclft(ifirst1-FACEG:ifirst1+mc+FACEG,NEQU),
     & ttracrgt(ifirst1-FACEG:ifirst1+mc+FACEG,NEQU)
c  cell variables ifirst1 to ifirst1+mc-1 plus ghost cells
      REAL 
     &  ttcelslp(ifirst1-CELLG:ifirst1+mc-1+CELLG,NEQU),
     &  ttsound(ifirst1-CELLG:ifirst1+mc-1+CELLG)
c***********************************************************************
      integer ic0,ic1,ic2,k,idir
c***********************************************************************
c trace higher order states at cell edges
c***********************************************************************
      idir = 1
!$OMP PARALLEL DEFAULT(none)
!$OMPc SHARED(dt,mc,dx,gamma,igdnv,sound,
!$OMPc        ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
!$OMPc        tracelft,tracergt,
!$OMPc        CELLG,FLUXG,FACEG,
!$OMPc        idir)
!$OMPc PRIVATE(ic0,ic1,ic2,
!$OMPc         ttsound,ttcelslp,ttedgslp,ttraclft,ttracrgt,
!$OMPc         k)
!$OMP DO SCHEDULE(DYNAMIC)
trace_call(1,2,0)dnl
!$OMP END DO
c
c***********************************************************************
!$OMP END PARALLEL
      return
      end
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
      subroutine chartracing3d2(dt,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  mc,dx,
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
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer mc,igdnv
      REAL dt 
c variables in 1d axis indexed
      REAL 
     &     dx, gamma
      REAL 
     &     sound(CELL3d(ifirst,ilast,CELLG))
      REAL 
     &     tracelft(FACE3d2(ifirst,ilast,FACEG),NEQU),
     &     tracergt(FACE3d2(ifirst,ilast,FACEG),NEQU)
      REAL
     & ttedgslp(ifirst2-FACEG:ifirst2+mc+FACEG,NEQU),
     & ttraclft(ifirst2-FACEG:ifirst2+mc+FACEG,NEQU),
     & ttracrgt(ifirst2-FACEG:ifirst2+mc+FACEG,NEQU)
c  cell variables ifirst2 to ifirst2+mc-1 plus ghost cells
      REAL 
     &  ttcelslp(ifirst2-CELLG:ifirst2+mc-1+CELLG,NEQU),
     &  ttsound(ifirst2-CELLG:ifirst2+mc-1+CELLG)
c***********************************************************************
      integer ic0,ic1,ic2,k,idir
c***********************************************************************
c trace higher order states at cell edges
c***********************************************************************
      idir = 2
!$OMP PARALLEL DEFAULT(none)
!$OMPc SHARED(dt,mc,dx,gamma,igdnv,sound,
!$OMPc        ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
!$OMPc        tracelft,tracergt,
!$OMPc        CELLG,FLUXG,FACEG,
!$OMPc        idir)
!$OMPc PRIVATE(ic0,ic1,ic2,
!$OMPc         ttsound,ttcelslp,ttedgslp,ttraclft,ttracrgt,
!$OMPc         k)
!$OMP DO SCHEDULE(DYNAMIC)
trace_call(2,0,1)dnl
!$OMP END DO
c
c***********************************************************************
!$OMP END PARALLEL
      return
      end

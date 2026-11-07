c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for updating 3d flux sums from fluxes.
c
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl
c
c***********************************************************************
c Add flux integrals to fluxsums
c***********************************************************************
c
define(upfluxsumface_case_3d,`dnl
      subroutine upfluxsumface3d$1(
     &  ilo0,ilo1,ilo2,ihi0,ihi1,ihi2,
     &  flxgc0,flxgc1,flxgc2,
     &  iface,
     &  flux,fluxsum)
c***********************************************************************
      implicit none
c
      integer
     &  ilo0,ilo1,ilo2,ihi0,ihi1,ihi2,
     &  flxgc0,flxgc1,flxgc2,
     &  iface
      REAL
     &  flux(FACE3d$1VECG(ilo,ihi,flxgc)),
     &  fluxsum(OUTERFACE3d$1(ilo,ihi,0))
      integer ie$1,ic$2,ic$3
c
c***********************************************************************
c
      if (iface.eq.0) then
        ie$1 = ilo$1
      else
        ie$1 = ihi$1+1
      endif 
 
      do ic$3=ilo$3,ihi$3
         do ic$2=ilo$2,ihi$2
            fluxsum(ic$2,ic$3)=fluxsum(ic$2,ic$3)+flux(ie$1,ic$2,ic$3)
         enddo
      enddo
c
      return
      end
')dnl
c
c
define(upfluxsumside_case_3d,`dnl
      subroutine upfluxsumside3d$1(
     &  ilo0,ilo1,ilo2,ihi0,ihi1,ihi2,
     &  flxgc0,flxgc1,flxgc2,
     &  iside,
     &  flux,fluxsum)
c***********************************************************************
      implicit none
c
      integer
     &  ilo0,ilo1,ilo2,ihi0,ihi1,ihi2,
     &  flxgc0,flxgc1,flxgc2,
     &  iside
      REAL
     &  flux(SIDE3d$1VECG(ilo,ihi,flxgc)),
     &  fluxsum(OUTERSIDE3d$1(ilo,ihi,0))
      integer ic$1,ic$2,ic$3
c
c***********************************************************************
c
      if (iside.eq.0) then
        ic$1 = ilo$1
      else
        ic$1 = ihi$1+1
      endif 
 
      do ic$3=ilo$3,ihi$3
         do ic$2=ilo$2,ihi$2
            fluxsum(ic$2,ic$3)=fluxsum(ic$2,ic$3)+flux(ic0,ic1,ic2)
         enddo
      enddo
c
      return
      end
')dnl
upfluxsumface_case_3d(0,1,2)dnl
c
upfluxsumface_case_3d(1,2,0)dnl
c
upfluxsumface_case_3d(2,0,1)dnl
c
upfluxsumside_case_3d(0,1,2)dnl
c
upfluxsumside_case_3d(1,0,2)dnl
c
upfluxsumside_case_3d(2,0,1)dnl

c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for updating 2d flux sums from fluxes.
c
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl
c
c***********************************************************************
c Add flux integrals to fluxsums
c***********************************************************************
c
define(upfluxsumface_case_2d,`dnl
      subroutine upfluxsumface2d$1(
     &  ilo0,ilo1,ihi0,ihi1,
     &  flxgc0,flxgc1,
     &  iface,
     &  flux,fluxsum)
c***********************************************************************
      implicit none
c
      integer
     &  ilo0,ilo1,ihi0,ihi1,
     &  flxgc0,flxgc1,
     &  iface
      REAL
     &  flux(FACE2d$1VECG(ilo,ihi,flxgc)),
     &  fluxsum(OUTERFACE2d$1(ilo,ihi,0))
      integer ie$1,ic$2
c
c***********************************************************************
c
      if (iface.eq.0) then
        ie$1 = ilo$1
      else 
        ie$1 = ihi$1+1
      endif 
 
      do ic$2=ilo$2,ihi$2
         fluxsum(ic$2)=fluxsum(ic$2)+flux(ie$1,ic$2)
      enddo
c
      return
      end
')dnl
define(upfluxsumside_case_2d,`dnl
      subroutine upfluxsumside2d$1(
     &  ilo0,ilo1,ihi0,ihi1,
     &  flxgc0,flxgc1,
     &  iface,
     &  flux,fluxsum)
c***********************************************************************
      implicit none
c
      integer
     &  ilo0,ilo1,ihi0,ihi1,
     &  flxgc0,flxgc1,
     &  iface
      REAL
     &  flux(SIDE2d$1VECG(ilo,ihi,flxgc)),
     &  fluxsum(OUTERFACE2d$1(ilo,ihi,0))
      integer ic$1,ic$2
c
c***********************************************************************
c
      if (iface.eq.0) then
        ic$1 = ilo$1
      else 
        ic$1 = ihi$1+1
      endif 
 
      do ic$2=ilo$2,ihi$2
         fluxsum(ic$2)=fluxsum(ic$2)+flux(ic0,ic1)
      enddo
c
      return
      end
')dnl
upfluxsumface_case_2d(0,1)dnl
c
upfluxsumface_case_2d(1,0)dnl
c
upfluxsumside_case_2d(0,1)dnl
c
upfluxsumside_case_2d(1,0)dnl
c

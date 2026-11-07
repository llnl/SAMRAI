c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine for trace in 3d.
c
define(NEQU,1)dnl
define(REAL,`double precision')dnl

      subroutine trace3d(dt,ifirst,ilast,mc,
     &  dx,dir,advecspeed,igdnv,
     &  tracelft,tracergt,
     &  celslope,edgslope)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(FORTDIR/const.i)dnl
include(FORTDIR/probparams.i)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL dt
      integer ifirst,ilast,mc,dir,igdnv
      REAL dx,advecspeed
c
      REAL
     &  celslope(ifirst-CELLG:ifirst+mc-1+CELLG),
c       side variables
     &  tracelft(ifirst-FACEG:ifirst+mc+FACEG),
     &  tracergt(ifirst-FACEG:ifirst+mc+FACEG),
     &  edgslope(ifirst-FACEG:ifirst+mc+FACEG)
c
      integer ie,ic
      REAL bound,coef4,slope2,dtdx,slope4
      REAL du,nu
c
c***********************************************************************
c     ******************************************************************
c     * check for inflection points in characteristic speeds
c     * compute slopes at cell edges
c     * zero slopes if neighboring cells have different loading
c     ******************************************************************
      do ie=ifirst+1-FACEG,ilast+FACEG
          edgslope(ie)=tracergt(ie)-tracelft(ie)
      enddo
c     ******************************************************************
c     * limit slopes
c     ******************************************************************
        do ic=ifirst+1-CELLG,ilast+CELLG-1
          celslope(ic)=zero
        enddo
      if (igdnv.eq.2) then
c       write(6,*) "second-order slopes"
c       ****************************************************************
c       * second-order slopes
c       ****************************************************************
          do ic=ifirst+1-CELLG,ilast+CELLG-1
            slope2=half*(edgslope(ic)+edgslope(ic+1))
            celslope(ic)=half*(edgslope(ic)+edgslope(ic+1))
            if (edgslope(ic)*edgslope(ic+1).le.zero) then
              celslope(ic)=zero
            else
              bound=min(abs(edgslope(ic)),abs(edgslope(ic+1)))
              celslope(ic)=sign(min(two*bound,abs(slope2)),slope2)
            endif
c            write(6,*) "ic,celslope= ", ic,celslope(ic)
c            call flush(6)
          enddo
      else if (igdnv.eq.4) then
c       write(6,*) "fourth-order slopes"
c       ****************************************************************
c       * fourth-order slopes
c       ****************************************************************
          do ic=ifirst+2-CELLG,ilast+CELLG-2
            slope4=fourth*(tracergt(ic+2)-tracergt(ic-2))
            celslope(ic)=half*(edgslope(ic)+edgslope(ic+1))
            coef4=third*(four*celslope(ic)-slope4)
            if (edgslope(ic)*edgslope(ic+1).le.zero .or.
     &      coef4*celslope(ic).lt.zero) then
              celslope(ic)=zero
            else
              bound=min(abs(edgslope(ic)),abs(edgslope(ic+1)))
              celslope(ic)=sign(min(two*bound,abs(coef4)),coef4)
            endif
          enddo
      endif
c     ******************************************************************
c     * characteristic projection
c     ******************************************************************
      dtdx=dt/dx
      nu = advecspeed*dtdx
      do ic=ifirst-1,ilast+1
        du = celslope(ic)
        if (nu.gt.0) then
          tracelft(ic+1) = tracelft(ic+1) + half*(one -nu)*du
          tracergt(ic) = tracergt(ic)
        else
          tracelft(ic+1) = tracelft(ic+1)
          tracergt(ic) = tracergt(ic) -half*(one + nu)*du
        endif
      enddo
      return
      end

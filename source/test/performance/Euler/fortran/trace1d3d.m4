c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine for 3d trace.
c
define(NEQU,5)dnl
define(REAL,`double precision')dnl

      subroutine trace3d(dt,ifirst,ilast,mc,
     &  dx,dir,igdnv,sound,
     &  tracelft,tracergt, 
     &  celslope,edgslope)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL dt
      integer ifirst,ilast,mc,dir,igdnv
      REAL dx, flattn
c
      REAL  
     &  celslope(ifirst-CELLG:ifirst+mc-1+CELLG,NEQU),
     &  sound  (ifirst-CELLG:ifirst+mc-1+CELLG),
c       side variables
     &  tracelft(ifirst-FACEG:ifirst+mc+FACEG,NEQU),
     &  tracergt(ifirst-FACEG:ifirst+mc+FACEG,NEQU),
     &  edgslope(ifirst-FACEG:ifirst+mc+FACEG,NEQU)
c
      integer ie,ic,i,k
      REAL bound,coef4,slope2,dtdx,slope4
      REAL csq, rho, u, v, w, p
      REAL drho, du, dv, dw, dp
      REAL spminus, spplus, spzero, smminus, smplus, smzero
      REAL alpham, alphap, alpha0r
      REAL alpha0v,alpha0w
      REAL apright, amright, azrright
      REAL azv1rght,azw1rght
      REAL apleft, amleft, azrleft
      REAL azv1left,azw1left,ceigv
      integer dir1,dir2,dir3

      if (dir.eq.0) then
        dir1 = 2
        dir2 = 3
        dir3 = 4
      else if (dir.eq.1) then
        dir1 = 3
        dir2 = 4
        dir3 = 2
      else if (dir.eq.2) then
        dir1 = 4
        dir2 = 2
        dir3 = 3
      endif
c
c***********************************************************************
c     ******************************************************************
c     * check for inflection points in characteristic speeds
c     * compute slopes at cell edges
c     * zero slopes if neighboring cells have different loading
c     ******************************************************************
      do ie=ifirst+1-FACEG,ilast+FACEG
        do k=1,NEQU
          edgslope(ie,k)=tracergt(ie,k)-tracelft(ie,k)
        enddo
      enddo
c     ******************************************************************
c     * limit slopes
c     ******************************************************************
      do i=1,NEQU
        do ic=ifirst-CELLG,ilast+CELLG
          celslope(ic,i)=zero
        enddo
      enddo
      if (igdnv.eq.2) then
c       ****************************************************************
c       * second-order slopes
c       ****************************************************************
        do ic=ifirst+1-CELLG,ilast+CELLG-1
c         call flaten3d(ifirst,ilast,ic,mc,dir,tracergt,sound,flattn)
          flattn=one
          do i=1,NEQU
            slope2=half*(edgslope(ic,i)+edgslope(ic+1,i))
            celslope(ic,i)=half*(edgslope(ic,i)+edgslope(ic+1,i))
            if (edgslope(ic,i)*edgslope(ic+1,i).le.zero) then
              celslope(ic,i)=zero
            else
              bound=min(abs(edgslope(ic,i)),abs(edgslope(ic+1,i)))
              celslope(ic,i)=sign(min(two*bound,abs(slope2)),slope2)
              celslope(ic,i)=flattn*celslope(ic,i)
            endif
          enddo
        enddo
      else if (igdnv.eq.4) then
c       ****************************************************************
c       * fourth-order slopes
c       ****************************************************************
        do ic=ifirst+2-CELLG,ilast+CELLG-2
c         call flaten3d(ifirst,ilast,ic,mc,dir,tracergt,sound,flattn)
          flattn=one
          do i=1,NEQU
            slope4=fourth*(tracergt(ic+2,i)-tracergt(ic-2,i))
            celslope(ic,i)=half*(edgslope(ic,i)+edgslope(ic+1,i))
            coef4=third*(four*celslope(ic,i)-slope4)
            if (edgslope(ic,i)*edgslope(ic+1,i).le.zero .or.
     &      coef4*celslope(ic,i).lt.zero) then
              celslope(ic,i)=zero
            else
              bound=min(abs(edgslope(ic,i)),abs(edgslope(ic+1,i)))
              celslope(ic,i)=sign(min(two*bound,abs(coef4)),coef4)
              celslope(ic,i)=flattn*celslope(ic,i)
            endif
          enddo
        enddo
      endif
c     ******************************************************************
c     * characteristic projection
c     ******************************************************************
      do ic=ifirst-FACEG+1,ilast+FACEG
        dtdx=dt/dx
        rho = max(smallr,tracergt(ic,1))
        u = tracergt(ic,dir1)
        v = tracergt(ic,dir2)
        p = tracergt(ic,NEQU)
        drho = celslope(ic,1)
        du = celslope(ic,dir1)
        dv = celslope(ic,dir2)
        dp = celslope(ic,NEQU)
        w = tracergt(ic,dir3)
        dw = celslope(ic,dir3)
        alpha0w = dw
 
        ceigv = sound(ic)
        csq = ceigv**2
 
        alpham  = half*(dp/(rho*ceigv) - du)*rho/ceigv
        alphap  = half*(dp/(rho*ceigv) + du)*rho/ceigv
        alpha0r = drho - dp/csq
        alpha0v = dv
 
        if ((u-ceigv).gt.0) then
           spminus = -one
           smminus = (u-ceigv)*dtdx
        else
           spminus = (u-ceigv)*dtdx
           smminus = one
        endif
        if ((u+ceigv).gt.0) then
           spplus = -one
           smplus = (u+ceigv)*dtdx
        else
           spplus = (u+ceigv)*dtdx
           smplus = one
        endif
        if ((u).gt.0) then
           spzero = -one
           smzero = u*dtdx
        else
           spzero = u*dtdx
           smzero = one
        endif
        apright  = half*(-one - spplus )*alphap
        amright  = half*(-one - spminus)*alpham
        azrright = half*(-one - spzero )*alpha0r
        azv1rght = half*(-one - spzero )*alpha0v
        azw1rght = half*(-one - spzero )*alpha0w
        tracergt(ic,1) = rho + apright + amright + azrright
        tracergt(ic,1) = max(smallr,tracergt(ic,1))
        tracergt(ic,dir1) = u + (apright - amright)*ceigv/rho
        tracergt(ic,dir2) = v + azv1rght
        tracergt(ic,dir3) = w + azw1rght
        tracergt(ic,NEQU) = p + (apright + amright)*csq
 
        apleft  = half*(one - smplus )*alphap
        amleft  = half*(one - smminus)*alpham
        azrleft = half*(one - smzero )*alpha0r
        azv1left = half*(one - smzero )*alpha0v
        azw1left = half*(one - smzero )*alpha0w
        tracelft(ic+1,1) = rho + apleft + amleft + azrleft
        tracelft(ic+1,1) = max(smallr,tracelft(ic+1,1))
        tracelft(ic+1,dir1) = u + (apleft - amleft)*ceigv/rho
        tracelft(ic+1,dir2) = v + azv1left
        tracelft(ic+1,dir3) = w + azw1left
        tracelft(ic+1,NEQU) = p + (apleft + amleft)*csq
      enddo
      return
      end

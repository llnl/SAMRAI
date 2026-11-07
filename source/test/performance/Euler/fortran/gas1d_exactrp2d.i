c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file defining routine to compute solution to
c                Riemann problem in 2d.
c
define(rarefvel,`dnl
      $3= ($1*two/(gamma-one))*(one- exp(log($2)*tau))
')dnl

define(shockvel,`dnl
      $3= $1*alpha*(one-$2)/sqrt(1+ $2*beta)
')dnl

define(derrvel,`dnl
      $3= (two*tau*$1/(one-gamma))*exp(log($2)*(tau-one))
')dnl

define(dersvel,`dnl
      $3=  (-1)*$1*half*alpha*(two+beta*(one+$2))
     &           /(sqrt(one+ $2*beta)*(one+ $2*beta))
')dnl

c***********************************************************************
c***********************************************************************
c***********************************************************************
      subroutine gas1dexactrp2d(gamma,minpres,wleft,wright,striem)
c***********************************************************************
c description of arguments:
c   input:
c     wleft  ==> left state
c     wright ==> right state
c                components ordered: density, vn, pressure
c   output:
c     striem (density,velocity,pressure)
c             <== solution to riemann problem)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(FORTDIR/const.i)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL wleft(3),wright(3),striem(3)
      REAL gamma,minpres
c
      REAL pstar,ustarl,ustarr,rstarl,rstarr
      REAL px,pxl,pxr,pstarl,pstarr
      REAL ustar,sein,ro,uo,po,co,wo
      REAL feval,dfeval
      REAL cxl,cxr
      REAL dummy, rstar,cstar,wso,wsi,ushok,t4,t3,t5
      REAL gp1g2i,tau,beta,alpha
      REAL lpres,lvel,lrho,rpres,rvel,rrho
      REAL fleft,fright,dfleft,dfright,temp
      integer it,leftwave,rightwave
c
      tau=half*(gamma-one)/gamma
      gp1g2i=half*(gamma+one)/gamma
      alpha=sqrt(two/(gamma*(gamma-one)))
      beta=(gamma+one)/(gamma-one)
      lrho  = wleft(1)
      lvel  = wleft(2)
      lpres = wleft(3)
      rrho  = wright(1)
      rvel  = wright(2)
      rpres = wright(3)

      temp = abs(lvel-rvel)
      temp = max(abs(lpres-rpres),temp)
      if (temp.lt.smallr) then
         striem(2)=half*(lvel+rvel)
         striem(3)=half*(lpres+rpres)
         if (striem(2).gt.zero) then
            striem(1)=lrho
         else
            striem(1)=rrho
         endif
         goto 666
      endif

      cxl = sqrt(gamma*lpres/lrho)
      cxr = sqrt(gamma*rpres/rrho)

c  TERMINAL VELOCITIES

      ustarl = lvel + two*cxl/(gamma-1)
      ustarr = rvel - two*cxr/(gamma-1)

      if (ustarl.lt.ustarr) then
         write(6,*)" VACUUM in RS: "
         write(6,*)"     left state  ",lrho,lvel,lpres
         write(6,*)"     right state ",rrho,rvel,rpres
         write(6,*)"     vacuum vels ",ustarl,ustarr
         if (ustarl.gt.zero) then
            ustar = ustarl
            ro=lrho
            uo=lvel
            po=lpres
         else if (ustarr.lt.zero) then
            ustar = ustarr
            ro=rrho
            uo=rvel
            po=rpres
         else
            ustar = 0
            striem(1)=0
            striem(2)=half*(ustarl+ustarr)
            striem(3)=0
            goto 666
         endif
         pstar = 0
         rstar  = 0
         rstarl = 0
         rstarr = 0
         cstar  = 0
         sein=-dsign( one ,ustar)
         co=sqrt(po*gamma/ro)
         wso=sein*uo+co
         wsi=sein*ustar
         if (wso.lt.zero) then
            striem(1)=ro
            striem(2)=uo
            striem(3)=po
         else
            t4=(wso+wsi)/max(wso-wsi,wso+wsi,smallr)
            t4= half *(t4+ one )
            t5= one -t4
            t3=t5*co
            striem(1)=t5*ro
            striem(2)=t4*ustar+t5*uo
            striem(3)=t3*t3*striem(1)/gamma
         endif
         return
      endif
      
c-----------------------------------------------------------------
c        CLASSIFY the problem by 1 and  3 waves
c     shock= 1, rarefaction=2

      px = lpres/rpres
      if (px .lt. one ) then
rarefvel(cxr,px,temp)dnl
      else
shockvel(cxr,px,temp)dnl
      endif
      if (lvel .gt. (rvel - temp)) then
         leftwave = 1
      else
         leftwave = 2
      endif


      px = rpres/lpres
      if (px .lt. one ) then
rarefvel(cxl,px,temp)
      else
shockvel(cxl,px,temp)
      endif
      if ((lvel + temp) .gt. rvel ) then
         rightwave = 1
      else
         rightwave = 2
      endif

c-----------------------------------------------------------------
c     SOLVE RIEMANN PROBLEM

      if ((leftwave.eq.2).and.(rightwave.eq.1)) then
c       LEFT RAREF , RIGHT SHOCK
         pstar = rpres
         pxl = pstar/lpres
         pxr = pstar/rpres
rarefvel(cxl,pxl,fleft)dnl
shockvel(cxr,pxr,fright)dnl
derrvel(cxl,pxl,dfleft)dnl
dersvel(cxr,pxr,dfright)dnl
         do it=1,9
            feval = lvel + fleft - rvel + fright
            dfeval = dfleft/lpres + dfright/rpres
            pstar = pstar -feval/dfeval
            pxl = pstar/lpres
            pxr = pstar/rpres
            if (pstar.lt.minpres) then
               write(6,*)" pressure is too small or negative in
     &                 exactsolver"
               write(6,*)"   left state  ",lrho,lvel,lpres
               write(6,*)"   right state ",rrho,rvel,rpres
               write(6,*)"   lwave rwave = ",leftwave,rightwave
            endif
rarefvel(cxl,pxl,fleft)dnl
shockvel(cxr,pxr,fright)dnl
derrvel(cxl,pxl,dfleft)dnl
dersvel(cxr,pxr,dfright)dnl
         enddo
         pstar=max(pstar,minpres  )
         ustarl = lvel + fleft
         ustarr = rvel - fright
         rstarl = lrho*exp(log(pxl)/gamma)
         rstarr = rrho*(1+beta*pxr)/(beta+pxr)
      elseif ((leftwave.eq.1).and.(rightwave.eq.2)) then
c       LEFT SHOCK ,RIGHT RAREF
         pstar = lpres
         pxl = pstar/lpres
         pxr = pstar/rpres
shockvel(cxl,pxl,fleft)dnl
rarefvel(cxr,pxr,fright)dnl
dersvel(cxl,pxl,dfleft)dnl
derrvel(cxr,pxr,dfright)dnl
         do it=1,5
            feval = lvel + fleft - rvel + fright
            dfeval = dfleft/lpres + dfright/rpres
            pstar = pstar -feval/dfeval
            pxl = pstar/lpres
            pxr = pstar/rpres
            if (pstar.lt.minpres) then
               write(6,*)" pressure is too small or negative in
     &                 exactsolver"
               write(6,*)"   left state  ",lrho,lvel,lpres
               write(6,*)"   right state ",rrho,rvel,rpres
               write(6,*)"   lwave rwave = ",leftwave,rightwave
            endif
shockvel(cxl,pxl,fleft)dnl
rarefvel(cxr,pxr,fright)dnl
dersvel(cxl,pxl,dfleft)dnl
derrvel(cxr,pxr,dfright)dnl
         enddo
         pstar=max(pstar,minpres  )
         ustarl = lvel + fleft
         ustarr = rvel - fright
         rstarl = lrho*(1+beta*pxl)/(beta+pxl)
         rstarr = rrho*exp(log(pxr)/gamma)
      elseif ((leftwave.eq.1).and.(rightwave.eq.1)) then
c       RIGHT SHOCK, LEFT SHOCK
         pstar = max(rpres,lpres)
         pxl = pstar/lpres
         pxr = pstar/rpres
shockvel(cxl,pxl,fleft)dnl
shockvel(cxr,pxr,fright)dnl
dersvel(cxl,pxl,dfleft)dnl
dersvel(cxr,pxr,dfright)dnl
         do it=1,5
            feval = lvel + fleft - rvel + fright
            dfeval = dfleft/lpres + dfright/rpres
            pstar = pstar -feval/dfeval
            pxl = pstar/lpres
            pxr = pstar/rpres
            if (pstar.lt.minpres) then
               write(6,*)" pressure is too small or negative in
     &                 exactsolver"
               write(6,*)"   left state  ",lrho,lvel,lpres
               write(6,*)"   right state ",rrho,rvel,rpres
               write(6,*)"   lwave rwave = ",leftwave,rightwave
            endif
shockvel(cxl,pxl,fleft)dnl
shockvel(cxr,pxr,fright)dnl
dersvel(cxl,pxl,dfleft)dnl
dersvel(cxr,pxr,dfright)dnl
         enddo
         pstar=max(pstar,minpres  )
         ustarl = lvel + fleft
         ustarr = rvel - fright
         rstarl = lrho*(1+beta*pxl)/(beta+pxl)
         rstarr = rrho*(1+beta*pxr)/(beta+pxr)
      else
c       RIGHT RAREF, LEFT RAREF
         pstarl =  exp(tau*log(lpres))
         pstarr =  exp(tau*log(rpres))
         pstar = (half*(gamma-1)*(lvel-rvel)+cxl+cxr)/
     &           (cxr/pstarr +cxl/pstarl)
         pstar = exp(log(pstar)/tau)
         if (pstar.lt.minpres) then
            write(6,*)" pressure is too small or negative in
     &              exactsolver"
            write(6,*)"   left state  ",lrho,lvel,lpres
            write(6,*)"   right state ",rrho,rvel,rpres
            write(6,*)"   lwave rwave = ",leftwave,rightwave
         endif
         pstar=max(pstar,minpres  )
         pxl = pstar/lpres
         pxr = pstar/rpres
rarefvel(cxl,pxl,fleft)dnl
rarefvel(cxr,pxr,fright)dnl
         ustarl = lvel + fleft
         ustarr = rvel - fright
         rstarl = lrho*exp(log(pxl)/gamma)
         rstarr = rrho*exp(log(pxr)/gamma)
      endif
      pstarl = pstar
      pstarr = pstar
 
      ustar=half*(ustarl+ustarr)

c-----------------------------------------------------------------
c     PICK STATE FOR FLUXES

      sein=-dsign( one ,ustar)
      if (sein.ge.zero) then
          ro=rrho
          uo=rvel
          po=rpres
          rstar =rstarr
      else
          ro=lrho
          uo=lvel
          po=lpres
          rstar =rstarl
      endif
      cstar=sqrt(gamma*pstar/rstar)
      dummy=pstar-po
      if (dummy.ge.zero) then
c               shock
         wo=sqrt(gamma*ro*po*( one +gp1g2i*dummy/po))
         ushok=sein*ustar+wo/rstar
         if (ushok.ge.zero) then
            striem(1)=rstar
            striem(2)=ustar
            striem(3)=pstar
         else 
            striem(1)=ro
            striem(2)=uo
            striem(3)=po
         endif
      else
c               rarefaction
         co=sqrt(po*gamma/ro)
         wso=sein*uo+co
         wsi=sein*ustar+cstar
         if (wsi.ge.zero) then
            striem(1)=rstar
            striem(2)=ustar
            striem(3)=pstar
         else if (wso.lt.zero) then
            striem(1)=ro
            striem(2)=uo
            striem(3)=po
         else
            t4=(wso+wsi)/max(wso-wsi,wso+wsi,smallr)
            t4= half *(t4+ one )
            t5= one -t4
            t3=t4*cstar+t5*co
            striem(1)=t4*rstar+t5*ro
            striem(2)=t4*ustar+t5*uo
            striem(3)=t3*t3*striem(1)/gamma
         endif
      endif
c     ------------------------------------------------------------------
c     write(6,*) "leaving riemnv"
c     ------------------------------------------------------------------
 666  return
      end
c***********************************************************************
c***********************************************************************
c***********************************************************************

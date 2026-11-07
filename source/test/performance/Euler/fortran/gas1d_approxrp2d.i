c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file defining routine to compute solution to
c                Riemann problem in 2d.
c
      subroutine gas1dapproxrp2d(eosgam,wleft,wright, striem)
c***********************************************************************
c description of arguments:
c   input:
c     eosgam ==> ratio of specific heats
c     wleft  ==> left state
c     wright ==> right state
c                components ordered: density, vn, pressure
c   output:
c     striem (density,velocity,pressure) 
c             <== solution to riemann problem
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(FORTDIR/const.i)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL wleft(3),wright(3),eosgam,striem(3)
c
      REAL clsqp,clsqm,pstar,zp,zm,ustarm,ustarp,zsum
      REAL wsqm,wsqp,wm,wp,ustar,sein,ro,uo,po,co,wo
      REAL dummy, rstar,cstar,wso,wsi,ushok,t4,t3,t5
      REAL gp1g2i, gp1o2gm1
      integer it
      REAL separ
c
c     ******************************************************************
c     * the entries of wleft, wright, and wriem are (in order):
c     *   density, tangential velocity, normal velocity, pressure
c     ******************************************************************

      separ = zero
      do it=1,3
         separ = separ+abs(wleft(it) - wright(it))
      enddo
      if (separ.lt.smallr) then
         do it=1,3
            striem(it) = half*(wleft(it)+wright(it))
         enddo
         return
      endif
 

      gp1g2i=half*(eosgam+one)/eosgam
      gp1o2gm1=half*(eosgam+one)/(eosgam-one)
      clsqp=eosgam*wleft(1)*max(smallr,wleft(3))
      wp=sqrt(clsqp)
      clsqm=eosgam*wright(1)*max(smallr,wright(3))
      wm=sqrt(clsqm)
      pstar=(wp*wright(3)+wm*wleft(3)-wm*wp*(wright(2)-wleft(2)))
     &          /(wp+wm)
      pstar=max(pstar,smallr  )
c
      do it=1,5
          wsqp=clsqp*( one +gp1g2i*(pstar/wleft(3)- one ))
          wp=sqrt(wsqp)
          wsqm=clsqm*( one +gp1g2i*(pstar/wright(3)- one ))
          wm=sqrt(wsqm)
          zp= two *wp*wsqp/(wsqp+clsqp)
          zm= two *wm*wsqm/(wsqm+clsqm)
          ustarm=wright(2)-(wright(3)-pstar)/wm
          ustarp=wleft(2)+(wleft(3)-pstar)/wp
          zsum=zp+zm
          pstar=pstar-zp*zm*(ustarm-ustarp)/zsum
          pstar=max(pstar,smallr)
      enddo
c
  
        ustar=(zp*ustarm+zm*ustarp)/zsum
        sein=-dsign( one ,ustar)
        if (sein.ge.zero) then
          ro=wright(1)
          uo=wright(2)
          po=wright(3)
        else
          ro=wleft(1)
          uo=wleft(2)
          po=wleft(3)
        endif
        po=max(smallr,po)
        co=sqrt(po*eosgam/ro)
        dummy=pstar-po
        wo=eosgam*ro*po*( one +gp1g2i*dummy/po)
        rstar=ro/( one -ro*dummy/wo)
        wo=sqrt(wo)
        cstar=sqrt(eosgam*pstar/rstar)
        wso=sein*uo+co
        wsi=sein*ustar+cstar
        ushok=sein*ustar+wo/rstar
        if (dummy.ge.zero) then
          wsi=ushok
          wso=ushok
        else
          wsi=wsi
          wso=wso
        endif
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
          striem(2) =t4*ustar+t5*uo
          striem(3)=t3*t3*striem(1)/eosgam
        endif
  
      return
      end

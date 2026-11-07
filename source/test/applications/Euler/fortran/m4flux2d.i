c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for 2d flux calculation.
c
define(riemann_solve,`dnl
c     write(6,*) "checking onedr sol in riemann solve  "
c     write(6,*) "         dt= ",dt
c     write(6,*) "  calculating flux$1, 1+extra_cell= ",$3
c     write(6,*) "  ic$2=",ifirst$2-1-$3,ilast$2+1+$3
c     write(6,*) "  ie$1=",ifirst$1-1-$3,ilast$1+1+1+$3

      if (rpchoice.eq.APPROX_RIEM_SOLVE
     &    .or. rpchoice.eq.EXACT_RIEM_SOLVE ) then

      do ic$2=ifirst$2-$3,    
     &        ilast$2+$3
         do ie$1=ifirst$1-$3,
     &           ilast$1+1+$3
 
c           ************************************************************
c           *  Assign left and right states.  Note only normal vel used.
c           ************************************************************
            stateL(1) = trlft$1(ie$1,ic$2,1)
            stateL(2) = trlft$1(ie$1,ic$2,2+$1)
            stateL(3) = trlft$1(ie$1,ic$2,NEQU)
  
            stateR(1) = trrgt$1(ie$1,ic$2,1)
            stateR(2) = trrgt$1(ie$1,ic$2,2+$1)
            stateR(3) = trrgt$1(ie$1,ic$2,NEQU)

            if (rpchoice.eq.APPROX_RIEM_SOLVE) then
               call gas1dapproxrp2d(gamma,stateL,stateR,riemst)
            else 
               call gas1dexactrp2d(gamma,smallr,stateL,stateR,riemst)
            endif

            if (riemst(2).le.zero) then
               vel_tan=trrgt$1(ie$1,ic$2,2+$2)
            else
               vel_tan=trlft$1(ie$1,ic$2,2+$2)
            endif

            mom$1=riemst(1)*riemst(2)
            v2norm = riemst(2)**2+vel_tan**2
            Hent = riemst(3)/gam_min_one + v2norm*riemst(1)/two

            flux$1(ie$1,ic$2,1)= dt*mom$1
            flux$1(ie$1,ic$2,2+$1)= dt*(riemst(3)+mom$1*riemst(2))
            flux$1(ie$1,ic$2,2+$2)= dt*mom$1*vel_tan 
            flux$1(ie$1,ic$2,NEQU)= dt*riemst(2)*(Hent+riemst(3))

         enddo
      enddo

      elseif (rpchoice.eq.HLLC_RIEM_SOLVE) then

      do ic$2=ifirst$2-$3,
     &        ilast$2+$3
         do ie$1=ifirst$1-$3,
     &           ilast$1+1+$3
 
c        ************************************************************
c        *  Assign left and right states.
c        *  Note all vel comps used for Roe average.
c        ************************************************************
         do j=1,NEQU
            stateL(j) = trlft$1(ie$1,ic$2,j)
            stateR(j) = trrgt$1(ie$1,ic$2,j)
         enddo

c        ************************************************************
c        *  Calculate bounding signal speeds. To do this, need the
c        *  Roe-average of the velocity and sound speed.
c        ************************************************************
         w    = one / ( one + sqrt( stateR(1)/stateL(1) ) )
         omw  = one - w
         hat(2+$1) = w*stateL(2+$1) + omw*stateR(2+$1)
         aLsq = gamma * stateL(NEQU) / stateL(1)
         aRsq = gamma * stateR(NEQU) / stateR(1)
         hat(NEQU+1) = sqrt( w*aLsq + omw*aRsq
     &               + half*gam_min_one*w*omw*(
     &                  (stateR(2)-stateL(2))**2
     &                + (stateR(3)-stateL(3))**2 ) )

         sL = min(stateL(2+$1) - sqrt(aLsq), hat(2+$1) - hat(NEQU+1))
         sR = max(stateR(2+$1) + sqrt(aRsq), hat(2+$1) + hat(NEQU+1))
         mfL = stateL(1) * ( sL - stateL(2+$1) )
         mfR = stateR(1) * ( sR - stateR(2+$1) )
         sM  = ( stateR(NEQU) - stateL(NEQU)
     &     + mfL*stateL(2+$1) - mfR*stateR(2+$1) ) / ( mfL - mfR )

c        ************************************************************
c        *  Calculate flux starting at upwind state.
c        ************************************************************
         if ( sM.gt.zero ) then

c           *********************************************************
c           * Flow is to the right; start from left state.
c           *********************************************************
            flux(1)   = stateL(1) * stateL(2+$1)
            flux(2+$1) = flux(1) * stateL(2+$1) + stateL(NEQU)
            flux(3-$1) = flux(1) * stateL(3-$1)
            keL = half * (stateL(2)**2 + stateL(3)**2)
            flux(NEQU) = flux(1) * (aLsq / gam_min_one + keL)

c           *********************************************************
c           * Check if flow is subsonic.
c           *********************************************************
            if ( sL.lt.zero ) then

c              ******************************************************
c              * Add contribution from left acoustic wave.
c              ******************************************************
               denom = one / (sL-sM)
               star(1)    = stateL(1) * (sL-stateL(2+$1)) * denom
               star(2+$1) = sM
               star(3-$1) = stateL(3-$1)
               star(NEQU) = stateL(NEQU)*( one + gamma * denom
     &                                   * ( sM-stateL(2+$1) ) )
     &                    + half * gam_min_one
     &                        * star(1) * ( sM-stateL(2+$1) )**2

               diff(1) = star(1) - stateL(1)
               diff(2) = star(1)*star(2) - stateL(1)*stateL(2)
               diff(3) = star(1)*star(3) - stateL(1)*stateL(3)
               diff(NEQU) = ( star(NEQU)-stateL(NEQU) )
     &                 / gam_min_one
     &                 + half * star(1) * (star(2)**2 + star(3)**2)
     &                         - stateL(1) * keL

               do j=1,NEQU
                  flux(j) = flux(j) + sL*diff(j)
               enddo

            endif

         else

c           *********************************************************
c           *  Flow is to the left; start from right state.
c           *********************************************************
            flux(1)   = stateR(1) * stateR(2+$1)
            flux(2+$1) = flux(1) * stateR(2+$1) + stateR(NEQU)
            flux(3-$1) = flux(1) * stateR(3-$1)
            keR = half * (stateR(2)**2 + stateR(3)**2)
            flux(NEQU)   = flux(1) * (aRsq / gam_min_one + keR)

c           *********************************************************
c           * Check if flow is subsonic.
c           *********************************************************
            if ( sR.gt.zero ) then

c              ******************************************************
c              * Add contribution from right acoustic wave.
c              ******************************************************
               denom = one / (sR-sM)
               star(1)    = stateR(1) * (sR-stateR(2+$1)) * denom
               star(2+$1) = sM
               star(3-$1) = stateR(3-$1)
               star(NEQU) = stateR(NEQU) * (1 + gamma * denom
     &                                   * ( sM-stateR(2+$1) ) )
     &                    + half * gam_min_one
     &                           * star(1) * ( sM-stateR(2+$1) )**2

               diff(1) = star(1) - stateR(1)
               diff(2) = star(1)*star(2) - stateR(1)*stateR(2)
               diff(3) = star(1)*star(3) - stateR(1)*stateR(3)
               diff(NEQU) = ( star(NEQU)-stateR(NEQU) )
     &                 / gam_min_one
     &                 + half * star(1) * (star(2)**2 + star(3)**2)
     &                      - stateR(1) * keR

               do j=1,NEQU
                 flux(j) = flux(j) + sR*diff(j)
               enddo

            endif

         endif

c        ************************************************************
c        *  Assign average interface fluxes.
c        ************************************************************
         do j=1,NEQU
            flux$1(ie$1,ic$2,j) = dt * flux(j)
         enddo

         enddo
      enddo

      endif
')dnl
define(artificial_viscosity1,`dnl
           do ic$2=ifirst$2-(FLUXG-1),ilast$2+(FLUXG-1)
             do ie$1=ifirst$1-(FLUXG-1),ilast$1+(FLUXG)
               maxeig =trrgt$1(ie$1,ic$2,4)-trlft$1(ie$1,ic$2,4)
               vcoef = tenth*abs(maxeig)
     
               mom0L=trlft$1(ie$1,ic$2,1)*trlft$1(ie$1,ic$2,2)
               mom1L=trlft$1(ie$1,ic$2,1)*trlft$1(ie$1,ic$2,3)
               v2norm = mom1L**2+mom0L**2
               enerL=trlft$1(ie$1,ic$2,NEQU)/gam_min_one+
     &                v2norm/(trlft$1(ie$1,ic$2,1)*two)
               mom0R=trrgt$1(ie$1,ic$2,1)*trrgt$1(ie$1,ic$2,2)
               mom1R=trrgt$1(ie$1,ic$2,1)*trrgt$1(ie$1,ic$2,3)
               v2norm = mom1R**2+mom0R**2
               enerR=trrgt$1(ie$1,ic$2,NEQU)/gam_min_one+
     &                  v2norm/(trrgt$1(ie$1,ic$2,1)*two)
  
               vcorr(1) = dt*vcoef*
     &           (trrgt$1(ie$1,ic$2,1)-trlft$1(ie$1,ic$2,1))
               vcorr(2) = dt*vcoef*(mom0R-mom0L)
               vcorr(3) = dt*vcoef*(mom1R-mom1L)
               vcorr(NEQU) = dt*vcoef*(enerR-enerL)
               do j=1,NEQU
                 flux$1(ie$1,ic$2,j)= flux$1(ie$1,ic$2,j) 
     &                                     -vcorr(j)
               enddo
             enddo
           enddo
')dnl
define(artificial_viscosity2,`dnl
        do ic1=ifirst1-(FLUXG-1),ilast1+(FLUXG-1)
          do ie0=ifirst0-(FLUXG-1),ilast0+(FLUXG)
            maxeig =pressure(ie0,ic1)-pressure(ie0-1,ic1)
            vcoef = tenth*abs(maxeig)
            mom0L=density(ie0-1,ic1,1)*velocity(ie0-1,ic1,0)
            mom1L=density(ie0-1,ic1,1)*velocity(ie0-1,ic1,1)
            mom0R=density(ie0  ,ic1,1)*velocity(ie0  ,ic1,0)
                mom1R=density(ie0  ,ic1,1)*velocity(ie0  ,ic1,1)
            v2norm = mom1L**2+mom0L**2
            enerL=pressure(ie0-1,ic1)/gam_min_one+
     &               v2norm/(density(ie0-1,ic1)*two)
            v2norm = mom1R**2+mom0R**2
            enerR=pressure(ie0,ic1)/gam_min_one+
     &             v2norm/(density(ie0,ic1)*two)
            vcorr1= dt*vcoef*(density(ie0,ic1)-density(ie0-1,ic1))
            vcorr2= dt*vcoef*(mom0R-mom0L)
            vcorr3= dt*vcoef*(mom1R-mom1L)
            vcorr4= dt*vcoef*(enerR-enerL)
            flux0(ie0,ic1,1)= flux0(ie0,ic1,1) - vcorr1
            flux0(ie0,ic1,2)= flux0(ie0,ic1,2) - vcorr2
            flux0(ie0,ic1,3)= flux0(ie0,ic1,3) - vcorr3
            flux0(ie0,ic1,4)= flux0(ie0,ic1,4) - vcorr4
          enddo
        enddo
        do ic0=ifirst0,ilast0
          do ie1=ifirst1,ilast1+1
            maxeig =pressure(ic0,ie1)-pressure(ic0-1,ie1)
            vcoef = tenth*abs(maxeig)
            mom0L=density(ic0,ie1-1,1)*velocity(ic0,ie1-1,0)
            mom1L=density(ic0,ie1-1,1)*velocity(ic0,ie1-1,1)
            mom0R=density(ic0,ie1  ,1)*velocity(ic0,ie1  ,0)
            mom1R=density(ic0,ie1  ,1)*velocity(ic0,ie1  ,1)
            v2norm = mom1L**2+mom0L**2
            enerL=pressure(ic0,ie1-1)/gam_min_one+
     &             v2norm/(density(ic0,ie1-1)*two)
            v2norm = mom1R**2+mom0R**2
            enerR=pressure(ic0,ie1)/gam_min_one+
     &             v2norm/(density(ic0,ie1)*two)
            vcorr1= dt*vcoef*(density(ic0,ie1)-density(ic0,ie1-1))
            vcorr2= dt*vcoef*(mom0R-mom0L)
            vcorr3= dt*vcoef*(mom1R-mom1L)
            vcorr4= dt*vcoef*(enerR-enerL)
            flux1(ie1,ic0,1)= flux1(ie1,ic0,1) - vcorr1
            flux1(ie1,ic0,2)= flux1(ie1,ic0,2) - vcorr2
            flux1(ie1,ic0,3)= flux1(ie1,ic0,3) - vcorr3
            flux1(ie1,ic0,4)= flux1(ie1,ic0,4) - vcorr4
          enddo
        enddo
')dnl

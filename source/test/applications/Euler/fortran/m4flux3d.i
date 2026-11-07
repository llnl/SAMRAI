c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for 3d flux calculation.
c
define(riemann_solve,`dnl
c     write(6,*) "checking onedr sol in riemann solve  "
c     write(6,*) "         dt= ",dt

      if (rpchoice.eq.APPROX_RIEM_SOLVE
     &    .or. rpchoice.eq.EXACT_RIEM_SOLVE ) then
c     ************************************************************
c     *  Approximate Riemann solver and exact Riemann solver have
c     *  identical setup and post-process phases.
c     ************************************************************

!$OMP DO SCHEDULE(DYNAMIC)
      do ic$3=ifirst$3-$6,ilast$3+$6
         do ic$2=ifirst$2-$5,ilast$2+$5
           do ie$1=ifirst$1-(FLUXG-1),ilast$1+1+(FLUXG-1)

c          ************************************************************
c          *  Assign left and right states.  Note only normal vel used.
c          ************************************************************
              stateL(1) = trlft$1(ie$1,$4,1)
              stateL(2) = trlft$1(ie$1,$4,2+$1)
              stateL(3) = trlft$1(ie$1,$4,NEQU)

              stateR(1) = trrgt$1(ie$1,$4,1)
              stateR(2) = trrgt$1(ie$1,$4,2+$1)
              stateR(3) = trrgt$1(ie$1,$4,NEQU)

              if (rpchoice.eq.APPROX_RIEM_SOLVE) then
                 call gas1dapproxrp3d(gamma,stateL,stateR,riemst)
              else if (rpchoice.eq.EXACT_RIEM_SOLVE) then
                 call gas1dexactrp3d(gamma,smallr,stateL,stateR,riemst)
              endif

              if (riemst(2).le.zero) then
                 vel($2)=trrgt$1(ie$1,$4,2+$2)
                 vel($3)=trrgt$1(ie$1,$4,2+$3)
              else
                 vel($2)=trlft$1(ie$1,$4,2+$2)
                 vel($3)=trlft$1(ie$1,$4,2+$3)
              endif
              vel($1)=riemst(2)

              mom$1 =riemst(1)*vel($1)
              v2norm = vel(0)**2+vel(1)**2+vel(2)**2
              Hent =riemst(3)/gam_min_one+v2norm*riemst(1)/two

              flux$1(ie$1,$4,1)= dt*mom$1
              flux$1(ie$1,$4,2+$1)= dt*(mom$1*vel($1)+riemst(3))
              flux$1(ie$1,$4,2+$2)= dt*mom$1*vel($2)
              flux$1(ie$1,$4,2+$3)= dt*mom$1*vel($3)
              flux$1(ie$1,$4,5)= dt*riemst(2)*(Hent+riemst(3))

           enddo
         enddo
      enddo
!$OMP END DO

      elseif (rpchoice.eq.HLLC_RIEM_SOLVE) then
c     ******************************************************************
c     *  HLLC Riemann Solver
c     ******************************************************************

!$OMP DO SCHEDULE(DYNAMIC)
      do ic$3=ifirst$3-$6,ilast$3+$6
         do ic$2=ifirst$2-$5,ilast$2+$5
            do ie$1=ifirst$1-(FLUXG-1),ilast$1+1+(FLUXG-1)

c           ************************************************************
c           *  Assign left and right states.
c           ************************************************************
            do j=1,NEQU
               stateL(j) = trlft$1(ie$1,$4,j)
               stateR(j) = trrgt$1(ie$1,$4,j)
            enddo

c           ************************************************************
c           *  Calculate bounding signal speeds. To do this, need the
c           *  Roe-average of the velocity and sound speed.
c           ************************************************************
            w    = one / ( one + sqrt( stateR(1)/stateL(1) ) )
            omw  = one - w
            hat(2+$1) = w*stateL(2+$1) + omw*stateR(2+$1)
            aLsq = gamma * stateL(NEQU) / stateL(1)
            aRsq = gamma * stateR(NEQU) / stateR(1)
            hat(NEQU+1) = sqrt( w*aLsq + omw*aRsq
     &                 + half*gam_min_one*w*omw*(
     &                     (stateR(2)-stateL(2))**2
     &                   + (stateR(3)-stateL(3))**2
     &                   + (stateR(4)-stateL(4))**2 ) )

            sL  = min( stateL(2+$1) - sqrt(aLsq),
     &                    hat(2+$1) - hat(NEQU+1) )
            sR  = max( stateR(2+$1) + sqrt(aRsq),
     &                    hat(2+$1) + hat(NEQU+1) )
            mfL = stateL(1) * ( sL - stateL(2+$1) )
            mfR = stateR(1) * ( sR - stateR(2+$1) )
            sM  = ( stateR(NEQU) - stateL(NEQU)
     &             + mfL*stateL(2+$1) - mfR*stateR(2+$1) )/(mfL-mfR)

c           ************************************************************
c           *  Calculate flux starting at upwind state.
c           ************************************************************
            if ( sM.gt.zero ) then

c              *********************************************************
c              * Flow is to the right; start from left state.
c              *********************************************************
               flux(1)   = stateL(1) * stateL(2+$1)
               flux(2+$1) = flux(1) * stateL(2+$1) + stateL(NEQU)
               flux(2+$2) = flux(1) * stateL(2+$2)
               flux(2+$3) = flux(1) * stateL(2+$3)
               keL = half * (stateL(2)**2 + stateL(3)**2 +
     &                               stateL(4)**2)
               flux(NEQU) = flux(1)
     &                    * (aLsq / gam_min_one + keL)

c              *********************************************************
c              * Check if flow is subsonic.
c              *********************************************************
               if ( sL.lt.zero ) then

c                 ******************************************************
c                 * Add contribution from left acoustic wave.
c                 ******************************************************
                  denom = one / (sL-sM)
                  star(1)    = stateL(1) * (sL-stateL(2+$1)) * denom
                  star(2+$1) = sM
                  star(2+$2) = stateL(2+$2)
                  star(2+$3) = stateL(2+$3)
                  star(NEQU) = stateL(NEQU)*( one + gamma * denom
     &                                        * ( sM-stateL(2+$1) ) )
     &                       + half * gam_min_one
     &                              * star(1) * ( sM-stateL(2+$1) )**2

                  diff(1) = star(1) - stateL(1)
                  diff(2) = star(1)*star(2) - stateL(1)*stateL(2)
                  diff(3) = star(1)*star(3) - stateL(1)*stateL(3)
                  diff(4) = star(1)*star(4) - stateL(1)*stateL(4)
                  diff(NEQU) = ( star(NEQU)-stateL(NEQU) )
     &               / gam_min_one
     &               + half * star(1) *
     &                   (star(2)**2 + star(3)**2+ star(4)**2)
     &               - stateL(1) * keL

                  do j=1,NEQU
                    flux(j) = flux(j) + sL*diff(j)
                  enddo

               endif

            else

c              *********************************************************
c              *  Flow is to the left; start from right state.
c              *********************************************************
               flux(1)   = stateR(1) * stateR(2+$1)
               flux(2+$1) = flux(1) * stateR(2+$1) + stateR(NEQU)
               flux(2+$2) = flux(1) * stateR(2+$2)
               flux(2+$3) = flux(1) * stateR(2+$3)
               keR = half * (stateR(2)**2 + stateR(3)**2 +
     &                               stateR(4)**2)
               flux(NEQU) = flux(1)
     &                    * (aRsq / gam_min_one + keR)

c              *********************************************************
c              * Check if flow is subsonic.
c              *********************************************************
               if ( sR.gt.zero ) then

c                 ******************************************************
c                 * Add contribution from right acoustic wave.
c                 ******************************************************
                  denom = one / (sR-sM)
                  star(1)    = stateR(1) * (sR-stateR(2+$1)) * denom
                  star(2+$1) = sM
                  star(2+$2) = stateR(2+$2)
                  star(2+$3) = stateR(2+$3)
                  star(NEQU) = stateR(NEQU)*(1 + gamma * denom
     &                                      * ( sM-stateR(2+$1) ) )
     &                       + half * gam_min_one
     &                              * star(1) * ( sM-stateR(2+$1) )**2

                  diff(1) = star(1) - stateR(1)
                  diff(2) = star(1)*star(2) - stateR(1)*stateR(2)
                  diff(3) = star(1)*star(3) - stateR(1)*stateR(3)
                  diff(4) = star(1)*star(4) - stateR(1)*stateR(4)
                  diff(NEQU) = ( star(NEQU)-stateR(NEQU) )
     &                       / gam_min_one
     &                       + half * star(1) *
     &                         (star(2)**2 + star(3)**2+ star(4)**2)
     &                       - stateR(1) * keR

                  do j=1,NEQU
                     flux(j) = flux(j) + sR*diff(j)
                  enddo

               endif

            endif

c        ************************************************************
c        *  Assign average interface fluxes.
c        ************************************************************
            do j=1,NEQU
               flux$1(ie$1,$4,j) = dt * flux(j)
            enddo

            enddo
         enddo
      enddo
!$OMP END DO

      endif
')dnl
define(correc_flux2d,`dnl
c   correct the $1-direction with $3-fluxes
!$OMP DO SCHEDULE(DYNAMIC)
      do ic$5=ifirst$5-(FLUXG),ilast$5+(FLUXG)
         do ic$3=ifirst$3-(FLUXG-1),ilast$3+(FLUXG-1)
           do ic$1=ifirst$1-(FLUXG),ilast$1+(FLUXG)
             rho  = density(ic0,ic1,ic2)
             vel0 = velocity(ic0,ic1,ic2,0)
             vel1 = velocity(ic0,ic1,ic2,1)
             vel2 = velocity(ic0,ic1,ic2,2)
             v2norm= vel0**2+vel1**2 +vel2**2
             do k=1,NEQU
               trnsvers(k)=
     &           (flux$3(ic$3+1,$4,k)-flux$3(ic$3,$4,k))*0.5/dx($3)
             enddo
c
             ttv(1)= trnsvers(1)
             ttv(2)= (trnsvers(2) - vel0*trnsvers(1))/rho
             ttv(3)= (trnsvers(3) - vel1*trnsvers(1))/rho
             ttv(4)= (trnsvers(4) - vel2*trnsvers(1))/rho
             ttv(NEQU)= (v2norm*half*trnsvers(1) -
     &                vel0*trnsvers(2) - vel1*trnsvers(3) -
     &                vel2*trnsvers(4) +
     &                trnsvers(NEQU))*gam_min_one
             do k=1,NEQU
               ttracelft$1(ic$1+1,$2,k)=tracelft$1(ic$1+1,$2,k)
     &                                    - ttv(k)
               ttracergt$1(ic$1  ,$2,k)=tracergt$1(ic$1  ,$2,k)
     &                                    - ttv(k)
             enddo
           enddo
         enddo
      enddo
!$OMP END DO
')dnl
define(correc_flux3d,`dnl
c   correct the $1-direction with $2$3-fluxes
!$OMP DO SCHEDULE(DYNAMIC)
      do ic$1=ifirst$1-FLUXG,ilast$1+FLUXG
         do ic$3=ifirst$3-(FLUXG-1),ilast$3+(FLUXG-1)
           do ic$2=ifirst$2-(FLUXG-1),ilast$2+(FLUXG-1)
             rho  = density(ic0,ic1,ic2)
             vel0 = velocity(ic0,ic1,ic2,0)
             vel1 = velocity(ic0,ic1,ic2,1)
             vel2 = velocity(ic0,ic1,ic2,2)
             v2norm= vel0**2+vel1**2 +vel2**2
             do k=1,NEQU
               trnsvers(k)=0.5*(
     &          (flux$4(ic$2+1,$6,k)-flux$4(ic$2,$6,k))/dx($2)+
     &          (flux$5(ic$3+1,$7,k)-flux$5(ic$3,$7,k))/dx($3))
             enddo

             ttv(1)= trnsvers(1)
             ttv(2)= (trnsvers(2) - vel0*trnsvers(1))/rho
             ttv(3)= (trnsvers(3) - vel1*trnsvers(1))/rho
             ttv(4)= (trnsvers(4) - vel2*trnsvers(1))/rho
             ttv(NEQU)= (v2norm*half*trnsvers(1) -
     &                vel0*trnsvers(2) - vel1*trnsvers(3) -
     &                vel2*trnsvers(4) +
     &                trnsvers(NEQU))*gam_min_one
             do k=1,NEQU
               tracelft$1(ic$1+1,ic$2,ic$3,k)=tracelft$1(ic$1+1,ic$2,ic$3,k)
     &                                           - ttv(k)
               tracergt$1(ic$1  ,ic$2,ic$3,k)=tracergt$1(ic$1  ,ic$2,ic$3,k)
     &                                           - ttv(k)
             enddo
           enddo
         enddo
      enddo
!$OMP END DO
')dnl
define(artificial_viscosity1,`dnl
!$OMP DO SCHEDULE(DYNAMIC)
      do ic$3=ifirst$3-(FLUXG-1),ilast$3+(FLUXG-1)
         do ic$2=ifirst$2-(FLUXG-1),ilast$2+(FLUXG-1)
           do ie$1=ifirst$1-(FLUXG-1),ilast$1+(FLUXG)
             maxeig =trrgt$1(ie$1,ic$2,ic$3,NEQU)-trlft$1(ie$1,ic$2,ic$3,NEQU)
             vcoef = tenth*abs(maxeig)

             mom0L=trlft$1(ie$1,ic$2,ic$3,1)*trlft$1(ie$1,ic$2,ic$3,2)
             mom1L=trlft$1(ie$1,ic$2,ic$3,1)*trlft$1(ie$1,ic$2,ic$3,3)
             mom2L=trlft$1(ie$1,ic$2,ic$3,1)*trlft$1(ie$1,ic$2,ic$3,4)
             v2norm = mom1L**2+mom0L**2+mom2L**2
             enerL=trlft$1(ie$1,ic$2,ic$3,NEQU)/gam_min_one+
     &                v2norm/(trlft$1(ie$1,ic$2,ic$3,1)*two)
             mom0R=trrgt$1(ie$1,ic$2,ic$3,1)*trrgt$1(ie$1,ic$2,ic$3,2)
             mom1R=trrgt$1(ie$1,ic$2,ic$3,1)*trrgt$1(ie$1,ic$2,ic$3,3)
             mom2R=trrgt$1(ie$1,ic$2,ic$3,1)*trrgt$1(ie$1,ic$2,ic$3,4)
             v2norm = mom1R**2+mom0R**2+mom2R**2
             enerR=trrgt$1(ie$1,ic$2,ic$3,NEQU)/gam_min_one+
     &             v2norm/(trrgt$1(ie$1,ic$2,ic$3,1)*two)

             vcorr(1) = dt*vcoef*
     &         (trrgt$1(ie$1,ic$2,ic$3,1)-trlft$1(ie$1,ic$2,ic$3,1))
             vcorr(2) = dt*vcoef*(mom0R-mom0L)
             vcorr(3) = dt*vcoef*(mom1R-mom1L)
             vcorr(4) = dt*vcoef*(mom2R-mom2L)
             vcorr(NEQU) = dt*vcoef*(enerR-enerL)
             do j=1,NEQU
                flux$1(ie$1,ic$2,ic$3,j)=flux$1(ie$1,ic$2,ic$3,j)
     &                                   -vcorr(j)
             enddo
           enddo
         enddo
      enddo
!$OMP END DO
')dnl
c
define(artificial_viscosity2,`dnl
!$OMP DO SCHEDULE(DYNAMIC)
      do ic1=ifirst1,ilast1
        do ie0=ifirst0,ilast0+1
          maxeig =pressure(ie0,ic1)-pressure(ie0-1,ic1)
          vcoef = tenth*abs(maxeig)
          mom0L=density(ie0-1,ic1,1)*velocity(ie0-1,ic1,0)
          mom1L=density(ie0-1,ic1,1)*velocity(ie0-1,ic1,1)
          mom0R=density(ie0  ,ic1,1)*velocity(ie0  ,ic1,0)
          mom1R=density(ie0  ,ic1,1)*velocity(ie0  ,ic1,1)
          v2norm = mom1L**2+mom0L**2
          enerL=pressure(ie0-1,ic1)/gam_min_one+
     &             v2norm/(density(ie0-1,ic1)*two)
          v2norm = mom1R**2+mom0R**2
          enerR=pressure(ie0,ic1)/gam_min_one+
     &             v2norm/(density(ie0,ic1)*two)
          vcorr(1)= dt*vcoef*(density(ie0,ic1)-density(ie0-1,ic1))
          vcorr(2)= dt*vcoef*(mom0R-mom0L)
          vcorr(3)= dt*vcoef*(mom1R-mom1L)
          vcorr(4)= dt*vcoef*(mom2R-mom2L)
          vcorr(NEQU)= dt*vcoef*(enerR-enerL)
             do j=1,NEQU
                flux0(ie0,ic1,ic2,j)=flux0(ie0,ic1,ic2,j)
     &                                   -vcorr(j)
             enddo
        enddo
      enddo
!$OMP END DO

!$OMP DO SCHEDULE(DYNAMIC)
      do ic0=ifirst0,ilast0
        do ie1=ifirst1,ilast1+1
          maxeig =pressure(ic0,ie1)-pressure(ic0-1,ie1)
          vcoef = 0.1*abs(maxeig)
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
!$OMP END DO
')dnl
c

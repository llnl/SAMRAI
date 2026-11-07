c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for corner flux computation in 3d.
c
define(st_third,`dnl
      do ic$3=ifirst$3-FLUXG,ilast$3+FLUXG
         do ic$2=ifirst$2-FLUXG,ilast$2+FLUXG
           do ic$1=ifirst$1-(FLUXG-1),ilast$1+(FLUXG-1)
             rho  = density(ic0,ic1,ic2)
             vel0 = velocity(ic0,ic1,ic2,0)
             vel1 = velocity(ic0,ic1,ic2,1)
             vel2 = velocity(ic0,ic1,ic2,2)
             pres = pressure(ic0,ic1,ic2)
             v2norm= vel0**2+vel1**2 +vel2**2
             do k=1,NEQU
               trnsvers(k)=
     &           (flux$1(ic$1+1,$4,k)-flux$1(ic$1,$4,k))/(3*dx($1))
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
             st3(ic0,ic1,ic2,1)=rho -ttv(1)
             st3(ic0,ic1,ic2,2)=vel0 -ttv(2)
             st3(ic0,ic1,ic2,3)=vel1 -ttv(3)
             st3(ic0,ic1,ic2,4)=vel2 -ttv(4)
             st3(ic0,ic1,ic2,5)=pres -ttv(5)
           enddo
         enddo
      enddo
')dnl
define(f_third,`dnl
c     write(6,*) "checking onedr sol in riemann solve  "
c     write(6,*) "         dt= ",dt
      do ic$3=ifirst$3-FLUXG,ilast$3+FLUXG
         do ic$2=ifirst$2-(FLUXG-1),ilast$2+(FLUXG-1)
           do ic$1=ifirst$1-(FLUXG-1),ilast$1+FLUXG
 
             stateL(1) = st3($5,1)
             stateL(2) = st3($5,2+$1)
             stateL(3) = st3($5,NEQU)
             stateL(4) = stateL(3)/gam_min_one

             stateR(1) = st3(ic0,ic1,ic2,1)
             stateR(2) = st3(ic0,ic1,ic2,2+$1)
             stateR(3) = st3(ic0,ic1,ic2,NEQU)
             stateR(4) = stateR(3)/gam_min_one

             if (rpchoice.eq.APPROX_RIEM_SOLVE) then
                call gas1dapproxrp3d(gamma,stateL,stateR,riemst)
             else if (rpchoice.eq.EXACT_RIEM_SOLVE) then
                call gas1dexactrp3d(gamma,smallr,stateL,stateR,riemst)
c            else if (rpchoice.eq.ROE_RIEM_SOLVE) then
c            else if (rpchoice.eq.HLLC_RIEM_SOLVE) then
             endif
 
             if (riemst(2).le.zero) then
               vel($2)=st3(ic0,ic1,ic2,2+$2)
               vel($3)=st3(ic0,ic1,ic2,2+$3)
             else
               vel($2)=st3($5,2+$2)
               vel($3)=st3($5,2+$3)
             endif
             vel($1)=riemst(2)

             mom$1 =riemst(1)*vel($1)
             v2norm = vel(0)**2+vel(1)**2+vel(2)**2
             Hent =riemst(3)/gam_min_one+v2norm*riemst(1)/two
             flux$1(ic$1,$4,1)= dt*mom$1
             flux$1(ic$1,$4,2+$1)= dt*(mom$1*vel($1)+riemst(3))
             flux$1(ic$1,$4,2+$2)= dt*mom$1*vel($2)
             flux$1(ic$1,$4,2+$3)= dt*mom$1*vel($3)
             flux$1(ic$1,$4,5)= dt*riemst(2)*(Hent+riemst(3))

           enddo
         enddo
      enddo
')dnl
define(correc_fluxjt,`dnl
c   correct the $2-direction with $3-fluxes
      do ic$3=ifirst$3-(FLUXG-1),ilast$3+(FLUXG-1)
        do ic$1=ifirst$1-(FLUXG-1),ilast$1+(FLUXG-1)
           ic$2=ifirst$2-FLUXG
           rho  = density(ic0,ic1,ic2)
           vel0 = velocity(ic0,ic1,ic2,0)
           vel1 = velocity(ic0,ic1,ic2,1)
           vel2 = velocity(ic0,ic1,ic2,2)
           v2norm= vel0**2+vel1**2 +vel2**2
           do k=1,NEQU
             trnsvers(k)=half*
     &        (flux$3(ic$3+1,$4,k)-flux$3(ic$3,$4,k))/dx($3)
           enddo
  
           ttv(1)= trnsvers(1)
           ttv(2)= (trnsvers(2) - vel0*trnsvers(1))/rho
           ttv(3)= (trnsvers(3) - vel1*trnsvers(1))/rho
           ttv(4)= (trnsvers(4) - vel2*trnsvers(1))/rho
           ttv(NEQU)= (v2norm*half*trnsvers(1) -
     &              vel0*trnsvers(2) - vel1*trnsvers(3) -
     &              vel2*trnsvers(4) +
     &              trnsvers(NEQU))*gam_min_one
           do k=1,NEQU
             tracelft$2(ic$2+1,ic$1,ic$3,k)=tracelft$2(ic$2+1,ic$1,ic$3,k)
     &                                         - ttv(k)
           enddo
           do ic$2=ifirst$2-(FLUXG-1),ilast$2+(FLUXG-1)
             rho  = density(ic0,ic1,ic2)
             vel0 = velocity(ic0,ic1,ic2,0)
             vel1 = velocity(ic0,ic1,ic2,1)
             vel2 = velocity(ic0,ic1,ic2,2)
             v2norm= vel0**2+vel1**2 +vel2**2
             do k=1,NEQU
               trnsvers(k)=half*
     &          (flux$3(ic$3+1,$4,k)-flux$3(ic$3,$4,k))/dx($3)
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
               tracelft$2(ic$2+1,ic$1,ic$3,k)=tracelft$2(ic$2+1,ic$1,ic$3,k)
     &                                           - ttv(k)
               tracergt$2(ic$2  ,ic$1,ic$3,k)=tracergt$2(ic$2  ,ic$1,ic$3,k)
     &                                           - ttv(k)
             enddo
           enddo
           ic$2=ilast$2+FLUXG
           rho  = density(ic0,ic1,ic2)
           vel0 = velocity(ic0,ic1,ic2,0)
           vel1 = velocity(ic0,ic1,ic2,1)
           vel2 = velocity(ic0,ic1,ic2,2)
           v2norm= vel0**2+vel1**2 +vel2**2
           do k=1,NEQU
             trnsvers(k)=half*
     &        (flux$3(ic$3+1,$4,k)-flux$3(ic$3,$4,k))/dx($3)
           enddo
  
           ttv(1)= trnsvers(1)
           ttv(2)= (trnsvers(2) - vel0*trnsvers(1))/rho
           ttv(3)= (trnsvers(3) - vel1*trnsvers(1))/rho
           ttv(4)= (trnsvers(4) - vel2*trnsvers(1))/rho
           ttv(NEQU)= (v2norm*half*trnsvers(1) -
     &              vel0*trnsvers(2) - vel1*trnsvers(3) -
     &              vel2*trnsvers(4) +
     &              trnsvers(NEQU))*gam_min_one
           do k=1,NEQU
             tracergt$2(ic$2  ,ic$1,ic$3,k)=tracergt$2(ic$2  ,ic$1,ic$3,k)
     &                                        - ttv(k)
           enddo
         enddo
      enddo
c
c   correct the $3-direction with $2-fluxes
      do ic$1=ifirst$1-(FLUXG-1),ilast$1+(FLUXG-1)
        do ic$2=ifirst$2-(FLUXG-1),ilast$2+(FLUXG-1)
           ic$3=ifirst$3-FLUXG
           rho  = density(ic0,ic1,ic2)
           vel0 = velocity(ic0,ic1,ic2,0)
           vel1 = velocity(ic0,ic1,ic2,1)
           vel2 = velocity(ic0,ic1,ic2,2)
           v2norm= vel0**2+vel1**2 +vel2**2
           do k=1,NEQU
             trnsvers(k)=half*
     &        (flux$2(ic$2+1,$5,k)-flux$2(ic$2,$5,k))/dx($2)
           enddo
  
           ttv(1)= trnsvers(1)
           ttv(2)= (trnsvers(2) - vel0*trnsvers(1))/rho
           ttv(3)= (trnsvers(3) - vel1*trnsvers(1))/rho
           ttv(4)= (trnsvers(4) - vel2*trnsvers(1))/rho
           ttv(NEQU)= (v2norm*half*trnsvers(1) -
     &              vel0*trnsvers(2) - vel1*trnsvers(3) -
     &              vel2*trnsvers(4) +
     &              trnsvers(NEQU))*gam_min_one
           do k=1,NEQU
             tracelft$3(ic$3+1,ic$2,ic$1,k)=tracelft$3(ic$3+1,ic$2,ic$1,k)
     &                                         - ttv(k)
           enddo
           do ic$3=ifirst$3-(FLUXG-1),ilast$3+(FLUXG-1)
             rho  = density(ic0,ic1,ic2)
             vel0 = velocity(ic0,ic1,ic2,0)
             vel1 = velocity(ic0,ic1,ic2,1)
             vel2 = velocity(ic0,ic1,ic2,2)
             v2norm= vel0**2+vel1**2 +vel2**2
             do k=1,NEQU
               trnsvers(k)=half*
     &          (flux$2(ic$2+1,$5,k)-flux$2(ic$2,$5,k))/dx($2)
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
               tracelft$3(ic$3+1,ic$2,ic$1,k)=tracelft$3(ic$3+1,ic$2,ic$1,k)
     &                                           - ttv(k)
               tracergt$3(ic$3  ,ic$2,ic$1,k)=tracergt$3(ic$3  ,ic$2,ic$1,k)
     &                                           - ttv(k)
             enddo
           enddo
           ic$3=ilast$3+FLUXG
           rho  = density(ic0,ic1,ic2)
           vel0 = velocity(ic0,ic1,ic2,0)
           vel1 = velocity(ic0,ic1,ic2,1)
           vel2 = velocity(ic0,ic1,ic2,2)
           v2norm= vel0**2+vel1**2 +vel2**2
           do k=1,NEQU
             trnsvers(k)=half*
     &        (flux$2(ic$2+1,$5,k)-flux$2(ic$2,$5,k))/dx($2)
           enddo
  
           ttv(1)= trnsvers(1)
           ttv(2)= (trnsvers(2) - vel0*trnsvers(1))/rho
           ttv(3)= (trnsvers(3) - vel1*trnsvers(1))/rho
           ttv(4)= (trnsvers(4) - vel2*trnsvers(1))/rho
           ttv(NEQU)= (v2norm*half*trnsvers(1) -
     &              vel0*trnsvers(2) - vel1*trnsvers(3) -
     &              vel2*trnsvers(4) +
     &              trnsvers(NEQU))*gam_min_one
           do k=1,NEQU
             tracergt$3(ic$3  ,ic$2,ic$1,k)=tracergt$3(ic$3  ,ic$2,ic$1,k)
     &                                        - ttv(k)
           enddo
         enddo
      enddo
')dnl

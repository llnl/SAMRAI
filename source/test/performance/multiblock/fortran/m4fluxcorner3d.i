c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for 3d corner flux calculation.
c
define(st_third,`dnl
      do ic$3=ifirst$3-1,ilast$3+1
         do ic$2=ifirst$2-1,ilast$2+1
           do ic$1=ifirst$1,ilast$1
             trnsvers=
     &           (flux$1(ic$1+1,$4)-flux$1(ic$1,$4))/(3*dx($1))
c
             st3(ic0,ic1,ic2)=uval(ic0,ic1,ic2) -trnsvers
           enddo
         enddo
      enddo
')dnl
define(f_third,`dnl
c     write(6,*) "checking onedr sol in riemann solve  "
c     write(6,*) "         dt= ",dt
      do ic$3=ifirst$3-1,ilast$3+1
         do ic$2=ifirst$2,ilast$2
           do ic$1=ifirst$1,ilast$1+1
 
             if (advecspeed($1).ge.zero) then
               riemst = st3($5)
             else
               riemst = st3(ic0,ic1,ic2)
             endif
             flux$1(ic$1,$4)= dt*riemst*advecspeed($1)
           enddo
         enddo
      enddo
')dnl
define(correc_fluxjt,`dnl
c   correct the $2-direction with $3-fluxes
      do ic$3=ifirst$3,ilast$3
        do ic$1=ifirst$1,ilast$1
           ic$2=ifirst$2-1
           trnsvers=half*
     &        (flux$3(ic$3+1,$4)-flux$3(ic$3,$4))/dx($3)
  
           tracelft$2(ic$2+1,ic$1,ic$3)=tracelft$2(ic$2+1,ic$1,ic$3)
     &                                         - trnsvers
           do ic$2=ifirst$2,ilast$2
             trnsvers=half*
     &          (flux$3(ic$3+1,$4)-flux$3(ic$3,$4))/dx($3)
  
             tracelft$2(ic$2+1,ic$1,ic$3)=tracelft$2(ic$2+1,ic$1,ic$3)
     &                                           - trnsvers
             tracergt$2(ic$2  ,ic$1,ic$3)=tracergt$2(ic$2  ,ic$1,ic$3)
     &                                           - trnsvers
           enddo
           ic$2=ilast$2+1
           trnsvers=half*
     &        (flux$3(ic$3+1,$4)-flux$3(ic$3,$4))/dx($3)
  
           tracergt$2(ic$2  ,ic$1,ic$3)=tracergt$2(ic$2  ,ic$1,ic$3)
     &                                        - trnsvers
         enddo
      enddo
c
c   correct the $3-direction with $2-fluxes
      do ic$1=ifirst$1,ilast$1
        do ic$2=ifirst$2,ilast$2
           ic$3=ifirst$3-1
           trnsvers=half*
     &        (flux$2(ic$2+1,$5)-flux$2(ic$2,$5))/dx($2)
  
           tracelft$3(ic$3+1,ic$2,ic$1)=tracelft$3(ic$3+1,ic$2,ic$1)
     &                                         - trnsvers
           do ic$3=ifirst$3,ilast$3
             trnsvers=half*
     &          (flux$2(ic$2+1,$5)-flux$2(ic$2,$5))/dx($2)
  
             tracelft$3(ic$3+1,ic$2,ic$1)=tracelft$3(ic$3+1,ic$2,ic$1)
     &                                           - trnsvers
             tracergt$3(ic$3  ,ic$2,ic$1)=tracergt$3(ic$3  ,ic$2,ic$1)
     &                                           - trnsvers
           enddo
           ic$3=ilast$3+1
           trnsvers=half*
     &        (flux$2(ic$2+1,$5)-flux$2(ic$2,$5))/dx($2)
  
           tracergt$3(ic$3  ,ic$2,ic$1)=tracergt$3(ic$3  ,ic$2,ic$1)
     &                                        - trnsvers
         enddo
      enddo
')dnl

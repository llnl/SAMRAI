c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for 2d trace calculation.
c
define(trace_init,`dnl
         do  ic$2=ifirst$2-FACEG,ilast$2+FACEG
           ie$1=ifirst$1-FACEG
           tracelft$1(ie$1,ic$2,1)=zero
           tracelft$1(ie$1,ic$2,2)=zero
           tracelft$1(ie$1,ic$2,3)=zero
           tracelft$1(ie$1,ic$2,NEQU)=zero
           tracergt$1(ie$1,ic$2,1)=density($5)
           tracergt$1(ie$1,ic$2,2)=velocity($5,1)
           tracergt$1(ie$1,ic$2,3)=velocity($5,2)
           tracergt$1(ie$1,ic$2,NEQU)=pressure($5)

           do  ie$1=ifirst$1+1-FACEG,ilast$1+FACEG
             tracelft$1(ie$1,ic$2,1)=density($4)
             tracelft$1(ie$1,ic$2,2)=velocity($4,1)
             tracelft$1(ie$1,ic$2,3)=velocity($4,2)
             tracelft$1(ie$1,ic$2,NEQU)=pressure($4)

             tracergt$1(ie$1,ic$2,1)=density($5)
             tracergt$1(ie$1,ic$2,2)=velocity($5,1)
             tracergt$1(ie$1,ic$2,3)=velocity($5,2)
             tracergt$1(ie$1,ic$2,NEQU)=pressure($5)

           enddo

           ie$1=ilast$1+FACEG+1
           tracelft$1(ie$1,ic$2,1)=density($4)
           tracelft$1(ie$1,ic$2,2)=velocity($4,1)
           tracelft$1(ie$1,ic$2,3)=velocity($4,2)
           tracelft$1(ie$1,ic$2,NEQU)=pressure($4)
           tracergt$1(ie$1,ic$2,1)=zero
           tracergt$1(ie$1,ic$2,2)=zero
           tracergt$1(ie$1,ic$2,3)=zero
           tracergt$1(ie$1,ic$2,NEQU)=zero
         enddo
 
')dnl
define(trace_call,`dnl
          do ic$2=ifirst$2-2,ilast$2+2
            do ic$1=ifirst$1-CELLG,ilast$1+CELLG
              ttsound(ic$1)= sound(ic0,ic1)
            enddo
            do ic$1=ifirst$1-FACEG,ilast$1+FACEG+1
              do k=1,NEQU
                ttraclft(ic$1,k) = tracelft(ic$1,ic$2,k)
                ttracrgt(ic$1,k) = tracergt(ic$1,ic$2,k)
              enddo
            enddo

            call trace2d(dt,ifirst$1,ilast$1,mc,
     &        dx,idir,igdnv,ttsound,
     &        ttraclft,ttracrgt,
     &        ttcelslp,ttedgslp)
            do ic$1=ifirst$1-FACEG,ilast$1+FACEG+1
              do k=1,NEQU
                tracelft(ic$1,ic$2,k) = ttraclft(ic$1,k)
                tracergt(ic$1,ic$2,k) = ttracrgt(ic$1,k)
              enddo
            enddo
          enddo
')dnl

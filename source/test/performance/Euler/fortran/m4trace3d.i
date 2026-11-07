c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for trace computation in 3d.
c
define(trace_init,`dnl
      do  ic$3=ifirst$3-FACEG,ilast$3+FACEG
         do  ic$2=ifirst$2-FACEG,ilast$2+FACEG
           ie$1=ifirst$1-FACEG
             tracelft$1(ie$1,ic$2,ic$3,1)=zero
             tracelft$1(ie$1,ic$2,ic$3,2)=zero
             tracelft$1(ie$1,ic$2,ic$3,3)=zero
             tracelft$1(ie$1,ic$2,ic$3,4)=zero
             tracelft$1(ie$1,ic$2,ic$3,5)=zero
             tracergt$1(ie$1,ic$2,ic$3,1)=density($5)
             tracergt$1(ie$1,ic$2,ic$3,2)=velocity($5,1)
             tracergt$1(ie$1,ic$2,ic$3,3)=velocity($5,2)
             tracergt$1(ie$1,ic$2,ic$3,4)=velocity($5,3)
             tracergt$1(ie$1,ic$2,ic$3,5)=pressure($5)

           do  ie$1=ifirst$1+1-FACEG,ilast$1+FACEG
             tracelft$1(ie$1,ic$2,ic$3,1)=density($4)
             tracelft$1(ie$1,ic$2,ic$3,2)=velocity($4,1)
             tracelft$1(ie$1,ic$2,ic$3,3)=velocity($4,2)
             tracelft$1(ie$1,ic$2,ic$3,4)=velocity($4,3)
             tracelft$1(ie$1,ic$2,ic$3,5)=pressure($4)
 
             tracergt$1(ie$1,ic$2,ic$3,1)=density($5)
             tracergt$1(ie$1,ic$2,ic$3,2)=velocity($5,1)
             tracergt$1(ie$1,ic$2,ic$3,3)=velocity($5,2)
             tracergt$1(ie$1,ic$2,ic$3,4)=velocity($5,3)
             tracergt$1(ie$1,ic$2,ic$3,5)=pressure($5)
    
           enddo

           ie$1=ilast$1+FACEG+1
             tracelft$1(ie$1,ic$2,ic$3,1)=density($4)
             tracelft$1(ie$1,ic$2,ic$3,2)=velocity($4,1)
             tracelft$1(ie$1,ic$2,ic$3,3)=velocity($4,2)
             tracelft$1(ie$1,ic$2,ic$3,4)=velocity($4,3)
             tracelft$1(ie$1,ic$2,ic$3,5)=pressure($4)
             tracergt$1(ie$1,ic$2,ic$3,1)=zero
             tracergt$1(ie$1,ic$2,ic$3,2)=zero
             tracergt$1(ie$1,ic$2,ic$3,3)=zero
             tracergt$1(ie$1,ic$2,ic$3,4)=zero
             tracergt$1(ie$1,ic$2,ic$3,5)=zero
         enddo
      enddo

')dnl
define(trace_call,`dnl
        do ic$3=ifirst$3-2,ilast$3+2
          do ic$2=ifirst$2-2,ilast$2+2
            do ic$1=ifirst$1-CELLG,ilast$1+CELLG
              ttsound(ic$1)= sound(ic0,ic1,ic2)
            enddo
            do k=1,NEQU
              do ic$1=ifirst$1-FACEG,ilast$1+FACEG+1
                ttraclft(ic$1,k) = tracelft(ic$1,ic$2,ic$3,k)
                ttracrgt(ic$1,k) = tracergt(ic$1,ic$2,ic$3,k)
              enddo
            enddo
   
            call trace3d(dt,ifirst$1,ilast$1,mc,
     &        dx,idir,igdnv,ttsound,
     &        ttraclft,ttracrgt,
     &        ttcelslp,ttedgslp)
            do k=1,NEQU
              do ic$1=ifirst$1-FACEG,ilast$1+FACEG+1
                tracelft(ic$1,ic$2,ic$3,k) = ttraclft(ic$1,k)
                tracergt(ic$1,ic$2,ic$3,k) = ttracrgt(ic$1,k)
              enddo
            enddo
          enddo
        enddo
')dnl

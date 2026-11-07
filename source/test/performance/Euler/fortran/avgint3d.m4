c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for conservative interlevel transfer of velocity 
c                and pressure for 3d euler equations.
c
define(NDIM,3)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl
include(FORTDIR/amrflaten3d.i)dnl
c
define(coarsen_index,`dnl
      if ($1.lt.0) then
         $2=($1+1)/$3-1
      else
         $2=$1/$3
      endif
')dnl
define(coarse_fine_cell_deltas,`dnl
      do ir$1=0,ratio($1)-1
         deltax(ir$1,$1)=(dble(ir$1)+half)*dxf($1)-dxc($1)*half
      enddo
')dnl
define(muscl_limited_conserved_slopes,`dnl
      do ie$1=ifirstc$1,ilastc$1+1
         diff$1(ie$1)=conservc($2)
     &               -conservc($3)
      enddo
      do ic$1=ifirstc$1,ilastc$1
         coef2=half*(diff$1(ic$1+1)+diff$1(ic$1))
         bound=two*min(abs(diff$1(ic$1+1)),abs(diff$1(ic$1)))
         if (diff$1(ic$1)*diff$1(ic$1+1).gt.zero) then
            slope$1($4)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc($1)
         else
            slope$1($4)=zero
         endif
         slope$1($4)=slope$1($4)*flat$1($4)
      enddo
')dnl
c
c***********************************************************************
c Conservative linear interpolation for 3d flux variables.
c***********************************************************************
c
      subroutine conservlinint3d(
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  gamma,
     &  densc,densf,
     &  velc,presc,
     &  velf,presf,
     &  conservc,
     &  tflat,tflat2,sound,mc,
     &  tdensc,tpresc,tvelc,
     &  flat0,flat1,flat2,
     &  diff0,slope0,diff1,slope1,diff2,slope2)
c***********************************************************************
      implicit none
      REAL zero,half,one,two
      parameter (zero=0.0d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:NDIM-1)
      REAL
     &  dxc(0:NDIM-1),
     &  dxf(0:NDIM-1),
     &  gamma
      REAL
     &  densc(CELL3d(cilo,cihi,0)),
     &  densf(CELL3d(filo,fihi,0)),
     &  velc(CELL3d(cilo,cihi,0),0:NDIM-1),
     &  presc(CELL3d(cilo,cihi,0))
      REAL 
     &  velf(CELL3d(filo,fihi,0),0:NDIM-1),
     &  presf(CELL3d(filo,fihi,0))
      REAL
     &  conservc(CELL3d(ifirstc,ilastc,1)),
     &  flat0(CELL3d(ifirstc,ilastc,0)),
     &  flat1(CELL3d(ifirstc,ilastc,0)),
     &  flat2(CELL3d(ifirstc,ilastc,0)),
     &  diff0(ifirstc0:ilastc0+1),
     &  slope0(CELL3d(ifirstc,ilastc,0)),
     &  diff1(ifirstc1:ilastc1+1),
     &  slope1(CELL3d(ifirstc,ilastc,0)),
     &  diff2(ifirstc2:ilastc2+1),
     &  slope2(CELL3d(ifirstc,ilastc,0))
      REAL
     &  deltax(0:15,0:NDIM-1)
      integer mc
      REAL
     &  tdensc(0:mc-1),tpresc(0:mc-1),tvelc(0:mc-1)
      REAL
     &  tflat(0:mc-1),tflat2(0:mc-1),sound(0:mc-1)

      integer id,ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,ir0,ir1,ir2,it
      REAL coef2,bound,val,valinv,
     &  deltax1,deltax2,v2norm
      logical presneg
c
c***********************************************************************
c

      presneg = .false.

coarse_fine_cell_deltas(0)dnl

coarse_fine_cell_deltas(1)dnl

coarse_fine_cell_deltas(2)dnl

c
c     compute the flatten coefficients to further limit the slopes
c
      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ic0=ifirstc0-1,ilastc0+1
               it = ic0-cilo0
               tdensc(it) = densc(ic0,ic1,ic2)
               tpresc(it) = presc(ic0,ic1,ic2)
               tvelc(it)  = velc(ic0,ic1,ic2,0)
            enddo
            call amrflaten3d(ifirstc0,ilastc0,cilo0,cihi0,mc,gamma,
     &                     tdensc,tpresc,tvelc,tflat,tflat2,sound)
            do ic0=ifirstc0,ilastc0
                 flat0(ic0,ic1,ic2) = tflat(ic0-ifirstc0)
            enddo
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         do ic2=ifirstc2,ilastc2
            do ic1=ifirstc1-1,ilastc1+1
               it = ic1-cilo1
               tdensc(it) = densc(ic0,ic1,ic2)
               tpresc(it) = presc(ic0,ic1,ic2)
               tvelc(it)  = velc(ic0,ic1,ic2,1)
            enddo
            call amrflaten3d(ifirstc1,ilastc1,cilo1,cihi1,mc,gamma,
     &                     tdensc,tpresc,tvelc,tflat,tflat2,sound)
            do ic1=ifirstc1,ilastc1
                 flat1(ic0,ic1,ic2) = tflat(ic1-ifirstc1)
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            do ic2=ifirstc2-1,ilastc2+1
               it = ic2-cilo2
               tdensc(it) = densc(ic0,ic1,ic2)
               tpresc(it) = presc(ic0,ic1,ic2)
               tvelc(it)  = velc(ic0,ic1,ic2,2)
            enddo
            call amrflaten3d(ifirstc2,ilastc2,cilo2,cihi2,mc,gamma,
     &                  tdensc,tpresc,tvelc,tflat,tflat2,sound)
            do ic2=ifirstc2,ilastc2
                 flat2(ic0,ic1,ic2) = tflat(ic2-ifirstc2)
            enddo
         enddo
      enddo

c
c

c
c construct fine velocity values using conservative linear 
c interpolation on momentum
c
      do id=0,NDIM-1

         do ic2=ifirstc2-1,ilastc2+1
            do ic1=ifirstc1-1,ilastc1+1
               do ic0=ifirstc0-1,ilastc0+1
                  conservc(ic0,ic1,ic2) = 
     &            densc(ic0,ic1,ic2)*velc(ic0,ic1,ic2,id) 
               enddo
            enddo
         enddo

         do ic2=ifirstc2,ilastc2
            do ic1=ifirstc1,ilastc1
muscl_limited_conserved_slopes(0,`ie0,ic1,ic2',`ie0-1,ic1,ic2',`ic0,ic1,ic2')dnl
            enddo
         enddo

         do ic2=ifirstc2,ilastc2
            do ic0=ifirstc0,ilastc0
muscl_limited_conserved_slopes(1,`ic0,ie1,ic2',`ic0,ie1-1,ic2',`ic0,ic1,ic2')dnl
            enddo
         enddo

         do ic1=ifirstc1,ilastc1
            do ic0=ifirstc0,ilastc0
muscl_limited_conserved_slopes(2,`ic0,ic1,ie2',`ic0,ic1,ie2-1',`ic0,ic1,ic2')dnl
            enddo
         enddo

         do if2=ifirstf2,ilastf2
coarsen_index(if2,ic2,ratio(2))dnl
            ir2=if2-ic2*ratio(2)
            deltax2=deltax(ir2,2)
            do if1=ifirstf1,ilastf1
coarsen_index(if1,ic1,ratio(1))dnl
               ir1=if1-ic1*ratio(1)
               deltax1=deltax(ir1,1)
               do if0=ifirstf0,ilastf0
coarsen_index(if0,ic0,ratio(0))dnl
                  ir0=if0-ic0*ratio(0)
                  velf(if0,if1,if2,id)=(conservc(ic0,ic1,ic2)
     &                        + slope0(ic0,ic1,ic2)*deltax(ir0,0)
     &                        + slope1(ic0,ic1,ic2)*deltax1
     &                        + slope2(ic0,ic1,ic2)*deltax2)
     &                        / densf(if0,if1,if2)
             enddo
           enddo
         enddo

      enddo

c
c construct fine pressure values using conservative linear 
c interpolation on energy
c
      val = (gamma-one) 
      valinv = one/(gamma-one) 

      do ic2=ifirstc2-1,ilastc2+1
         do ic1=ifirstc1-1,ilastc1+1
            do ic0=ifirstc0-1,ilastc0+1
               v2norm = velc(ic0,ic1,ic2,0)**2
     &                + velc(ic0,ic1,ic2,1)**2
     &                + velc(ic0,ic1,ic2,2)**2
               conservc(ic0,ic1,ic2) = presc(ic0,ic1,ic2)*valinv
     &                               + half*densc(ic0,ic1,ic2)*v2norm
            enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
muscl_limited_conserved_slopes(0,`ie0,ic1,ic2',`ie0-1,ic1,ic2',`ic0,ic1,ic2')dnl
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic0=ifirstc0,ilastc0
muscl_limited_conserved_slopes(1,`ic0,ie1,ic2',`ic0,ie1-1,ic2',`ic0,ic1,ic2')dnl
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
muscl_limited_conserved_slopes(2,`ic0,ic1,ie2',`ic0,ic1,ie2-1',`ic0,ic1,ic2')dnl
         enddo
      enddo

      do if2=ifirstf2,ilastf2
coarsen_index(if2,ic2,ratio(2))dnl
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)
         do if1=ifirstf1,ilastf1
coarsen_index(if1,ic1,ratio(1))dnl
            ir1=if1-ic1*ratio(1)
            deltax1=deltax(ir1,1)
            do if0=ifirstf0,ilastf0
coarsen_index(if0,ic0,ratio(0))dnl
               ir0=if0-ic0*ratio(0)
               v2norm = velf(if0,if1,if2,0)**2
     &                + velf(if0,if1,if2,1)**2
     &                + velf(if0,if1,if2,2)**2
               presf(if0,if1,if2) = ((conservc(ic0,ic1,ic2)
     &                            + slope0(ic0,ic1,ic2)*deltax(ir0,0)
     &                            + slope1(ic0,ic1,ic2)*deltax1
     &                            + slope2(ic0,ic1,ic2)*deltax2)
     &                     - half*densf(if0,if1,if2)*v2norm) * val
            if (presf(if0,if1,if2).lt.zero) then
c             write(6,*) "IN conservlinint3d"
c             write(6,*) "gamma,val,valinv = ",
c    &                      gamma,val,valinv
c               write(6,*) "ifirstf0,ilastf0 ",ifirstf0,ilastf0,
c    &                     ", ifirstf1,ilastf1 ",ifirstf1,ilastf1
c    &                     ", ifirstf2,ilastf2 ",ifirstf2,ilastf2
c               write(6,*) "if0,if1,if2 ",if0,if,if21
c               write(6,*) "ic0,ic1,if2 ",ic0,ic,if21
c               write(6,*) "fine energy = ",
c    &   conservc(ic0,ic1,ic2) + slope0(ic0,ic1,ic2)*deltax(ir0,0)
c    &                     +  slope1(ic0,ic1,ic2)*deltax1
c    &                     +  slope2(ic0,ic1,ic2)*deltax2
c               write(6,*) "densf(if0,if1,if2) ",densf(if0,if1,if2)
c               write(6,*) "vel() ",velf(if0,if1,if2,0),
c    &                              velf(if0,if1,if2,1),
c    &                              velf(if0,if1,if2,2)
c               write(6,*) "presf(if0,if1,if2) ",presf(if0,if1,if2)
c               write(6,*) "conservc(",ic0,ic1,ic2,") ",
c    &                      conservc(ic0,ic1,ic2)
c               write(6,*) "slope0(",ic0,ic1,ic2,") ",
c    &                      slope0(ic0,ic1,ic2)
c               write(6,*) "deltax(",ir0,0,") ",deltax(ir0,0)
c               write(6,*) "slope1(",ic0,ic1,ic2,") ",
c    &                      slope1(ic0,ic1,ic2)
c               write(6,*) "deltax1(",ir1,1,") ",deltax1
c               write(6,*) "slope2(",ic0,ic1,ic2,") ",
c    &                      slope2(ic0,ic1,ic2)
c               write(6,*) "deltax2(",ir2,2,") ",deltax2
c               write(6,*) "val ",val,", v2norm ",v2norm
c               write(6,*) "presc(",ic0,ic1,ic2,") ",presc(ic0,ic1,ic2)
c               write(6,*) "densc(",ic0,ic1,ic2,") ",densc(ic0,ic1,ic2)
c               write(6,*) "velc() ",velc(ic0,ic1,ic2,0),
c    &                               velc(ic0,ic1,ic2,1),
c    &                               velc(ic0,ic1,ic2,2)
                presneg = .true.
            endif
          enddo
        enddo
      enddo

      if (presneg) then
        write(6,*) "negative pressure reported in conservlinint3d"
c       write(6,*) "coarse conserved values"
c       do ic2=ifirstc2-1,ilastc2+1
c          do ic1=ifirstc1-1,ilastc1+1
c             do ic0=ifirstc0-1,ilastc0+1
c                write(6,*) "ic0,ic1,ic2,conservc = ",
c    &                       ic0,ic1,ic2,conservc(ic0,ic1,ic2)
c             enddo
c          enddo
c       enddo
c       write(6,*) "coarse density values"
c       do ic2=ifirstc2-1,ilastc2+1
c          do ic1=ifirstc1-1,ilastc1+1
c             do ic0=ifirstc0-1,ilastc0+1
c                write(6,*) "ic0,ic1,ic2,densc = ",
c    &                       ic0,ic1,ic2,densc(ic0,ic1,ic2)
c             enddo
c          enddo
c       enddo
c       write(6,*) "coarse velocity0 values"
c       do ic2=ifirstc2-1,ilastc2+1
c          do ic1=ifirstc1-1,ilastc1+1
c             do ic0=ifirstc0-1,ilastc0+1
c                write(6,*) "ic0,ic1,ic2,velc = ",
c    &                       ic0,ic1,ic2,velc(ic0,ic1,ic2,0)
c             enddo
c          enddo
c       enddo
c       write(6,*) "coarse velocity1 values"
c       do ic2=ifirstc2-1,ilastc2+1
c          do ic1=ifirstc1-1,ilastc1+1
c             do ic0=ifirstc0-1,ilastc0+1
c                write(6,*) "ic0,ic1,ic2,velc = ",
c    &                       ic0,ic1,ic2,velc(ic0,ic1,ic2,1)
c             enddo
c          enddo
c       enddo
c       write(6,*) "coarse velocity2 values"
c       do ic2=ifirstc2-1,ilastc2+1
c          do ic1=ifirstc1-1,ilastc1+1
c             do ic0=ifirstc0-1,ilastc0+1
c                write(6,*) "ic0,ic1,ic2,velc = ",
c    &                       ic0,ic1,ic2,velc(ic0,ic1,ic2,2)
c             enddo
c          enddo
c       enddo
c       write(6,*) "coarse pressure values"
c       do ic2=ifirstc2-1,ilastc2+1
c          do ic1=ifirstc1-1,ilastc1+1
c             do ic0=ifirstc0-1,ilastc0+1
c                write(6,*) "ic0,ic1,ic2,presc = ",
c    &                       ic0,ic1,ic2,presc(ic0,ic1,ic2)
c             enddo
c          enddo
c       enddo
c       write(6,*) "fine density values"
c       do if2=ifirstf2,ilastf2
c          do if1=ifirstf1,ilastf1
c             do if0=ifirstf0,ilastf0
c                write(6,*) "if0,if1,if2,densf = ",
c    &           if0,if1,if2,densf(if0,if1,if2)
c             enddo
c          enddo
c       enddo
c       write(6,*) "fine velocity0 values"
c       do if2=ifirstf2,ilastf2
c          do if1=ifirstf1,ilastf1
c             do if0=ifirstf0,ilastf0
c                write(6,*) "if0,if1,if2,velf = ",
c    &           if0,if1,if2,velf(if0,if1,if2,0)
c             enddo
c          enddo
c       enddo
c       write(6,*) "fine velocity1 values"
c       do if2=ifirstf2,ilastf2
c          do if1=ifirstf1,ilastf1
c             do if0=ifirstf0,ilastf0
c                write(6,*) "if0,if1,if2,velf = ",
c    &           if0,if1,if2,velf(if0,if1,if2,1)
c             enddo
c          enddo
c       enddo
c       write(6,*) "fine velocity2 values"
c       do if2=ifirstf2,ilastf2
c          do if1=ifirstf1,ilastf1
c             do if0=ifirstf0,ilastf0
c                write(6,*) "if0,if1,if2,velf = ",
c    &           if0,if1,if2,velf(if0,if1,if2,2)
c             enddo
c          enddo
c       enddo
        stop
      endif
c
      return
      end
c
c
c***********************************************************************
c Volume weighted averaging for 3d flux variables.
c***********************************************************************
c
      subroutine conservavg3d(
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
     &  gamma,
     &  densf,densc,
     &  velf,presf,
     &  velc,presc,
     &  conservf)
c***********************************************************************
      implicit none
      REAL zero,half,one 
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
c
      integer
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2
      integer ratio(0:NDIM-1)
      REAL
     &  dxf(0:NDIM-1),
     &  dxc(0:NDIM-1),
     &  gamma
      REAL
     &  densf(CELL3d(filo,fihi,0)),
     &  densc(CELL3d(cilo,cihi,0)),
     &  velf(CELL3d(filo,fihi,0),0:NDIM-1),
     &  presf(CELL3d(filo,fihi,0))
      REAL
     &  velc(CELL3d(cilo,cihi,0),0:NDIM-1),
     &  presc(CELL3d(cilo,cihi,0))
      REAL
     &  conservf(CELL3d(ifirstf,ilastf,0))
      integer ic0,ic1,ic2,if0,if1,if2,ir0,ir1,ir2,id
      REAL dVf,dVc,val,valinv,v2norm
c
c***********************************************************************
c
      dVf = dxf(0)*dxf(1)*dxf(2)
      dVc = dxc(0)*dxc(1)*dxc(2)
c
c construct coarse velocity values using conservative average
c of momentum
c
      do id=0,NDIM-1

         do if2=ifirstf2,ilastf2
            do if1=ifirstf1,ilastf1
               do if0=ifirstf0,ilastf0
                  conservf(if0,if1,if2) = 
     &            densf(if0,if1,if2)*velf(if0,if1,if2,id)
               enddo
            enddo
         enddo

         do ic2=ifirstc2,ilastc2
            do ic1=ifirstc1,ilastc1
               do ic0=ifirstc0,ilastc0
                  velc(ic0,ic1,ic2,id) = zero
               enddo
            enddo
         enddo

         do ir2=0,ratio(2)-1
            do ir1=0,ratio(1)-1
               do ir0=0,ratio(0)-1
                  do ic2=ifirstc2,ilastc2
                     if2=ic2*ratio(2)+ir2
                     do ic1=ifirstc1,ilastc1
                        if1=ic1*ratio(1)+ir1
                        do ic0=ifirstc0,ilastc0
                           if0=ic0*ratio(0)+ir0
                           velc(ic0,ic1,ic2,id) = velc(ic0,ic1,ic2,id)
     &                     + conservf(if0,if1,if2)*dVf
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo

         do ic2=ifirstc2,ilastc2
            do ic1=ifirstc1,ilastc1
               do ic0=ifirstc0,ilastc0
                  velc(ic0,ic1,ic2,id) = velc(ic0,ic1,ic2,id)
     &                                 / (densc(ic0,ic1,ic2)*dVc)
               enddo
            enddo
         enddo

      enddo

c
c construct coarse velocity values using conservative average
c of energy
c
      val = (gamma-one)
      valinv = one/(gamma-one)

      do if2=ifirstf2,ilastf2
         do if1=ifirstf1,ilastf1
            do if0=ifirstf0,ilastf0
               v2norm = velf(if0,if1,if2,0)**2
     &                + velf(if0,if1,if2,1)**2
     &                + velf(if0,if1,if2,2)**2
               conservf(if0,if1,if2) = presf(if0,if1,if2)*valinv
     &                               + half*densf(if0,if1,if2)*v2norm
            enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ic0=ifirstc0,ilastc0
               presc(ic0,ic1,ic2) = zero
            enddo
         enddo
      enddo

      do ir2=0,ratio(2)-1
         do ir1=0,ratio(1)-1
            do ir0=0,ratio(0)-1
               do ic2=ifirstc2,ilastc2
                  if2=ic2*ratio(2)+ir2
                  do ic1=ifirstc1,ilastc1
                     if1=ic1*ratio(1)+ir1
                     do ic0=ifirstc0,ilastc0
                        if0=ic0*ratio(0)+ir0
                        presc(ic0,ic1,ic2) = presc(ic0,ic1,ic2)
     &                  + conservf(if0,if1,if2)*dVf
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ic0=ifirstc0,ilastc0
               v2norm = velc(ic0,ic1,ic2,0)**2
     &                + velc(ic0,ic1,ic2,1)**2
     &                + velc(ic0,ic1,ic2,2)**2
               presc(ic0,ic1,ic2) = ((presc(ic0,ic1,ic2)/dVc)
     &         - half*densc(ic0,ic1,ic2)*v2norm) * val 
            enddo
         enddo
      enddo


c
      return
      end

c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for interlevel transfer of velocity and pressure
c                for 2d euler equations.
c
define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl
include(FORTDIR/amrflaten2d.i)dnl
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
c Conservative linear interpolation for 2d flux variables.
c***********************************************************************
c
      subroutine conservlinint2d(
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  gamma,
     &  densc,densf,
     &  velc,presc,
     &  velf,presf,
     &  conservc,
     &  tflat,tflat2,sound,mc,
     &  tdensc,tpresc,tvelc,
     &  flat0,flat1,
     &  diff0,slope0,diff1,slope1)
c***********************************************************************
      implicit none
include(FORTDIR/const.i)dnl
c***********************************************************************
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:NDIM-1)
      REAL
     &  dxc(0:NDIM-1),
     &  dxf(0:NDIM-1),
     &  gamma
      REAL
     &  densc(CELL2d(cilo,cihi,0)),
     &  densf(CELL2d(filo,fihi,0)),
     &  velc(CELL2d(cilo,cihi,0),0:NDIM-1),
     &  presc(CELL2d(cilo,cihi,0))
      REAL 
     &  velf(CELL2d(filo,fihi,0),0:NDIM-1),
     &  presf(CELL2d(filo,fihi,0))
      REAL
     &  conservc(CELL2d(ifirstc,ilastc,1)),
     &  flat0(CELL2d(ifirstc,ilastc,0)),
     &  flat1(CELL2d(ifirstc,ilastc,0)),
     &  diff0(ifirstc0:ilastc0+1),
     &  slope0(CELL2d(ifirstc,ilastc,0)),
     &  diff1(ifirstc1:ilastc1+1),
     &  slope1(CELL2d(ifirstc,ilastc,0))
      REAL
     &  deltax(0:15,0:NDIM-1)
      integer mc
      REAL
     &  tdensc(0:mc-1),tpresc(0:mc-1),tvelc(0:mc-1)
      REAL
     &  tflat(0:mc-1),tflat2(0:mc-1),sound(0:mc-1)

      integer id,ic0,ic1,ie0,ie1,if0,if1,ir0,ir1,it
      REAL coef2,bound,val,valinv,deltax1,v2norm
      logical presneg
c
c***********************************************************************
c

      presneg = .false.

coarse_fine_cell_deltas(0)dnl

coarse_fine_cell_deltas(1)dnl

c
c     compute the flatten coefficients to further limit the slopes
c
      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0-1,ilastc0+1
            it = ic0-cilo0
            tdensc(it) = densc(ic0,ic1)
            tpresc(it) = presc(ic0,ic1)
            tvelc(it)  = velc(ic0,ic1,0)
c           write(6,*) "presc(",ic0,ic1,") ",presc(ic0,ic1)
c           write(6,*) "tpresc(",it,") ",tpresc(it)
         enddo
         call amrflaten2d(ifirstc0,ilastc0,cilo0,cihi0,mc,gamma,
     &                  tdensc,tpresc,tvelc,tflat,tflat2,sound)
         do ic0=ifirstc0,ilastc0
              flat0(ic0,ic1) = tflat(ic0-ifirstc0)
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         do ic1=ifirstc1-1,ilastc1+1
            it = ic1-cilo1
            tdensc(it) = densc(ic0,ic1)
            tpresc(it) = presc(ic0,ic1)
            tvelc(it)  = velc(ic0,ic1,1)
         enddo
         call amrflaten2d(ifirstc1,ilastc1,cilo1,cihi1,mc,gamma,
     &                  tdensc,tpresc,tvelc,tflat,tflat2,sound)
         do ic1=ifirstc1,ilastc1
              flat1(ic0,ic1) = tflat(ic1-ifirstc1)
         enddo
      enddo

c
c

c
c construct fine velocity values using conservative linear 
c interpolation on momentum
c
      do id=0,NDIM-1

         do ic1=ifirstc1-1,ilastc1+1
            do ic0=ifirstc0-1,ilastc0+1
               conservc(ic0,ic1) = densc(ic0,ic1)*velc(ic0,ic1,id) 
            enddo
         enddo

         do ic1=ifirstc1,ilastc1
muscl_limited_conserved_slopes(0,`ie0,ic1',`ie0-1,ic1',`ic0,ic1')dnl
         enddo

         do ic0=ifirstc0,ilastc0
muscl_limited_conserved_slopes(1,`ic0,ie1',`ic0,ie1-1',`ic0,ic1')dnl
         enddo

         do if1=ifirstf1,ilastf1
coarsen_index(if1,ic1,ratio(1))dnl
            ir1=if1-ic1*ratio(1)
            deltax1=deltax(ir1,1)
            do if0=ifirstf0,ilastf0
coarsen_index(if0,ic0,ratio(0))dnl
               ir0=if0-ic0*ratio(0)
               velf(if0,if1,id) = (conservc(ic0,ic1)
     &                          + slope0(ic0,ic1)*deltax(ir0,0)
     &                          + slope1(ic0,ic1)*deltax1)
     &                          / densf(if0,if1)
           enddo
         enddo

      enddo

c
c construct fine pressure values using conservative linear 
c interpolation on energy
c
      val = (gamma-one) 
      valinv = one/(gamma-one) 

      do ic1=ifirstc1-1,ilastc1+1
         do ic0=ifirstc0-1,ilastc0+1
            v2norm = velc(ic0,ic1,0)**2+velc(ic0,ic1,1)**2 
            conservc(ic0,ic1) = presc(ic0,ic1)*valinv
     &                        + half*densc(ic0,ic1)*v2norm
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
muscl_limited_conserved_slopes(0,`ie0,ic1',`ie0-1,ic1',`ic0,ic1')dnl
      enddo

      do ic0=ifirstc0,ilastc0
muscl_limited_conserved_slopes(1,`ic0,ie1',`ic0,ie1-1',`ic0,ic1')dnl
      enddo

      do if1=ifirstf1,ilastf1
coarsen_index(if1,ic1,ratio(1))dnl
         ir1=if1-ic1*ratio(1)
         deltax1=deltax(ir1,1)
         do if0=ifirstf0,ilastf0
coarsen_index(if0,ic0,ratio(0))dnl
            ir0=if0-ic0*ratio(0)
            v2norm = velf(if0,if1,0)**2+velf(if0,if1,1)**2
            presf(if0,if1) = ((conservc(ic0,ic1)
     &                     +  slope0(ic0,ic1)*deltax(ir0,0)
     &                     +  slope1(ic0,ic1)*deltax1)
     &                     - half*densf(if0,if1)*v2norm) * val
            presf(if0,if1) = max(presf(if0,if1),smallr)
            if (presf(if0,if1).lt.zero) then
c               write(6,*) "IN conservlinint2d"
c               write(6,*) "gamma,val,valinv = ",
c    &                      gamma,val,valinv
c               write(6,*) "ifirstf0,ilastf0 ",ifirstf0,ilastf0,
c    &                     ", ifirstf1,ilastf1 ",ifirstf1,ilastf1
c               write(6,*) "if0,if1 ",if0,if1
c               write(6,*) "ic0,ic1 ",ic0,ic1
c               write(6,*) "fine energy = ",
c    &   conservc(ic0,ic1) + slope0(ic0,ic1)*deltax(ir0,0)
c    &                     +  slope1(ic0,ic1)*deltax1
c               write(6,*) "densf(if0,if1) ",densf(if0,if1)
c               write(6,*) "vel() ",velf(if0,if1,0),velf(if0,if1,1)
c               write(6,*) "presf(if0,if1) ",presf(if0,if1)
c               write(6,*) "conservc(",ic0,ic1,") ",conservc(ic0,ic1)
c               write(6,*) "slope0(",ic0,ic1,") ",slope0(ic0,ic1)
c               write(6,*) "deltax(",ir0,0,") ",deltax(ir0,0)
c               write(6,*) "slope1(",ic0,ic1,") ",slope1(ic0,ic1)
c               write(6,*) "deltax1(",ir1,1,") ",deltax1
c               write(6,*) "val ",val,", v2norm ",v2norm
c               write(6,*) "presc(",ic0,ic1,") ",presc(ic0,ic1)
c               write(6,*) "densc(",ic0,ic1,") ",densc(ic0,ic1)
c               write(6,*) "velc() ",velc(ic0,ic1,0),velc(ic0,ic1,1)
                presneg = .true.
            endif
        enddo
      enddo

      if (presneg) then
        write(6,*) "negative pressure reported in conservlinint2d"
c       write(6,*) "coarse conserved values"
c       do ic1=ifirstc1-1,ilastc1+1
c          do ic0=ifirstc0-1,ilastc0+1
c             write(6,*) "ic0,ic1,conservc = ",
c    &        ic0,ic1,conservc(ic0,ic1)
c          enddo
c       enddo
c       write(6,*) "coarse density values"
c       do ic1=ifirstc1-1,ilastc1+1
c          do ic0=ifirstc0-1,ilastc0+1
c             write(6,*) "ic0,ic1,densc = ",
c    &        ic0,ic1,densc(ic0,ic1)
c          enddo
c       enddo
c       write(6,*) "coarse velocity0 values"
c       do ic1=ifirstc1-1,ilastc1+1
c          do ic0=ifirstc0-1,ilastc0+1
c             write(6,*) "ic0,ic1,velc0 = ",
c    &        ic0,ic1,velc(ic0,ic1,0)
c          enddo
c       enddo
c       write(6,*) "coarse velocity1 values"
c       do ic1=ifirstc1-1,ilastc1+1
c          do ic0=ifirstc0-1,ilastc0+1
c             write(6,*) "ic0,ic1,velc1 = ",
c    &        ic0,ic1,velc(ic0,ic1,1)
c          enddo
c       enddo
c       write(6,*) "coarse pressure values"
c       do ic1=ifirstc1-1,ilastc1+1
c          do ic0=ifirstc0-1,ilastc0+1
c             write(6,*) "ic0,ic1,presc = ",
c    &        ic0,ic1,presc(ic0,ic1)
c          enddo
c       enddo
c       write(6,*) "fine density values"
c       do if1=ifirstf1,ilastf1
c          do if0=ifirstf0,ilastf0
c              write(6,*) "if0,if1,densf = ",
c    &         if0,if1,densf(if0,if1)
c          enddo
c        enddo
c       write(6,*) "fine velocity0 values"
c       do if1=ifirstf1,ilastf1
c          do if0=ifirstf0,ilastf0
c              write(6,*) "if0,if1,velf0 = ",
c    &         if0,if1,velf(if0,if1,0)
c          enddo
c       enddo
c       write(6,*) "fine velocity1 values"
c       do if1=ifirstf1,ilastf1
c          do if0=ifirstf0,ilastf0
c              write(6,*) "if0,if1,velf1 = ",
c    &         if0,if1,velf(if0,if1,1)
c          enddo
c       enddo
c        call flush(6)
        stop
      endif
c
      return
      end
c
c
c***********************************************************************
c Volume weighted averaging for 2d flux variables.
c***********************************************************************
c
      subroutine conservavg2d(
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
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
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:NDIM-1)
      REAL
     &  dxf(0:NDIM-1),
     &  dxc(0:NDIM-1),
     &  gamma
      REAL
     &  densf(CELL2d(filo,fihi,0)),
     &  densc(CELL2d(cilo,cihi,0)),
     &  velf(CELL2d(filo,fihi,0),0:NDIM-1),
     &  presf(CELL2d(filo,fihi,0))
      REAL
     &  velc(CELL2d(cilo,cihi,0),0:NDIM-1),
     &  presc(CELL2d(cilo,cihi,0))
      REAL
     &  conservf(CELL2d(ifirstf,ilastf,0)) 
      integer ic0,ic1,if0,if1,ir0,ir1,id
      REAL dVf,dVc,val,valinv,v2norm
c
c***********************************************************************
c
      dVf = dxf(0)*dxf(1)
      dVc = dxc(0)*dxc(1)
c
c construct coarse velocity values using conservative average
c of momentum
c
      do id=0,NDIM-1 

         do if1=ifirstf1,ilastf1
            do if0=ifirstf0,ilastf0
               conservf(if0,if1) = densf(if0,if1)*velf(if0,if1,id)
            enddo
         enddo

         do ic1=ifirstc1,ilastc1
            do ic0=ifirstc0,ilastc0
               velc(ic0,ic1,id) = zero
            enddo
         enddo

         do ir1=0,ratio(1)-1
            do ir0=0,ratio(0)-1
               do ic1=ifirstc1,ilastc1
                  if1=ic1*ratio(1)+ir1
                  do ic0=ifirstc0,ilastc0
                     if0=ic0*ratio(0)+ir0
                     velc(ic0,ic1,id) = velc(ic0,ic1,id)
     &                                + conservf(if0,if1)*dVf
                  enddo
               enddo
            enddo
         enddo

         do ic1=ifirstc1,ilastc1
            do ic0=ifirstc0,ilastc0
               velc(ic0,ic1,id) = velc(ic0,ic1,id)
     &                          / (densc(ic0,ic1)*dVc)
           enddo
         enddo

      enddo

c
c construct coarse velocity values using conservative average
c of energy
c
      val = (gamma-one)
      valinv = one/(gamma-one)

      do if1=ifirstf1,ilastf1
         do if0=ifirstf0,ilastf0
            v2norm = velf(if0,if1,0)**2 + velf(if0,if1,1)**2
            conservf(if0,if1) = presf(if0,if1)*valinv
     &                        + half*densf(if0,if1)*v2norm
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            presc(ic0,ic1) = zero
         enddo
      enddo

      do ir1=0,ratio(1)-1
         do ir0=0,ratio(0)-1
            do ic1=ifirstc1,ilastc1
               if1=ic1*ratio(1)+ir1
               do ic0=ifirstc0,ilastc0
                  if0=ic0*ratio(0)+ir0
                  presc(ic0,ic1) = presc(ic0,ic1)
     &                           + conservf(if0,if1)*dVf
               enddo
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            v2norm = velc(ic0,ic1,0)**2 + velc(ic0,ic1,1)**2
            presc(ic0,ic1) = ((presc(ic0,ic1)/dVc)
     &                     - half*densc(ic0,ic1)*v2norm) * val
        enddo
      enddo
c
      return
      end

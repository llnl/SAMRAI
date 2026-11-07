c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file defining routine to compute flattening
c                coefficient in 2d.
c
c
c***********************************************************************
c flatten routine
c***********************************************************************
c
      subroutine amrflaten2d(ifirst,ilast,
     &                  cilo,cihi,
     &                  mc,gamma,
     &                  densc,presc,velc,flattn,flat2,sound)
c***********************************************************************
c description of arguments:
c   input:
c     p        ==> pressure            tracest(NEQU)
c     vn       ==> normal velocity     tracest(2+idir)
c   output:
c     flattn   ==> flattening coefficient
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
cnclude(FORTDIR/../const.i)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer ifirst,ilast,
     &  cilo,cihi,
     &  mc
      REAL gamma
      REAL
     &  densc(cilo:cilo+mc),
     &  velc(cilo:cilo+mc),
     &  presc(cilo:cilo+mc)
      REAL  flattn(ifirst:ifirst+mc)
      REAL  flat2(cilo:cilo+mc),sound(cilo:cilo+mc)

      REAL zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)

c
c***********************************************************************
c


      REAL shktst, zcut1, zcut2, dzcut, smallp
      REAL denom, zeta, tst, tmp, dp,zi
      REAL rhoc1,rhoc2,chi
      integer ic

c
c***********************************************************************
c
      shktst =0.33d0
      zcut1 =0.75d0
      zcut2 =0.85d0
      smallp =1.0d-6
      dzcut = one/(zcut2-zcut1)

      do ic=ifirst,ilast    
         flattn(ic) = one
      enddo
      do ic=ifirst-1,ilast+1   
         sound(ic) = sqrt(gamma*presc(ic)/densc(ic))
      enddo

      do ic=ifirst+1,ilast-1
         dp = presc(ic+1) - presc(ic-1)
         denom = max(smallp,abs(presc(ic+2)-presc(ic-2)))
         zeta = abs(dp)/denom
         zi = min( one, max( zero, dzcut*(zeta - zcut1) ) )
         if ((velc(ic-1) - velc(ic+1)).gt.zero) then
            tst = one
         else
            tst = zero
         endif
         rhoc1 = densc(ic+1)*sound(ic+1)**2
         rhoc2 = densc(ic-1)*sound(ic-1)**2
         tmp = min(rhoc1,rhoc2)
         if ((abs(dp)/tmp).gt.shktst) then
            chi = tst
         else
            chi = zero
         endif
         flat2(ic) = chi*zi
      enddo
      do ic=ifirst+2,ilast-2
         flattn(ic) = one - max(flat2(ic-1),flat2(ic),flat2(ic+1))
      enddo
      
      return
      end

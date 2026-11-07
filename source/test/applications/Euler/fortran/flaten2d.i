c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file defining routine to compute flattening
c                coefficient for 2d euler equation.
c
c***********************************************************************
      subroutine flaten2d(ifirst,ilast,i,mc,idir,
     &                  tracest,sound,flattn)
c***********************************************************************
c Compute "flattn" (flattening coefficient) for cell i
c description of arguments:
c   input:
c     ifirst,ilast,mc  ==> array dimensions 
c     i                ==> cell index flattn is being computed on
c     idir             ==> coordinate direction
c     tracest          ==> array of traced states
c     sound            ==> array of sound speeds
c   output:
c     flattn   ==> flattening coefficient
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
include(FORTDIR/probparams.i)dnl
include(FORTDIR/const.i)dnl
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer ifirst,ilast,i,mc,idir
      REAL tracest(ifirst-FACEG:ifirst+mc+FACEG,NEQU)
      REAL  sound(ifirst-CELLG:ifirst+mc-1+CELLG)
      REAL flat2(3)
      REAL flattn

      REAL shktst, zcut1, zcut2, dzcut, smallp
      REAL denom, zeta, tst, tmp, dp,zi
      REAL rhoc1,rhoc2,chi
      integer j,ii
c
c***********************************************************************
c
      shktst =0.33d0
      zcut1 =0.75d0
      zcut2 =0.85d0
      smallp =1.0d-6
      dzcut = one/(zcut2-zcut1)
      
      flattn = one 

c
c  Leave flattn = 1. at bounds.
c
       if ((i.lt.ifirst-1) .or. (i.gt.ilast+1)) return
c
c  Compute "flat2" at i-1,i,i+1
c
      do j=1,3
         ii = (i-2) + j  
         dp = tracest(ii+1,NEQU) - tracest(ii-1,NEQU)
         denom = max(smallp,abs(tracest(ii+2,NEQU)-tracest(ii-2,NEQU)))
         zeta = abs(dp)/denom
         zi = min( one, max( zero, dzcut*(zeta - zcut1) ) )
         if ((tracest(ii-1,idir+2) - tracest(ii+1,idir+2)).gt.zero) then
            tst = one
         else 
            tst = zero
         endif
         rhoc1 = tracest(ii+1,1)*sound(ii+1)**2
         rhoc2 = tracest(ii-1,1)*sound(ii-1)**2
         tmp = min(rhoc1,rhoc2)
         if ((abs(dp)/tmp).gt.shktst) then
            chi = tst
         else 
            chi = zero
         endif
         flat2(j) = chi*zi
      enddo
c
c  Compute "flattn" at cell i, using flat2 at i-1,i,i+1
c
      flattn = one - max(flat2(1),flat2(2),flat2(3))

      return
      end

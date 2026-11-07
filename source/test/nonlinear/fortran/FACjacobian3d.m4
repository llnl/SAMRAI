c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines in 3d.
c
define(NDIM,3)dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine compfacdiag3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  dt,cellvol,
     &  expu,dsrcdu,
     &  diag)
c***********************************************************************
      implicit none
      double precision one
      parameter(one=1.d0)
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      double precision dt, cellvol
      double precision 
     &  expu(CELL3d(ifirst,ilast,0)),
     &  dsrcdu(CELL3d(ifirst,ilast,0))
c ouput arrays:
      double precision 
     &  diag(CELL3d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer ic0,ic1,ic2
      double precision dtfrac

      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            do ic0=ifirst0,ilast0
               diag(ic0,ic1,ic2) = one - 
     &                      dt*(expu(ic0,ic1,ic2) + dsrcdu(ic0,ic1,ic2))
            enddo
         enddo
      enddo

      return
      end
c
c
c
      subroutine compfacoffdiag3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  dt, cellvol,
     &  sidediff0,sidediff1,sidediff2,
     &  offdiag0,offdiag1,offdiag2)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      double precision dt, cellvol
      double precision
     &  sidediff0(SIDE3d0(ifirst,ilast,0)),
     &  sidediff1(SIDE3d1(ifirst,ilast,0)),
     &  sidediff2(SIDE3d2(ifirst,ilast,0))
c ouput arrays:
      double precision
     &  offdiag0(FACE3d0(ifirst,ilast,0)),
     &  offdiag1(FACE3d1(ifirst,ilast,0)),
     &  offdiag2(FACE3d2(ifirst,ilast,0))
c
c***********************************************************************
c
      integer ic0,ic1,ic2,ie0,ie1,ie2
      double precision factor

      factor = dt
c
      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            do ie0=ifirst0,ilast0+1
               offdiag0(ie0,ic1,ic2) = factor*sidediff0(ie0,ic1,ic2)
            enddo
         enddo
      enddo
c
      do ic0=ifirst0,ilast0
         do ic2=ifirst2,ilast2
            do ie1=ifirst1,ilast1+1
               offdiag1(ie1,ic2,ic0) = factor*sidediff1(ic0,ie1,ic2)
            enddo
         enddo
      enddo
c
      do ic1=ifirst1,ilast1
         do ic0=ifirst0,ilast0
            do ie2=ifirst2,ilast2+1
               offdiag2(ie2,ic0,ic1) = factor*sidediff2(ic0,ic1,ie2)
            enddo
         enddo
      enddo
c
      return
      end


      subroutine compjv3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  gcw,
     &  diag,
     &  flux0,flux1,flux2,
     &  v,
     &  dx, dt,
     &  jv
     &  )
c
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer gcw
      double precision 
     &  v(CELL3d(ifirst,ilast,gcw))
      double precision 
     &  diag(CELL3d(ifirst,ilast,0))
      double precision
     &  flux0(SIDE3d0(ifirst,ilast,0)),
     &  flux1(SIDE3d1(ifirst,ilast,0)),
     &  flux2(SIDE3d2(ifirst,ilast,0))
      double precision
     &  jv(CELL3d(ifirst,ilast,0))
      double precision dx(0:NDIM), dt
      integer ic0,ic1,ic2
c
      do ic2=ifirst2,ilast2
         do ic1=ifirst1,ilast1
            do ic0=ifirst0,ilast0
               jv(ic0,ic1,ic2) = v(ic0,ic1,ic2) - 
     &              dt*( diag(ic0,ic1,ic2)*v(ic0,ic1,ic2)
     &              + (flux0(ic0+1,ic1,ic2)-flux0(ic0,ic1,ic2))/dx(0)
     &              + (flux1(ic0,ic1+1,ic2)-flux1(ic0,ic1,ic2))/dx(1)
     &              + (flux2(ic0,ic1,ic2+1)-flux2(ic0,ic1,ic2))/dx(2)
     &              )
            enddo
         enddo
      enddo
c
      return
      end

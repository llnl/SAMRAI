c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines in 2d.
c
define(NDIM,2)dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine compfacdiag2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  dt,cellvol,
     &  expu,dsrcdu,
     &  diag)
c***********************************************************************
      implicit none
      double precision one
      parameter(one=1.d0)
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      double precision dt, cellvol
      double precision 
     &  expu(CELL2d(ifirst,ilast,0)),
     &  dsrcdu(CELL2d(ifirst,ilast,0))
c ouput arrays:
      double precision 
     &  diag(CELL2d(ifirst,ilast,0))
c
c***********************************************************************
c
      integer ic0,ic1

      do ic1=ifirst1,ilast1
         do ic0=ifirst0,ilast0
            diag(ic0,ic1) = one - dt*(expu(ic0,ic1) + dsrcdu(ic0,ic1))
         enddo
      enddo

      return
      end
c
c
c
      subroutine compfacoffdiag2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  dt, cellvol,
     &  sidediff0,sidediff1,
     &  offdiag0,offdiag1)
c***********************************************************************
      implicit none
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      double precision dt, cellvol
      double precision
     &  sidediff0(SIDE2d0(ifirst,ilast,0)),
     &  sidediff1(SIDE2d1(ifirst,ilast,0))
c ouput arrays:
      double precision
     &  offdiag0(FACE2d0(ifirst,ilast,0)),
     &  offdiag1(FACE2d1(ifirst,ilast,0))
c
c***********************************************************************
c
      integer ic0,ic1,ie0,ie1
      double precision factor

      factor =  dt
c
      do ic1=ifirst1,ilast1
         do ie0=ifirst0,ilast0+1
            offdiag0(ie0,ic1) = factor*sidediff0(ie0,ic1)
         enddo
      enddo
c
      do ic0=ifirst0,ilast0
         do ie1=ifirst1,ilast1+1
            offdiag1(ie1,ic0) = factor*sidediff1(ic0,ie1)
         enddo
      enddo
c
      return
      end


      subroutine compjv2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw,
     &  diag,
     &  flux0,flux1,
     &  v,
     &  dx, dt,
     &  jv
     &  )
c
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw
      double precision 
     &  v(CELL2d(ifirst,ilast,gcw))
      double precision 
     &  diag(CELL2d(ifirst,ilast,0))
      double precision
     &  flux0(SIDE2d0(ifirst,ilast,0)),
     &  flux1(SIDE2d1(ifirst,ilast,0))
      double precision
     &  jv(CELL2d(ifirst,ilast,0))
      double precision dx(0:NDIM), dt
      integer ic0,ic1
c
      do ic1=ifirst1,ilast1
         do ic0=ifirst0,ilast0
            jv(ic0,ic1) = v(ic0,ic1) - 
     &           dt*( diag(ic0,ic1)*v(ic0,ic1)
     &           + ( flux0(ic0+1,ic1) - flux0(ic0,ic1) )/dx(0)
     &           + ( flux1(ic0,ic1+1) - flux1(ic0,ic1) )/dx(1)
     &           )
         enddo
      enddo
c
      return
      end

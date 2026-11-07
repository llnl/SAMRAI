c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to compute error in solution at time t in 3d.
c
define(NDIM,3)dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine error3d(
     & lo0, hi0, lo1, hi1, lo2, hi2,
     & u, w,
     & lambda,
     & xlo, xhi, dx,
     & t,
     & maxerror,
     & l2error )

c  Compute error in solution at time t.

      implicit none

      integer lo0
      integer lo1
      integer lo2

      integer hi0
      integer hi1
      integer hi2

      double precision lambda
      double precision t
      double precision l2error
      double precision maxerror

      double precision u(CELL3d(lo,hi,0))
      double precision w(CELL3d(lo,hi,0))

      double precision xlo(0:NDIM-1)
      double precision xhi(0:NDIM-1)
      double precision dx(0:NDIM-1)

      integer i
      integer j
      integer k

      double precision diff
      double precision error_ijk
      double precision localerror

      double precision xi, xterm
      double precision yj, yterm
      double precision zk, zterm

      double precision zero,       half,       one
      parameter      ( zero=0.0d0, half=0.5d0, one=1.0d0 )

      intrinsic DABS, MAX
      double precision DABS, MAX

      localerror = zero
      l2error = zero
      zk = xlo(2) + dx(2)*half
      do k = lo2, hi2
         zterm = zk*(one - zk)
         yj = xlo(1) + dx(1)*half
         do j = lo1, hi1
            yterm = yj*(one - yj)
            xi = xlo(0) + dx(0)*half
            do i = lo0, hi0
               xterm = xi*(one - xi)
               error_ijk = ABS(u(i,j,k) - t*xterm*yterm*zterm)
               l2error = l2error + w(i,j,k)*error_ijk**2
               diff = w(i,j,k)*error_ijk
               localerror = MAX(localerror, diff)
               xi = xi + dx(0)
            end do
            yj = yj + dx(1)
         end do
         zk = zk + dx(2)
      end do

c  Since the weight was used to mask coarse cells that have been refined,
c  the error in those that aren't covered are scaled by the cell volume.

      localerror = localerror/(dx(0)*dx(1)*dx(2))
      l2error = sqrt(l2error)

      maxerror = MAX(localerror, maxerror)

      return
      end

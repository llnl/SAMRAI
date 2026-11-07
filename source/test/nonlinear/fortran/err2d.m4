c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to compute error in solution at time t in 2d.
c
define(NDIM,2)dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine error2d(
     & lo0, hi0, lo1, hi1,
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

      integer hi0
      integer hi1

      double precision lambda
      double precision t
      double precision l2error
      double precision maxerror

      double precision u(CELL2d(lo,hi,0))
      double precision w(CELL2d(lo,hi,0))

      double precision xlo(0:NDIM-1)
      double precision xhi(0:NDIM-1)
      double precision dx(0:NDIM-1)

      integer i
      integer j

      double precision diff
      double precision error_ij
      double precision localerror
      double precision solmax

      double precision xi, xterm
      double precision yj, yterm

      double precision zero,       half,       one
      parameter      ( zero=0.0d0, half=0.5d0, one=1.0d0 )

      intrinsic DABS, MAX
      double precision DABS, MAX

      localerror = zero
      l2error = zero
      solmax = zero
      yj = xlo(1) + dx(1)*half
      do j = lo1, hi1
         yterm = yj*(one - yj)
         xi = xlo(0) + dx(0)*half
         do i = lo0, hi0
            xterm = xi*(one - xi)
            solmax = MAX(solmax,t*xterm*yterm)
            error_ij = ABS(u(i,j) - t*xterm*yterm)
            l2error = l2error + w(i,j)*error_ij**2
            diff = w(i,j)*error_ij
            localerror = MAX(localerror, diff)
            xi = xi + dx(0)
         end do
         yj = yj + dx(1)
      end do

c  Since the weight was used to mask coarse cells that have been refined,
c  the error in those that aren't covered are scaled by the cell volume.
c  The max error is scaled by the max of the exact solution.

      localerror = localerror/(dx(0)*dx(1)*solmax)
      l2error = sqrt(l2error)

      maxerror = MAX(localerror, maxerror)

      return
      end

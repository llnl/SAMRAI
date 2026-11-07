c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to take an RK step in 2d.
c
define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine rkstep2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,
     &  dt, alpha_1, alpha_2, beta,
     &  conv_coeff,
     &  diff_coeff,
     &  src_coeff,
     &  soln_updated,
     &  soln_fixed,
     &  rhs,
     &  nequ)
c***********************************************************************
      implicit none
include(FORTDIR/const.i)dnl
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0,gcw1
      integer nequ

      REAL    dt, alpha_1, alpha_2, beta  
      REAL    conv_coeff(0:NDIM-1)
      REAL    diff_coeff, src_coeff
c
c variables in 2d cell indexed
      REAL
     &     soln_updated(CELL2dVECG(ifirst,ilast,gcw),0:nequ-1),
     &     soln_fixed(CELL2d(ifirst,ilast,0),0:nequ-1),
     &     rhs(CELL2dVECG(ifirst,ilast,gcw),0:nequ-1)
c
c***********************************************************************
c***********************************************************************
c
      integer ic0,ic1,ineq
c
      do ic1=ifirst1,ilast1
         do ic0=ifirst0,ilast0
            do ineq=0,nequ-1
               soln_updated(ic0,ic1,ineq) = 
     &            alpha_1*soln_fixed(ic0,ic1,ineq) + 
     &            alpha_2*soln_updated(ic0,ic1,ineq) +
     &            beta*dt*rhs(ic0,ic1,ineq)
            end do
         end do
      end do

      return
      end

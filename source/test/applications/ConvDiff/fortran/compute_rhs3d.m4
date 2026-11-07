c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to compute RHS of 3d convection diffusion
c                equation.
c
define(NDIM,3)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine computerhs3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  gcw0,gcw1,gcw2,
     &  dx,
     &  conv_coeff,
     &  diff_coeff,
     &  src_coeff,
     &  var,
     &  rhs,
     &  nequ)
c***********************************************************************
      implicit none
include(FORTDIR/const.i)dnl
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer gcw0,gcw1,gcw2
      integer nequ

      REAL dx(0:NDIM-1)
      REAL conv_coeff(0:NDIM-1)
      REAL diff_coeff, src_coeff
c
c variables in 3d cell indexed         
      REAL
     &     var(CELL3dVECG(ifirst,ilast,gcw),0:nequ-1),
     &     rhs(CELL3dVECG(ifirst,ilast,gcw),0:nequ-1)
c
c***********************************************************************     
c
      integer ic0,ic1,ic2,ineq
      REAL conv_term_x, conv_term_y, conv_term_z, conv_term
      REAL diff_term_x, diff_term_y, diff_term_z, diff_term
c
c Compute the RHS for the convection-diffusion equation:
c    du/dt + div(a*u) = mu div^2(u) + gamma
c
c Compute as:
c    du/dt = RHS = -div(a*u) + mu div^2(u) + gamma
c
      do ic2=ifirst2,ilast2
        do ic1=ifirst1,ilast1
          do ic0=ifirst0,ilast0
            do ineq=0,nequ-1

c
c  2nd order accurate difference for convective terms (div(au))
c
             conv_term_x =
     &         ( var(ic0+1,ic1,ic2,ineq) - var(ic0-1,ic1,ic2,ineq) ) /
     &         ( 2.*dx(0) )

             conv_term_y =
     &         ( var(ic0,ic1+1,ic2,ineq) - var(ic0,ic1-1,ic2,ineq) ) /
     &         ( 2.*dx(1) )

             conv_term_z =
     &         ( var(ic0,ic1,ic2+1,ineq) - var(ic0,ic1,ic2-1,ineq) ) /
     &         ( 2.*dx(2) )

             conv_term = conv_coeff(0)*conv_term_x
     &                   + conv_coeff(1)*conv_term_y
     &                   + conv_coeff(2)*conv_term_z

c
c  2nd order accurate difference for diffusive terms (div^2(u))
c
             diff_term_x = (   var(ic0+1,ic1,ic2,ineq)
     &                     - 2*var(ic0,ic1,ic2,ineq)
     &                     +   var(ic0-1,ic1,ic2,ineq) )
     &                     / dx(0)**2

             diff_term_y = (   var(ic0,ic1+1,ic2,ineq)
     &                     - 2*var(ic0,ic1,ic2,ineq)
     &                     +   var(ic0,ic1-1,ic2,ineq) )
     &                     / dx(1)**2

             diff_term_z = (   var(ic0,ic1,ic2+1,ineq)
     &                     - 2*var(ic0,ic1,ic2,ineq)
     &                     +   var(ic0,ic1,ic2-1,ineq) )
     &                     / dx(2)**2

             diff_term   =  diff_coeff*
     &                       (diff_term_x + diff_term_y + diff_term_z)

c
c  function evaluation (RHS) 
c
             rhs(ic0,ic1,ic2,ineq) = -1.*(conv_term) 
     &          + diff_term 
     &          + src_coeff

            end do
          end do
        end do
      end do

      return
      end

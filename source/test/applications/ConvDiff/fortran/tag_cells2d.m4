c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routine to tag cells in 2d.
c
define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl

      subroutine tagcells2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,
     &  tags,
     &  var,
     &  refine_tag_val,
     &  tolerance,
     &  nequ)
c***********************************************************************
c***********************************************************************     
      implicit none
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0,gcw1
      integer refine_tag_val
      integer nequ
      REAL    tolerance(0:nequ-1)
c
c variables in 2d cell indexed         
      integer 
     &     tags(CELL2d(ifirst,ilast,0))
      REAL
     &     var(CELL2dVECG(ifirst,ilast,gcw),0:nequ-1)
c
      integer ic0,ic1,ineq
c
c***********************************************************************     

c//// Constant value - tag wherever var > tol //////
       do ic1=ifirst1,ilast1
          do ic0=ifirst0,ilast0
             do ineq=0,nequ-1
                tags(ic0,ic1) = 0
                if (var(ic0,ic1,ineq) .gt. tolerance(ineq)) 
     &            tags(ic0,ic1) = refine_tag_val
             end do
          end do
       end do

c//// Gradient, double sided  - tag wherever grad > tol //////
c      do ic1=ifirst1,ilast1
c         do ic0=ifirst0,ilast0
c            do ineq=0,nequ-1
c               tags(ic0,ic1) = 0
c               gradxp = abs(var(ic0+1,ic1,ineq) - var(ic0,ic1,ineq))
c               gradxm = abs(var(ic0,ic1,ineq)   - var(ic0-1,ic1,ineq))
c               gradyp = abs(var(ic0,ic1+1,ineq) - var(ic0,ic1,ineq))
c               gradym = abs(var(ic0,ic1,ineq)   - var(ic0,ic1-1,ineq))
c               if (gradxp .gt. tolerance(ineq))
c     &            tags(ic0,ic1) = refine_tag_val
c               if (gradxm .gt. tolerance(ineq))
c     &            tags(ic0,ic1) = refine_tag_val
c               if (gradyp .gt. tolerance(ineq))
c     &            tags(ic0,ic1) = refine_tag_val
c               if (gradym .gt. tolerance(ineq))
c     &            tags(ic0,ic1) = refine_tag_val
c            end do
c        end do
c      end do

      return
      end


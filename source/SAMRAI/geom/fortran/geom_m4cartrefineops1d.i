c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for 1d Cartesian refine operators.
c
define(NDIM,1)dnl
include(PDAT_FORTDIR/pdat_m4arrdim1d.i)dnl
include(FORTDIR/geom_m4cartopstuff.i)dnl
c
define(cart_refine_op_subroutine_head_1d,`dnl
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0
      integer ratio(0:NDIM-1)
      double precision
     &  dxc(0:NDIM-1),
     &  dxf(0:NDIM-1)
')dnl
c
define(cart_clinrefine_op_subroutine_head_1d,`dnl
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff0,slope0)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0
      integer ratio(0:NDIM-1)
      double precision
     &  dxc(0:NDIM-1),
     &  dxf(0:NDIM-1),
     &  deltax(0:15,0:NDIM-1)
')dnl
c
define(cart_linref_cell_body_1d,`dnl
c
c***********************************************************************
c
coarse_fine_cell_deltas(0)dnl

      do if0=ifirstf0,ilastf0
coarsen_index(if0,ic0,ratio(0))dnl
        ir0=if0-ic0*ratio(0)
        x=deltax(ir0,0)/dxc(0)
        if( x .lt. 0.d0 ) then
           ic0 = ic0-1
           x = x + one
        endif
        arrayf(if0)=arrayc(ic0)+(arrayc(ic0+1)-arrayc(ic0))*x
      enddo
c
      return
      end
')dnl
c
define(cart_clinref_cell_body_1d,`dnl
c
c***********************************************************************
c
coarse_fine_cell_deltas(0)dnl

muscl_limited_cell_slopes$1(0,ie0,ie0-1,ic0)dnl

      do if0=ifirstf0,ilastf0
coarsen_index(if0,ic0,ratio(0))dnl
         ir0=if0-ic0*ratio(0)
         arrayf(if0)=arrayc(ic0)+slope0(ic0)*deltax(ir0,0)
      enddo
c
      return
      end
')dnl
c
define(cart_clinref_edge_body_1d,`dnl
c
c***********************************************************************
c
coarse_fine_cell_deltas(0)dnl

muscl_limited_cell_slopes(0,ie0,ie0-1,ic0)dnl

      do if0=ifirstf0,ilastf0
coarsen_index(if0,ic0,ratio(0))dnl
         ir0=if0-ic0*ratio(0)
         arrayf(if0)=arrayc(ic0)+slope0(ic0)*deltax(ir0,0)
      enddo
c
      return
      end
')dnl
c
define(cart_clinref_face_body_1d,`dnl
c
c***********************************************************************
c

coarse_fine_face_deltas(0)dnl

muscl_limited_face_slopes(0,ic0+1,ic0,ie0)dnl

      do if0=ifirstf0,ilastf0+1
coarsen_index(if0,ic0,ratio(0))dnl
         ir0=if0-ic0*ratio(0)
         arrayf(if0)=arrayc(ic0)+slope0(ic0)*deltax(ir0,0)
      enddo
c
      return
      end
')dnl
c
define(cart_linref_node_body_1d,`dnl
c
c***********************************************************************
c
      realrat=one/dble(ratio(0))

      do ic0=ifirstc0,ilastc0
         if0=ic0*ratio(0)
         if (if0.ge.filo0.and.if0.le.fihi0+1) then 
            do ir0=0,ratio(0)-1
               ie0=if0+ir0
               x = dble(ir0)*realrat
               if (ie0.ge.filo0.and.ie0.le.fihi0+1) then
                  arrayf(ie0) = arrayc(ic0)*(one-x)
     &                              + arrayc(ic0+1)*x
               endif
            enddo
         endif
      enddo
c
      ic0 = ilastc0+1
      if0 = ic0*ratio(0)
      if (if0.ge.filo0.and.if0.le.fihi0+1) then
         arrayf(if0) = arrayc(ic0)
      endif
c
      return
      end
')dnl
c 
define(cart_linref_op_cell_1d,`dnl
cart_refine_op_subroutine_head_1d()dnl
      $1
     &  arrayc(CELL1d(cilo,cihi,0)),
     &  arrayf(CELL1d(filo,fihi,0))
      double precision deltax(0:15,0:NDIM-1),x
      integer ic0,if0,ir0
cart_linref_cell_body_1d()dnl
')dnl
define(cart_clinref_op_cell_1d,`dnl
cart_clinrefine_op_subroutine_head_1d()dnl
      $1
     &  arrayc(CELL1d(cilo,cihi,0)),
     &  arrayf(CELL1d(filo,fihi,0)),
     &  diff0(cilo0:cihi0+1),
     &  slope0(CELL1d(cilo,cihi,0))
      integer ic0,ie0,if0,ir0
ifelse($1,`double complex',`dnl
      double precision
     &  coef2real,coef2imag,boundreal,boundimag,
     &  diff0real,diff0imag,diff1real,diff1imag,
     &  slopereal,slopeimag
',`dnl
      $1
     &  coef2,bound
')dnl
ifelse($1,`double complex',`dnl
cart_clinref_cell_body_1d(`_complex')dnl
',`dnl
cart_clinref_cell_body_1d(`')dnl
')dnl
')dnl
define(cart_clinref_op_edge_1d,`dnl
cart_clinrefine_op_subroutine_head_1d()dnl
      $1
     &  arrayc(EDGE1d(cilo,cihi,0)),
     &  arrayf(EDGE1d(filo,fihi,0)),
     &  diff0(cilo0:cihi0+1),
     &  slope0(EDGE1d(cilo,cihi,0))
      integer ic0,ie0,if0,ir0
      $1
     &  coef2,bound
cart_clinref_edge_body_1d()dnl
')dnl
define(cart_clinref_op_face_1d,`dnl
cart_clinrefine_op_subroutine_head_1d()dnl
      $1
     &  arrayc(FACE1d(cilo,cihi,0)),
     &  arrayf(FACE1d(filo,fihi,0)),
     &  diff0(cilo0:cihi0),
     &  slope0(FACE1d(cilo,cihi,0))
      integer ic0,ie0,if0,ir0
      $1
     &  coef2,bound
cart_clinref_face_body_1d()dnl
')dnl
define(cart_linref_op_node_1d,`dnl
cart_refine_op_subroutine_head_1d()dnl
      $1
     &  arrayc(NODE1d(cilo,cihi,0)),
     &  arrayf(NODE1d(filo,fihi,0))
      double precision realrat,x
      integer i,ic0,if0,ie0,ir0
cart_linref_node_body_1d()dnl
')dnl
define(cart_clinref_op_side_1d,`dnl
cart_clinrefine_op_subroutine_head_1d()dnl
      $1
     &  arrayc(SIDE1d(cilo,cihi,0)),
     &  arrayf(SIDE1d(filo,fihi,0)),
     &  diff0(cilo0:cihi0),
     &  slope0(SIDE1d(cilo,cihi,0))
      integer ic0,ie0,if0,ir0
      $1
     &  coef2,bound
cart_clinref_face_body_1d()dnl
')dnl
define(cart_clinref_op_outerface_1d,`dnl
cart_refine_op_subroutine_head_1d()dnl
      $1
     &  arrayc(OUTERFACE1d(cilo,cihi,0)),
     &  arrayf(OUTERFACE1d(filo,fihi,0))
c
c***********************************************************************
c
      arrayf(1)=arrayc(1)
c
      return
      end
')dnl

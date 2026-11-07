c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for 2d Cartesian refine operators.
c
define(NDIM,2)dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl
include(FORTDIR/geom_m4cartopstuff.i)dnl
c
define(cart_refine_op_subroutine_head_2d,`dnl
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:NDIM-1)
      double precision
     &  dxc(0:NDIM-1),
     &  dxf(0:NDIM-1)
')dnl
c
define(cart_clinrefine_op_subroutine_head_2d,`dnl
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff$1,slope$1,diff$2,slope$2)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:NDIM-1)
      double precision
     &  dxc(0:NDIM-1),
     &  dxf(0:NDIM-1),
     &  deltax(0:15,0:NDIM-1)
')dnl
c
define(cart_linref_cell_body_2d,`dnl
c
c***********************************************************************
c
coarse_fine_cell_deltas(0)dnl

coarse_fine_cell_deltas(1)dnl

      do if1=ifirstf1,ilastf1
coarsen_index(if1,ic1,ratio(1))dnl
         ir1=if1-ic1*ratio(1)
         y=deltax(ir1,1)/dxc(1)
         if( y .lt. 0.d0 ) then
            ic1 = ic1-1
            y = y + one
         endif
         do if0=ifirstf0,ilastf0
coarsen_index(if0,ic0,ratio(0))dnl
            ir0=if0-ic0*ratio(0)
            x=deltax(ir0,0)/dxc(0)
            if( x .lt. 0.d0 ) then
               ic0 = ic0-1
               x = x + one
            endif
            arrayf(if0,if1)=
     &      (arrayc(ic0,ic1)+(arrayc(ic0+1,ic1)-arrayc(ic0,ic1))*x)
     &        *(one-y)
     &      +(arrayc(ic0,ic1+1)
     &        +(arrayc(ic0+1,ic1+1)-arrayc(ic0,ic1+1))*x)*y
         enddo
      enddo
c
      return
      end
')dnl
c
define(cart_clinref_cell_body_2d,`dnl
c
c***********************************************************************
c
coarse_fine_cell_deltas(0)dnl

coarse_fine_cell_deltas(1)dnl

      do ic1=ifirstc1,ilastc1
muscl_limited_cell_slopes$1(0,`ie0,ic1',`ie0-1,ic1',`ic0,ic1')dnl
      enddo

      do ic0=ifirstc0,ilastc0
muscl_limited_cell_slopes$1(1,`ic0,ie1',`ic0,ie1-1',`ic0,ic1')dnl
      enddo

      do if1=ifirstf1,ilastf1
coarsen_index(if1,ic1,ratio(1))dnl
         ir1=if1-ic1*ratio(1)
         deltax1=deltax(ir1,1)
         do if0=ifirstf0,ilastf0
coarsen_index(if0,ic0,ratio(0))dnl
            ir0=if0-ic0*ratio(0)
            arrayf(if0,if1)=arrayc(ic0,ic1)
     &                      +slope0(ic0,ic1)*deltax(ir0,0)
     &                      +slope1(ic0,ic1)*deltax1
         enddo
      enddo
c
      return
      end
')dnl
c
define(cart_clinref_edge_body_2d,`dnl
c
c***********************************************************************
c

ifelse($1,`0',`
coarse_fine_cell_deltas(0)dnl
',`
coarse_fine_face_deltas(0)dnl
')dnl

ifelse($1,`1',`
coarse_fine_cell_deltas(1)dnl
',`
coarse_fine_face_deltas(1)dnl
')dnl

      do ic1=ifirstc1,ilastc1+$2
ifelse($1,`0',`
muscl_limited_cell_slopes(0,`ie0,ic1',`ie0-1,ic1',`ic0,ic1')dnl
',`
muscl_limited_face_slopes(0,`ic0+1,ic1',`ic0,ic1',`ie0,ic1')dnl
')dnl
      enddo

      do ic0=ifirstc0,ilastc0+$1
ifelse($1,`1',`
muscl_limited_cell_slopes(1,`ic0,ie1',`ic0,ie1-1',`ic0,ic1')dnl
',`
muscl_limited_face_slopes(1,`ic0,ic1+1',`ic0,ic1',`ic0,ie1')dnl
')dnl
      enddo

      do if1=ifirstf1,ilastf1+$2
coarsen_index(if1,ic1,ratio(1))dnl
         ir1=if1-ic1*ratio(1)
         deltax1=deltax(ir1,1)
         do if0=ifirstf0,ilastf0+$1
coarsen_index(if0,ie0,ratio(0))dnl
            ir0=if0-ie0*ratio(0)
            arrayf(if0,if1)=arrayc(ie0,ic1)
     &                +slope0(ie0,ic1)*deltax(ir0,0)
     &                +slope1(ie0,ic1)*deltax1
         enddo
      enddo
c
      return
      end
')dnl
c
define(cart_clinref_face_body_2d,`dnl
c
c***********************************************************************
c

coarse_fine_cell_deltas($2)dnl

coarse_fine_face_deltas($1)dnl

      do ic$2=ifirstc$2,ilastc$2
muscl_limited_face_slopes($1,`ic$1+1,ic$2',`ic$1,ic$2',`ie$1,ic$2')dnl
      enddo

      do ie$1=ifirstc$1,ilastc$1+1
muscl_limited_cell_slopes($2,`ie$1,ie$2',`ie$1,ie$2-1',`ie$1,ic$2')dnl
      enddo

      do if$2=ifirstf$2,ilastf$2
coarsen_index(if$2,ic$2,ratio($2))dnl
         ir$2=if$2-ic$2*ratio($2)
         deltax$2=deltax(ir$2,$2)
         do if$1=ifirstf$1,ilastf$1+1
coarsen_index(if$1,ie$1,ratio($1))dnl
            ir$1=if$1-ie$1*ratio($1)
            arrayf(if$1,if$2)=arrayc(ie$1,ic$2)
     &                +slope$1(ie$1,ic$2)*deltax(ir$1,$1)
     &                +slope$2(ie$1,ic$2)*deltax$2
         enddo
      enddo
c
      return
      end
')dnl
c
define(cart_linref_node_body_2d,`dnl
c
c***********************************************************************
c
      realrat0=one/dble(ratio(0))
      realrat1=one/dble(ratio(1))

node_refloop_pre(1)dnl
node_refloop_pre(0)dnl
               x = dble(ir0)*realrat0
               y = dble(ir1)*realrat1
               arrayf(ie0,ie1)=
     &              (arrayc(ic0,ic1)*(one-x) + 
     &               arrayc(ic0+1,ic1)*x)*(one-y) +
     &              (arrayc(ic0,ic1+1)*(one-x) + 
     &               arrayc(ic0+1,ic1+1)*x)*y
node_refloop_post()dnl
node_refloop_post()dnl
c
      return
      end
')dnl
c
define(cart_clinref_side_body_2d,`dnl
c
c***********************************************************************
c

ifelse($1,`0',`
coarse_fine_face_deltas(0)dnl
',`
coarse_fine_cell_deltas(0)dnl
')dnl

ifelse($1,`1',`
coarse_fine_face_deltas(1)dnl
',`
coarse_fine_cell_deltas(1)dnl
')dnl

      do ic1=ifirstc1,ilastc1+$1
ifelse($1,`0',`
muscl_limited_face_slopes(0,`ic0+1,ic1',`ic0,ic1',`ie0,ic1')dnl
',`
muscl_limited_cell_slopes(0,`ie0,ic1',`ie0-1,ic1',`ic0,ic1')dnl
')dnl
      enddo

      do ic0=ifirstc0,ilastc0+$2
ifelse($1,`1',`
muscl_limited_face_slopes(1,`ic0,ic1+1',`ic0,ic1',`ic0,ie1')dnl
',`
muscl_limited_cell_slopes(1,`ic0,ie1',`ic0,ie1-1',`ic0,ic1')dnl
')dnl
      enddo

      do if1=ifirstf1,ilastf1+$1
coarsen_index(if1,ic1,ratio(1))dnl
         ir1=if1-ic1*ratio(1)
         deltax1=deltax(ir1,1)
         do if0=ifirstf0,ilastf0+$2
coarsen_index(if0,ic0,ratio(0))dnl
            ir0=if0-ic0*ratio(0)
            arrayf(if0,if1)=arrayc(ic0,ic1)
     &                +slope0(ic0,ic1)*deltax(ir0,0)
     &                +slope1(ic0,ic1)*deltax1
         enddo
      enddo
c
      return
      end
')dnl
c
define(cart_clinref_outerface_body_2d,`dnl
c
c***********************************************************************
c

coarse_fine_cell_deltas($1)dnl

muscl_limited_cell_slopes($1,`ie$1',`ie$1-1',`ic$1')dnl

      do if$1=ifirstf$1,ilastf$1
coarsen_index(if$1,ic$1,ratio($1))dnl
         ir$1=if$1-ic$1*ratio($1)
         arrayf(if$1)=arrayc(ic$1)
     &                +slope$1(ic$1)*deltax(ir$1,$1)
      enddo
c
      return
      end
')dnl
c
define(cart_linref_op_cell_2d,`dnl
cart_refine_op_subroutine_head_2d()dnl
      $1
     &  arrayc(CELL2d(cilo,cihi,0)),
     &  arrayf(CELL2d(filo,fihi,0))
      double precision deltax(0:15,0:NDIM-1),x,y
      integer ic0,ic1,if0,if1,ir0,ir1
cart_linref_cell_body_2d()dnl
')dnl
define(cart_clinref_op_cell_2d,`dnl
cart_clinrefine_op_subroutine_head_2d(0,1)dnl
      $1
     &  arrayc(CELL2d(cilo,cihi,0)),
     &  arrayf(CELL2d(filo,fihi,0)),
     &  diff0(cilo0:cihi0+1),
     &  slope0(CELL2d(cilo,cihi,0)),
     &  diff1(cilo1:cihi1+1),
     &  slope1(CELL2d(cilo,cihi,0))
      integer ic0,ic1,ie0,ie1,if0,if1,ir0,ir1
ifelse($1,`double complex',`dnl
      double precision
     &  coef2real,coef2imag,boundreal,boundimag,
     &  diff0real,diff0imag,diff1real,diff1imag,
     &  slopereal,slopeimag
',`dnl
      $1
     &  coef2,bound
')dnl
      double precision deltax1
ifelse($1,`double complex',`dnl
cart_clinref_cell_body_2d(`_complex')dnl
',`dnl
cart_clinref_cell_body_2d(`')dnl
')dnl
')dnl
define(cart_clinref_op_edge_2d,`dnl
cart_clinrefine_op_subroutine_head_2d($2,$3)dnl
      $1
     &  arrayc(EDGE2d$2(cilo,cihi,0)),
     &  arrayf(EDGE2d$2(filo,fihi,0)),
     &  diff$2(cilo$2:cihi$2+1),
     &  slope$2(EDGE2d$2(cilo,cihi,0)),
     &  diff$3(cilo$3:cihi$3+1),
     &  slope$3(EDGE2d$2(cilo,cihi,0))
      integer ic0,ic1,ie0,ie1,if0,if1,ir0,ir1
      $1
     &  coef2,bound
      double precision deltax1
cart_clinref_edge_body_2d($2,$3)dnl
')dnl
define(cart_clinref_op_face_2d,`dnl
cart_clinrefine_op_subroutine_head_2d($2,$3)dnl
      $1
     &  arrayc(FACE2d$2(cilo,cihi,0)),
     &  arrayf(FACE2d$2(filo,fihi,0)),
     &  diff$2(cilo$2:cihi$2+1),
     &  slope$2(FACE2d$2(cilo,cihi,0)),
     &  diff$3(cilo$3:cihi$3+1),
     &  slope$3(FACE2d$2(cilo,cihi,0))
      integer ic$2,ic$3,ie$2,ie$3,if$2,if$3,ir$2,ir$3
      $1
     &  coef2,bound
      double precision deltax$3
cart_clinref_face_body_2d($2,$3)dnl
')dnl
define(cart_linref_op_node_2d,`dnl
cart_refine_op_subroutine_head_2d()dnl
      $1
     &  arrayc(NODE2d(cilo,cihi,0)),
     &  arrayf(NODE2d(filo,fihi,0))
      double precision x,y,realrat0,realrat1
      integer ic0,ic1,if0,if1,ie0,ie1,ir0,ir1,i,j
cart_linref_node_body_2d()dnl
')dnl
define(cart_clinref_op_side_2d,`dnl
cart_clinrefine_op_subroutine_head_2d($2,$3)dnl
      $1
     &  arrayc(SIDE2d$2(cilo,cihi,0)),
     &  arrayf(SIDE2d$2(filo,fihi,0)),
     &  diff$2(cilo$2:cihi$2+1),
     &  slope$2(SIDE2d$2(cilo,cihi,0)),
     &  diff$3(cilo$3:cihi$3+1),
     &  slope$3(SIDE2d$2(cilo,cihi,0))
      integer ic0,ic1,ie0,ie1,if0,if1,ir0,ir1
      $1
     &  coef2,bound
      double precision deltax1
cart_clinref_side_body_2d($2,$3)dnl
')dnl
define(cart_clinref_op_outerface_2d,`dnl
cart_clinrefine_op_subroutine_head_2d($2,$3)dnl
      $1
     &  arrayc(OUTERFACE2d$2(cilo,cihi,0)),
     &  arrayf(OUTERFACE2d$2(filo,fihi,0)),
     &  diff$2(cilo$2:cihi$2+1),
     &  slope$2(OUTERFACE2d$2(cilo,cihi,0)),
     &  diff$3(cilo$3:cihi$3+1),
     &  slope$3(OUTERFACE2d$2(cilo,cihi,0))
      integer ic$3,ie$3,if$3,ir$3
      $1
     &  coef2,bound
cart_clinref_outerface_body_2d($3)dnl
')dnl

c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for 3d Cartesian refine operators.
c
define(NDIM,3)dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl
include(FORTDIR/geom_m4cartopstuff.i)dnl
c
define(cart_refine_op_subroutine_head_3d,`dnl
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf)
c***********************************************************************
      implicit none
      double precision half,one
      parameter (half=0.5d0)
      parameter (one=1.0d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:NDIM-1)
      double precision
     &  dxc(0:NDIM-1),
     &  dxf(0:NDIM-1)
')dnl
c
define(cart_clinrefine_op_subroutine_head_3d,`dnl
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,dxc,dxf,
     &  arrayc,arrayf,
     &  diff$1,slope$1,diff$2,slope$2,diff$3,slope$3)
c***********************************************************************
      implicit none
      double precision zero,half,one,two
      parameter (zero=0.d0)
      parameter (half=0.5d0)
      parameter (one=1.d0)
      parameter (two=2.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2
      integer ratio(0:NDIM-1)
      double precision
     &  dxc(0:NDIM-1),
     &  dxf(0:NDIM-1),
     &  deltax(0:15,0:NDIM-1)
')dnl
c
define(cart_linref_cell_body_3d,`dnl
c
c***********************************************************************
c
coarse_fine_cell_deltas(0)dnl

coarse_fine_cell_deltas(1)dnl

coarse_fine_cell_deltas(2)dnl

      do if2=ifirstf2,ilastf2
coarsen_index(if2,ic2,ratio(2))dnl
         ir2=if2-ic2*ratio(2)
         z=deltax(ir2,2)/dxc(2)
         if( z .lt. 0.d0 ) then
           ic2 = ic2-1
           z = z + one
         endif
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
               arrayf(if0,if1,if2)=
     &          ( (arrayc(ic0,ic1,ic2)
     &      +(arrayc(ic0+1,ic1,ic2)-arrayc(ic0,ic1,ic2))*x)*(one-y)
     &      +(arrayc(ic0,ic1+1,ic2)
     &      +(arrayc(ic0+1,ic1+1,ic2)-arrayc(ic0,ic1+1,ic2))*x)*y )
     &       *(one-z)
     &         +( (arrayc(ic0,ic1,ic2+1)
     &      +(arrayc(ic0+1,ic1,ic2+1)-arrayc(ic0,ic1,ic2+1))*x)
     &       *(one-y)
     &      +(arrayc(ic0,ic1+1,ic2+1)
     &      +(arrayc(ic0+1,ic1+1,ic2+1)-arrayc(ic0,ic1+1,ic2+1))*x)*y)
     &       *z
          enddo
        enddo
      enddo
c
      return
      end
')dnl
c
define(cart_clinref_cell_body_3d,`dnl
c
c***********************************************************************
c
coarse_fine_cell_deltas(0)dnl

coarse_fine_cell_deltas(1)dnl

coarse_fine_cell_deltas(2)dnl

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
muscl_limited_cell_slopes$1(0,`ie0,ic1,ic2',`ie0-1,ic1,ic2',`ic0,ic1,ic2')dnl
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic0=ifirstc0,ilastc0
muscl_limited_cell_slopes$1(1,`ic0,ie1,ic2',`ic0,ie1-1,ic2',`ic0,ic1,ic2')dnl
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
muscl_limited_cell_slopes$1(2,`ic0,ic1,ie2',`ic0,ic1,ie2-1',`ic0,ic1,ic2')dnl
         enddo
      enddo

      do if2=ifirstf2,ilastf2
coarsen_index(if2,ic2,ratio(2))dnl
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)
         do if1=ifirstf1,ilastf1
coarsen_index(if1,ic1,ratio(1))dnl
            ir1=if1-ic1*ratio(1)
            deltax1=deltax(ir1,1)
            do if0=ifirstf0,ilastf0
coarsen_index(if0,ic0,ratio(0))dnl
               ir0=if0-ic0*ratio(0)
               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
     &                     +slope0(ic0,ic1,ic2)*deltax(ir0,0)
     &                     +slope1(ic0,ic1,ic2)*deltax1
     &                     +slope2(ic0,ic1,ic2)*deltax2
          enddo
        enddo
      enddo
c
      return
      end
')dnl
c
define(cart_clinref_edge_body_3d,`dnl
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

ifelse($1,`2',`
coarse_fine_cell_deltas(2)dnl
',`
coarse_fine_face_deltas(2)dnl
')dnl

ifelse($1,`2',`
      do ic2=ifirstc2,ilastc2
',`
      do ic2=ifirstc2,ilastc2+1
')dnl
ifelse($1,`1',`
         do ic1=ifirstc1,ilastc1
',`
         do ic1=ifirstc1,ilastc1+1
')dnl
ifelse($1,`0',`
muscl_limited_cell_slopes(0,`ie0,ic1,ic2',`ie0-1,ic1,ic2',`ic0,ic1,ic2')dnl
',`
muscl_limited_face_slopes(0,`ic0+1,ic1,ic2',`ic0,ic1,ic2',`ie0,ic1,ic2')dnl
')dnl
         enddo
      enddo

ifelse($1,`2',`
      do ic2=ifirstc2,ilastc2
',`
      do ic2=ifirstc2,ilastc2+1
')dnl
ifelse($1,`0',`
            do ic0=ifirstc0,ilastc0
',`
            do ic0=ifirstc0,ilastc0+1
')dnl

ifelse($1,`1',`
muscl_limited_cell_slopes(1,`ic0,ie1,ic2',`ic0,ie1-1,ic2',`ic0,ic1,ic2')dnl
',`
muscl_limited_face_slopes(1,`ic0,ic1+1,ic2',`ic0,ic1,ic2',`ic0,ie1,ic2')dnl
')dnl
         enddo
      enddo

ifelse($1,`1',`
         do ic1=ifirstc1,ilastc1
',`
         do ic1=ifirstc1,ilastc1+1
')dnl
ifelse($1,`0',`
            do ic0=ifirstc0,ilastc0
',`
            do ic0=ifirstc0,ilastc0+1
')dnl
ifelse($1,`2',`
muscl_limited_cell_slopes(2,`ic0,ic1,ie2',`ic0,ic1,ie2-1',`ic0,ic1,ic2')dnl
',`
muscl_limited_face_slopes(2,`ic0,ic1,ic2+1',`ic0,ic1,ic2',`ic0,ic1,ie2')dnl
')dnl
         enddo
      enddo

ifelse($1,`2',`
      do if2=ifirstf2,ilastf2
',`
      do if2=ifirstf2,ilastf2+1
')dnl
coarsen_index(if2,ic2,ratio(2))dnl
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)

ifelse($1,`1',`
         do if1=ifirstf1,ilastf1
',`
         do if1=ifirstf1,ilastf1+1
')dnl
coarsen_index(if1,ic1,ratio(1))dnl
            ir1=if1-ic1*ratio(1)
            deltax1=deltax(ir1,1)

ifelse($1,`0',`
            do if0=ifirstf0,ilastf0
',`
            do if0=ifirstf0,ilastf0+1
')dnl
coarsen_index(if0,ic0,ratio(0))dnl
               ir0=if0-ic0*ratio(0)
               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
     &                +slope0(ic0,ic1,ic2)*deltax(ir0,0)
     &                +slope1(ic0,ic1,ic2)*deltax1
     &                +slope2(ic0,ic1,ic2)*deltax2
          enddo
        enddo
      enddo
c
      return
      end
')dnl
c
define(cart_clinref_face_body_3d,`dnl
c
c***********************************************************************
c

coarse_fine_cell_deltas($3)dnl

coarse_fine_cell_deltas($2)dnl

coarse_fine_face_deltas($1)dnl

      do ic$3=ifirstc$3,ilastc$3
         do ic$2=ifirstc$2,ilastc$2
muscl_limited_face_slopes($1,`ic$1+1,ic$2,ic$3',`ic$1,ic$2,ic$3',`ie$1,ic$2,ic$3')dnl
         enddo
      enddo

      do ic$3=ifirstc$3,ilastc$3
         do ic$1=ifirstc$1,ilastc$1+1
muscl_limited_cell_slopes($2,`ic$1,ie$2,ic$3',`ic$1,ie$2-1,ic$3',`ic$1,ic$2,ic$3')dnl
         enddo
      enddo

      do ic$2=ifirstc$2,ilastc$2
         do ic$1=ifirstc$1,ilastc$1+1
muscl_limited_cell_slopes($3,`ic$1,ic$2,ie$3',`ic$1,ic$2,ie$3-1',`ic$1,ic$2,ic$3')dnl
         enddo
      enddo

      do if$3=ifirstf$3,ilastf$3
coarsen_index(if$3,ic$3,ratio($3))dnl
         ir$3=if$3-ic$3*ratio($3)
         deltax$3=deltax(ir$3,$3)
         do if$2=ifirstf$2,ilastf$2
coarsen_index(if$2,ic$2,ratio($2))dnl
            ir$2=if$2-ic$2*ratio($2)
            deltax$2=deltax(ir$2,$2)
            do if$1=ifirstf$1,ilastf$1+1
coarsen_index(if$1,ie$1,ratio($1))dnl
               ir$1=if$1-ie$1*ratio($1)
               arrayf(if$1,if$2,if$3)=arrayc(ie$1,ic$2,ic$3)
     &                +slope$1(ie$1,ic$2,ic$3)*deltax(ir$1,$1)
     &                +slope$2(ie$1,ic$2,ic$3)*deltax$2
     &                +slope$3(ie$1,ic$2,ic$3)*deltax$3
          enddo
        enddo
      enddo
c
      return
      end
')dnl

c
define(cart_linref_node_body_3d,`dnl
c
c***********************************************************************
c
      realrat0=one/dble(ratio(0))
      realrat1=one/dble(ratio(1))
      realrat2=one/dble(ratio(2))

node_refloop_pre(2)dnl
node_refloop_pre(1)dnl
node_refloop_pre(0)dnl
               x = dble(ir0)*realrat0
               y = dble(ir1)*realrat1
               z = dble(ir2)*realrat2
               arrayf(ie0,ie1,ie2)=
     &            ( (arrayc(ic0,ic1,ic2)*(one-x) +
     &               arrayc(ic0+1,ic1,ic2)*x)*(one-y)
     &            + (arrayc(ic0,ic1+1,ic2)*(one-x) +
     &               arrayc(ic0+1,ic1+1,ic2)*x)*y ) * (one-z) + 
     &            ( (arrayc(ic0,ic1,ic2+1)*(one-x) +
     &               arrayc(ic0+1,ic1,ic2+1)*x)*(one-y)
     &            + (arrayc(ic0,ic1+1,ic2+1)*(one-x) +
     &               arrayc(ic0+1,ic1+1,ic2+1)*x)*y ) * z
node_refloop_post()dnl 
node_refloop_post()dnl 
node_refloop_post()dnl 

      return
      end
')dnl
c
define(cart_clinref_outerface_body_3d,`dnl
c
c***********************************************************************
c

coarse_fine_cell_deltas($2)dnl

coarse_fine_cell_deltas($1)dnl

      do ic$2=ifirstc$2,ilastc$2
muscl_limited_cell_slopes($1,`ie$1,ic$2',`ie$1-1,ic$2',`ic$1,ic$2')dnl
      enddo

      do ic$1=ifirstc$1,ilastc$1
muscl_limited_cell_slopes($2,`ic$1,ie$2',`ic$1,ie$2-1',`ic$1,ic$2')dnl
      enddo

      do if$2=ifirstf$2,ilastf$2
coarsen_index(if$2,ic$2,ratio($2))dnl
         ir$2=if$2-ic$2*ratio($2)
         deltax$2=deltax(ir$2,$2)
         do if$1=ifirstf$1,ilastf$1
coarsen_index(if$1,ic$1,ratio($1))dnl
            ir$1=if$1-ic$1*ratio($1)
            arrayf(if$1,if$2)=arrayc(ic$1,ic$2)
     &            +slope$1(ic$1,ic$2)*deltax(ir$1,$1)
     &            +slope$2(ic$1,ic$2)*deltax$2
         enddo
      enddo
c
      return
      end
')dnl
c
define(cart_clinref_side_body_3d,`dnl
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

ifelse($1,`2',`
coarse_fine_face_deltas(2)dnl
',`
coarse_fine_cell_deltas(2)dnl
')dnl

ifelse($1,`2',`
      do ic2=ifirstc2,ilastc2+1
',`
      do ic2=ifirstc2,ilastc2
')dnl
ifelse($1,`1',`
         do ic1=ifirstc1,ilastc1+1
',`
         do ic1=ifirstc1,ilastc1
')dnl
ifelse($1,`0',`
muscl_limited_face_slopes(0,`ic0+1,ic1,ic2',`ic0,ic1,ic2',`ie0,ic1,ic2')dnl
',`
muscl_limited_cell_slopes(0,`ie0,ic1,ic2',`ie0-1,ic1,ic2',`ic0,ic1,ic2')dnl
')dnl
         enddo
      enddo

ifelse($1,`2',`
            do ic2=ifirstc2,ilastc2+1
',`
            do ic2=ifirstc2,ilastc2
')dnl
ifelse($1,`0',`
            do ic0=ifirstc0,ilastc0+1
',`
            do ic0=ifirstc0,ilastc0
')dnl
ifelse($1,`1',`
muscl_limited_face_slopes(1,`ic0,ic1+1,ic2',`ic0,ic1,ic2',`ic0,ie1,ic2')dnl
',`
muscl_limited_cell_slopes(1,`ic0,ie1,ic2',`ic0,ie1-1,ic2',`ic0,ic1,ic2')dnl
')dnl
         enddo
      enddo

ifelse($1,`1',`
            do ic1=ifirstc1,ilastc1+1
',`
            do ic1=ifirstc1,ilastc1
')dnl
ifelse($1,`0',`
            do ic0=ifirstc0,ilastc0+1
',`
            do ic0=ifirstc0,ilastc0
')dnl
ifelse($1,`2',`
muscl_limited_face_slopes(2,`ic0,ic1,ic2+1',`ic0,ic1,ic2',`ic0,ic1,ie2')dnl
',`
muscl_limited_cell_slopes(2,`ic0,ic1,ie2',`ic0,ic1,ie2-1',`ic0,ic1,ic2')dnl
')dnl
         enddo
      enddo

ifelse($1,`2',`
      do if2=ifirstf2,ilastf2+1
',`
      do if2=ifirstf2,ilastf2
')dnl
coarsen_index(if2,ic2,ratio(2))dnl
         ir2=if2-ic2*ratio(2)
         deltax2=deltax(ir2,2)

ifelse($1,`1',`
         do if1=ifirstf1,ilastf1+1
',`
         do if1=ifirstf1,ilastf1
')dnl
coarsen_index(if1,ic1,ratio(1))dnl
            ir1=if1-ic1*ratio(1)
            deltax1=deltax(ir1,1)

ifelse($1,`0',`
            do if0=ifirstf0,ilastf0+1
',`
            do if0=ifirstf0,ilastf0
')dnl
coarsen_index(if0,ic0,ratio(0))dnl
               ir0=if0-ic0*ratio(0)

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
     &                +slope0(ic0,ic1,ic2)*deltax(ir0,0)
     &                +slope1(ic0,ic1,ic2)*deltax1
     &                +slope2(ic0,ic1,ic2)*deltax2
          enddo
        enddo
      enddo
c
      return
      end
')dnl
c
define(cart_linref_op_cell_3d,`dnl
cart_refine_op_subroutine_head_3d()dnl
      $1
     &  arrayc(CELL3d(cilo,cihi,0)),
     &  arrayf(CELL3d(filo,fihi,0))
      double precision deltax(0:15,0:NDIM-1),x,y,z
      integer ic0,ic1,ic2,if0,if1,if2,ir0,ir1,ir2
cart_linref_cell_body_3d()dnl
')dnl
define(cart_clinref_op_cell_3d,`dnl
cart_clinrefine_op_subroutine_head_3d(0,1,2)dnl
      $1
     &  arrayc(CELL3d(cilo,cihi,0)),
     &  arrayf(CELL3d(filo,fihi,0)),
     &  diff0(cilo0:cihi0+1),
     &  slope0(CELL3d(cilo,cihi,0)),
     &  diff1(cilo1:cihi1+1),
     &  slope1(CELL3d(cilo,cihi,0)),
     &  diff2(cilo2:cihi2+1),
     &  slope2(CELL3d(cilo,cihi,0))
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,ir0,ir1,ir2
ifelse($1,`double complex',`dnl
      double precision
     &  coef2real,coef2imag,boundreal,boundimag,
     &  diff0real,diff0imag,diff1real,diff1imag,
     &  slopereal,slopeimag
',`dnl
      $1
     &  coef2,bound
')dnl
      double precision deltax1,deltax2
ifelse($1,`double complex',`dnl
cart_clinref_cell_body_3d(`_complex')dnl
',`dnl
cart_clinref_cell_body_3d(`')dnl
')dnl
')dnl
define(cart_clinref_op_edge_3d,`dnl
cart_clinrefine_op_subroutine_head_3d($2,$3,$4)dnl
      $1
     &  arrayc(EDGE3d$2(cilo,cihi,0)),
     &  arrayf(EDGE3d$2(filo,fihi,0)),
     &  diff$2(cilo$2:cihi$2+1),
     &  slope$2(EDGE3d$2(cilo,cihi,0)),
     &  diff$3(cilo$3:cihi$3+1),
     &  slope$3(EDGE3d$2(cilo,cihi,0)),
     &  diff$4(cilo$4:cihi$4+1),
     &  slope$4(EDGE3d$2(cilo,cihi,0))
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,
     &        ir0,ir1,ir2
      $1
     &  coef2,bound
      double precision deltax1,deltax2
cart_clinref_edge_body_3d($2,$3,$4)dnl
')dnl
define(cart_clinref_op_face_3d,`dnl
cart_clinrefine_op_subroutine_head_3d($2,$3,$4)dnl
      $1
     &  arrayc(FACE3d$2(cilo,cihi,0)),
     &  arrayf(FACE3d$2(filo,fihi,0)),
     &  diff$2(cilo$2:cihi$2+1),
     &  slope$2(FACE3d$2(cilo,cihi,0)),
     &  diff$3(cilo$3:cihi$3+1),
     &  slope$3(FACE3d$2(cilo,cihi,0)),
     &  diff$4(cilo$4:cihi$4+1),
     &  slope$4(FACE3d$2(cilo,cihi,0))
      integer ic$2,ic$3,ic$4,ie$2,ie$3,ie$4,if$2,if$3,if$4,
     &        ir$2,ir$3,ir$4
      $1
     &  coef2,bound
      double precision deltax$3,deltax$4
cart_clinref_face_body_3d($2,$3,$4)dnl
')dnl
define(cart_linref_op_node_3d,`dnl
cart_refine_op_subroutine_head_3d()dnl
      $1
     &  arrayc(NODE3d(cilo,cihi,0)),
     &  arrayf(NODE3d(filo,fihi,0))
      double precision x,y,z,realrat0,realrat1,realrat2
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,ir0,ir1,ir2,i,j,k
cart_linref_node_body_3d()dnl
')dnl
define(cart_clinref_op_outerface_3d,`dnl
cart_clinrefine_op_subroutine_head_3d($2,$3,$4)dnl
      $1
     &  arrayc(OUTERFACE3d$2(cilo,cihi,0)),
     &  arrayf(OUTERFACE3d$2(filo,fihi,0)),
     &  diff$2(cilo$2:cihi$2+1),
     &  slope$2(OUTERFACE3d$2(cilo,cihi,0)),
     &  diff$3(cilo$3:cihi$3+1),
     &  slope$3(OUTERFACE3d$2(cilo,cihi,0)),
     &  diff$4(cilo$4:cihi$4+1),
     &  slope$4(OUTERFACE3d$2(cilo,cihi,0))
      integer ic$3,ic$4,ie$3,ie$4,if$3,if$4,ir$3,ir$4
      $1
     &  coef2,bound
      double precision deltax$4
cart_clinref_outerface_body_3d($3,$4)dnl
')dnl
define(cart_clinref_op_side_3d,`dnl
cart_clinrefine_op_subroutine_head_3d($2,$3,$4)dnl
      $1
     &  arrayc(SIDE3d$2(cilo,cihi,0)),
     &  arrayf(SIDE3d$2(filo,fihi,0)),
     &  diff$2(cilo$2:cihi$2+1),
     &  slope$2(SIDE3d$2(cilo,cihi,0)),
     &  diff$3(cilo$3:cihi$3+1),
     &  slope$3(SIDE3d$2(cilo,cihi,0)),
     &  diff$4(cilo$4:cihi$4+1),
     &  slope$4(SIDE3d$2(cilo,cihi,0))
      integer ic0,ic1,ic2,ie0,ie1,ie2,if0,if1,if2,
     &        ir0,ir1,ir2
      $1
     &  coef2,bound
      double precision deltax2,deltax1
cart_clinref_side_body_3d($2,$3,$4)dnl
')dnl

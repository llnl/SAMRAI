c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for 3d Cartesian coarsen operators.
c
define(NDIM,3)dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl
define(cart_coarsen_op_subroutine_head_3d,`dnl
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2
      integer ratio(0:NDIM-1)
      double precision
     &  dxf(0:NDIM-1),
     &  dxc(0:NDIM-1)
')dnl
c
define(cart_wgtavg_cell_body_3d,`dnl
c
c***********************************************************************
c
      dVf = dxf(0)*dxf(1)*dxf(2)
      dVc = dxc(0)*dxc(1)*dxc(2)

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ic0=ifirstc0,ilastc0
ifelse($1,`double complex',`dnl
               arrayc(ic0,ic1,ic2)=cmplx(zero,zero)
',`dnl
               arrayc(ic0,ic1,ic2)=zero
')dnl
            enddo
         enddo
      enddo

      do ir2=0,ratio(2)-1
         do ir1=0,ratio(1)-1
            do ir0=0,ratio(0)-1
               do ic2=ifirstc2,ilastc2
                  if2=ic2*ratio(2)+ir2
                  do ic1=ifirstc1,ilastc1
                     if1=ic1*ratio(1)+ir1
                     do ic0=ifirstc0,ilastc0
                        if0=ic0*ratio(0)+ir0
                        arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                                      +arrayf(if0,if1,if2)*dVf
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic2=ifirstc2,ilastc2
         do ic1=ifirstc1,ilastc1
            do ic0=ifirstc0,ilastc0
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/dVc
            enddo
         enddo
      enddo
c
      return
      end
')dnl
c
define(cart_wgtavg_edge_body_3d,`dnl
c
c***********************************************************************
c
      lengthf=dxf($2)
      lengthc=dxc($2)

ifelse($2,`2',`
      do ic2=ifirstc2,ilastc2 
',`
      do ic2=ifirstc2,ilastc2+1
')dnl
ifelse($2,`1',`
         do ic1=ifirstc1,ilastc1
',`
         do ic1=ifirstc1,ilastc1+1
')dnl
ifelse($2,`0',`
            do ic0=ifirstc0,ilastc0
',`
            do ic0=ifirstc0,ilastc0+1
')dnl

ifelse($1,`double complex',`dnl
               arrayc(ic0,ic1,ic2)=cmplx(zero,zero)
',`dnl
               arrayc(ic0,ic1,ic2)=zero
')dnl
            enddo
         enddo
      enddo

ifelse($2,`2',`
      do ic2=ifirstc2,ilastc2 
         do ir=0,ratio($2)-1
            if2=ic2*ratio(2)+ir
',`
      do ic2=ifirstc2,ilastc2+1
            if2=ic2*ratio(2)
')dnl

ifelse($2,`1',`
         do ic1=ifirstc1,ilastc1
            do ir=0,ratio($2)-1
               if1=ic1*ratio(1)+ir
',`
         do ic1=ifirstc1,ilastc1+1
               if1=ic1*ratio(1)
')dnl

ifelse($2,`0',`
            do ic0=ifirstc0,ilastc0
               do ir=0,ratio($2)-1
                  if0=ic0*ratio(0)+ir
',`
            do ic0=ifirstc0,ilastc0+1
                  if0=ic0*ratio(0)
')dnl
                  arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                             +arrayf(if0,if1,if2)*lengthf
               enddo
            enddo
         enddo
      enddo

ifelse($2,`2',`
      do ic2=ifirstc2,ilastc2 
',`
      do ic2=ifirstc2,ilastc2+1
')dnl
ifelse($2,`1',`
         do ic1=ifirstc1,ilastc1
',`
         do ic1=ifirstc1,ilastc1+1
')dnl
ifelse($2,`0',`
            do ic0=ifirstc0,ilastc0
',`
            do ic0=ifirstc0,ilastc0+1
')dnl
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/lengthc
            enddo
         enddo
      enddo
c
      return
      end
')dnl
c
define(cart_wgtavg_face_body_3d,`dnl
c
c***********************************************************************
c
      areaf=dxf($3)*dxf($4)
      areac=dxc($3)*dxc($4)

      do ic$4=ifirstc$4,ilastc$4
         do ic$3=ifirstc$3,ilastc$3
            do ie$2=ifirstc$2,ilastc$2+1
ifelse($1,`double complex',`dnl
               arrayc(ie$2,ic$3,ic$4)=cmplx(zero,zero)
',`dnl
               arrayc(ie$2,ic$3,ic$4)=zero
')dnl
            enddo
         enddo
      enddo

      do ir$4=0,ratio($4)-1
         do ir$3=0,ratio($3)-1
            do ic$4=ifirstc$4,ilastc$4
               if$4=ic$4*ratio($4)+ir$4
               do ic$3=ifirstc$3,ilastc$3
                  if$3=ic$3*ratio($3)+ir$3
                  do ie$2=ifirstc$2,ilastc$2+1
                     if$2=ie$2*ratio($2)
                     arrayc(ie$2,ic$3,ic$4)=arrayc(ie$2,ic$3,ic$4)
     &                                +arrayf(if$2,if$3,if$4)*areaf
                  enddo
               enddo
            enddo
         enddo
      enddo

      do ic$4=ifirstc$4,ilastc$4
         do ic$3=ifirstc$3,ilastc$3
            do ie$2=ifirstc$2,ilastc$2+1
               arrayc(ie$2,ic$3,ic$4)=arrayc(ie$2,ic$3,ic$4)/areac
            enddo
         enddo
      enddo
c
      return
      end
')dnl
c
define(cart_wgtavg_side_body_3d,`dnl
c
c***********************************************************************
c
      areaf=dxf($3)*dxf($4)
      areac=dxc($3)*dxc($4)

ifelse($2,`2',`
      do ic2=ifirstc2,ilastc2+1
',`
      do ic2=ifirstc2,ilastc2 
')dnl
ifelse($2,`1',`
         do ic1=ifirstc1,ilastc1+1
',`
         do ic1=ifirstc1,ilastc1
')dnl
ifelse($2,`0',`
            do ic0=ifirstc0,ilastc0+1
',`
            do ic0=ifirstc0,ilastc0
')dnl

ifelse($1,`double complex',`dnl
               arrayc(ic0,ic1,ic2)=cmplx(zero,zero)
',`dnl
               arrayc(ic0,ic1,ic2)=zero
')dnl
            enddo
         enddo
      enddo

ifelse($2,`2',`
      do ic2=ifirstc2,ilastc2+1
               if2=ic2*ratio(2)
',`
      do ic2=ifirstc2,ilastc2 
         do ir2=0,ratio(2)-1
               if2=ic2*ratio(2)+ir2
')dnl
ifelse($2,`1',`
         do ic1=ifirstc1,ilastc1+1
                  if1=ic1*ratio(1)
',`
         do ic1=ifirstc1,ilastc1
            do ir1=0,ratio(1)-1
                  if1=ic1*ratio(1)+ir1
')dnl
ifelse($2,`0',`
            do ic0=ifirstc0,ilastc0+1
                  if0=ic0*ratio(0)
',`
            do ic0=ifirstc0,ilastc0
               do ir0=0,ratio(0)-1
                  if0=ic0*ratio(0)+ir0
')dnl
                     arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)
     &                                +arrayf(if0,if1,if2)*areaf
                  enddo
               enddo
            enddo
         enddo
      enddo

ifelse($2,`2',`
      do ic2=ifirstc2,ilastc2+1
',`
      do ic2=ifirstc2,ilastc2 
')dnl
ifelse($2,`1',`
         do ic1=ifirstc1,ilastc1+1
',`
         do ic1=ifirstc1,ilastc1
')dnl
ifelse($2,`0',`
            do ic0=ifirstc0,ilastc0+1
',`
            do ic0=ifirstc0,ilastc0
')dnl
               arrayc(ic0,ic1,ic2)=arrayc(ic0,ic1,ic2)/areac
            enddo
         enddo
      enddo
c
      return
      end
')dnl
c
define(cart_wgtavg_outerface_body_3d,`dnl
c
c***********************************************************************
c
      areaf=dxf($3)*dxf($4)
      areac=dxc($3)*dxc($4)

      do ic$4=ifirstc$4,ilastc$4
         do ic$3=ifirstc$3,ilastc$3
ifelse($1,`double complex',`dnl
            arrayc(ic$3,ic$4)=cmplx(zero,zero)
',`dnl
            arrayc(ic$3,ic$4)=zero
')dnl
         enddo
      enddo

      do ir$4=0,ratio($4)-1
         do ir$3=0,ratio($3)-1
            do ic$4=ifirstc$4,ilastc$4
               if$4=ic$4*ratio($4)+ir$4
               do ic$3=ifirstc$3,ilastc$3
                  if$3=ic$3*ratio($3)+ir$3
                  arrayc(ic$3,ic$4)=arrayc(ic$3,ic$4)
     &                              +arrayf(if$3,if$4)*areaf
               enddo
            enddo
         enddo
      enddo

      do ic$4=ifirstc$4,ilastc$4
         do ic$3=ifirstc$3,ilastc$3
            arrayc(ic$3,ic$4)=arrayc(ic$3,ic$4)/areac
         enddo
      enddo
c
      return
      end
')dnl
define(cart_wgtavg_outerside_body_3d,`dnl
c
c***********************************************************************
c
      areaf=dxf($3)*dxf($4)
      areac=dxc($3)*dxc($4)

ifelse($2,`0',`
      do ic_outer=ifirstc2,ilastc2
         do ic_inner=ifirstc1,ilastc1
',`')dnl
ifelse($2,`1',`
      do ic_outer=ifirstc2,ilastc2
         do ic_inner=ifirstc0,ilastc0
',`')dnl
ifelse($2,`2',`
      do ic_outer=ifirstc1,ilastc1
         do ic_inner=ifirstc0,ilastc0
',`')dnl

ifelse($1,`double complex',`dnl
            arrayc(ic_inner,ic_outer)=cmplx(zero,zero)
',`dnl
            arrayc(ic_inner,ic_outer)=zero
')dnl
         enddo
      enddo

ifelse($2,`0',`
      do ic_outer=ifirstc2,ilastc2
         do ir_outer=0,ratio(2)-1
            if_outer=ic_outer*ratio(2)+ir_outer
            do ic_inner=ifirstc1,ilastc1
               do ir_inner=0,ratio(1)-1
               if_inner=ic_inner*ratio(1)+ir_inner
',`')dnl
ifelse($2,`1',`
      do ic_outer=ifirstc2,ilastc2
         do ir_outer=0,ratio(2)-1
            if_outer=ic_outer*ratio(2)+ir_outer
            do ic_inner=ifirstc0,ilastc0
               do ir_inner=0,ratio(0)-1
               if_inner=ic_inner*ratio(0)+ir_inner
',`')dnl
ifelse($2,`2',`
      do ic_outer=ifirstc1,ilastc1
         do ir_outer=0,ratio(1)-1
            if_outer=ic_outer*ratio(1)+ir_outer
            do ic_inner=ifirstc0,ilastc0
               do ir_inner=0,ratio(0)-1
               if_inner=ic_inner*ratio(0)+ir_inner
',`')dnl
 
                  arrayc(ic_inner,ic_outer)=
     &               arrayc(ic_inner,ic_outer)
     &               +arrayf(if_inner,if_outer)*areaf
               enddo
            enddo
         enddo
      enddo

ifelse($2,`0',`
      do ic_outer=ifirstc2,ilastc2
         do ic_inner=ifirstc1,ilastc1
',`')dnl
ifelse($2,`1',`
      do ic_outer=ifirstc2,ilastc2
         do ic_inner=ifirstc0,ilastc0
',`')dnl
ifelse($2,`2',`
      do ic_outer=ifirstc1,ilastc1
         do ic_inner=ifirstc0,ilastc0
',`')dnl
            arrayc(ic_inner,ic_outer)=
     &         arrayc(ic_inner,ic_outer)/areac
         enddo
      enddo
c
      return
      end
')dnl
c
define(cart_wgtavg_op_cell_3d,`dnl
cart_coarsen_op_subroutine_head_3d()dnl
      $1
     &  arrayf(CELL3d(filo,fihi,0)),
     &  arrayc(CELL3d(cilo,cihi,0))
      double precision dVf,dVc
      integer ic0,ic1,ic2,if0,if1,if2,ir0,ir1,ir2
cart_wgtavg_cell_body_3d($1)dnl
')dnl
c
define(cart_wgtavg_op_edge_3d,`dnl
cart_coarsen_op_subroutine_head_3d()dnl
      $1
     &  arrayf(EDGE3d$2(filo,fihi,0)),
     &  arrayc(EDGE3d$2(cilo,cihi,0))
      double precision lengthf,lengthc
      integer ic0,ic1,ic2,if0,if1,if2,ir
cart_wgtavg_edge_body_3d($1,$2,$3,$4)dnl
')dnl
define(cart_wgtavg_op_face_3d,`dnl
cart_coarsen_op_subroutine_head_3d()dnl
      $1
     &  arrayf(FACE3d$2(filo,fihi,0)),
     &  arrayc(FACE3d$2(cilo,cihi,0))
      double precision areaf,areac
      integer ie$2,ic$3,ic$4,if$2,if$3,if$4,ir$3,ir$4
cart_wgtavg_face_body_3d($1,$2,$3,$4)dnl
')dnl
define(cart_wgtavg_op_side_3d,`dnl
cart_coarsen_op_subroutine_head_3d()dnl
      $1
     &  arrayf(SIDE3d$2(filo,fihi,0)),
     &  arrayc(SIDE3d$2(cilo,cihi,0))
      double precision areaf,areac
      integer ic0,ic1,ic2,if0,if1,if2,
ifelse($2,`0',`
     &  ir1,ir2
',`')dnl
ifelse($2,`1',`
     &  ir0,ir2
',`')dnl
ifelse($2,`2',`
     &  ir0,ir1
',`')dnl
cart_wgtavg_side_body_3d($1,$2,$3,$4)dnl
')dnl
define(cart_wgtavg_op_outerface_3d,`dnl
cart_coarsen_op_subroutine_head_3d()dnl
      $1
     &  arrayf(OUTERFACE3d$2(filo,fihi,0)),
     &  arrayc(OUTERFACE3d$2(cilo,cihi,0))
      double precision areaf,areac
      integer ic$3,ic$4,if$3,if$4,ir$3,ir$4
cart_wgtavg_outerface_body_3d($1,$2,$3,$4)dnl
')dnl
define(cart_wgtavg_op_outerside_3d,`dnl
cart_coarsen_op_subroutine_head_3d()dnl
      $1
     &  arrayf(OUTERSIDE3d$2(filo,fihi,0)),
     &  arrayc(OUTERSIDE3d$2(cilo,cihi,0))
      double precision areaf,areac
      integer ic_outer,ic_inner,if_outer,if_inner,
     &  ir_outer,ir_inner
cart_wgtavg_outerside_body_3d($1,$2,$3,$4)dnl
')dnl

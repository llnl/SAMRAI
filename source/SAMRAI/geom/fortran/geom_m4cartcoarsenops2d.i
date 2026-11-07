c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for 2d Cartesian coarsen operators.
c
define(NDIM,2)dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl
define(cart_coarsen_op_subroutine_head_2d,`dnl
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1
      integer ratio(0:NDIM-1)
      double precision
     &  dxf(0:NDIM-1),
     &  dxc(0:NDIM-1)
')dnl
c
define(cart_wgtavg_cell_body_2d,`dnl
c
c***********************************************************************
c
      dVf = dxf(0)*dxf(1)
      dVc = dxc(0)*dxc(1)

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
ifelse($1,`double complex',`dnl
            arrayc(ic0,ic1)=cmplx(zero,zero)
',`dnl
            arrayc(ic0,ic1)=zero
')dnl
         enddo
      enddo

      do ir1=0,ratio(1)-1
         do ir0=0,ratio(0)-1
            do ic1=ifirstc1,ilastc1
               if1=ic1*ratio(1)+ir1
               do ic0=ifirstc0,ilastc0
                  if0=ic0*ratio(0)+ir0
                  arrayc(ic0,ic1)=arrayc(ic0,ic1)
     &                            +arrayf(if0,if1)*dVf
               enddo
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1
         do ic0=ifirstc0,ilastc0
            arrayc(ic0,ic1)=arrayc(ic0,ic1)/dVc
        enddo
      enddo
c
      return
      end
')dnl
c
define(cart_wgtavg_edge_body_2d,`dnl
c
c***********************************************************************
c
      lengthf=dxf($2)
      lengthc=dxc($2)

      do ic1=ifirstc1,ilastc1+$3
         do ic0=ifirstc0,ilastc0+$2
ifelse($1,`double complex',`dnl
            arrayc(ic0,ic1)=cmplx(zero,zero)
',`dnl
            arrayc(ic0,ic1)=zero
')dnl
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+$3
ifelse($2,`1',`
         do ir=0,ratio(1)-1
            if1=ic1*ratio(1)+ir
',`
            if1=ic1*ratio(1)
')dnl
            do ic0=ifirstc0,ilastc0+$2
ifelse($2,`0',`
               do ir=0,ratio(0)-1
                  if0=ic0*ratio(0)+ir
',`
                  if0=ic0*ratio(0)
')dnl
                  arrayc(ic0,ic1)=arrayc(ic0,ic1)
     &                           +arrayf(if0,if1)*lengthf
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+$3
         do ic0=ifirstc0,ilastc0+$2
            arrayc(ic0,ic1)=arrayc(ic0,ic1)/lengthc
        enddo
      enddo
c
      return
      end
')dnl
c
define(cart_wgtavg_face_body_2d,`dnl
c
c***********************************************************************
c
      lengthf=dxf($3)
      lengthc=dxc($3)

      do ic$3=ifirstc$3,ilastc$3
         do ie$2=ifirstc$2,ilastc$2+1
ifelse($1,`double complex',`dnl
            arrayc(ie$2,ic$3)=cmplx(zero,zero)
',`dnl
            arrayc(ie$2,ic$3)=zero
')dnl
         enddo
      enddo

      do ir$3=0,ratio($3)-1
         do ic$3=ifirstc$3,ilastc$3
            if$3=ic$3*ratio($3)+ir$3
            do ie$2=ifirstc$2,ilastc$2+1
               if$2=ie$2*ratio($2)
               arrayc(ie$2,ic$3)=arrayc(ie$2,ic$3)
     &                           +arrayf(if$2,if$3)*lengthf
            enddo
         enddo
      enddo

      do ic$3=ifirstc$3,ilastc$3
         do ie$2=ifirstc$2,ilastc$2+1
            arrayc(ie$2,ic$3)=arrayc(ie$2,ic$3)/lengthc
        enddo
      enddo
c
      return
      end
')dnl
c
define(cart_wgtavg_side_body_2d,`dnl
c
c***********************************************************************
c
      lengthf=dxf($3)
      lengthc=dxc($3)

      do ic1=ifirstc1,ilastc1+$2
         do ic0=ifirstc0,ilastc0+$3
ifelse($1,`double complex',`dnl
            arrayc(ic0,ic1)=cmplx(zero,zero)
',`dnl
            arrayc(ic0,ic1)=zero
')dnl
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+$2
ifelse($2,`1',`
            if1=ic1*ratio(1)
',`
         do ir=0,ratio(1)-1
            if1=ic1*ratio(1)+ir
')dnl
            do ic0=ifirstc0,ilastc0+$3
ifelse($2,`0',`
                  if0=ic0*ratio(0)
',`
               do ir=0,ratio(0)-1
                  if0=ic0*ratio(0)+ir
')dnl
               arrayc(ic0,ic1)=arrayc(ic0,ic1)
     &                           +arrayf(if0,if1)*lengthf
            enddo
         enddo
      enddo

      do ic1=ifirstc1,ilastc1+$2
         do ic0=ifirstc0,ilastc0+$3
            arrayc(ic0,ic1)=arrayc(ic0,ic1)/lengthc
        enddo
      enddo
c
      return
      end
')dnl
c
define(cart_wgtavg_outerface_body_2d,`dnl
c
c***********************************************************************
c
      lengthf=dxf($3)
      lengthc=dxc($3)

      do ic$3=ifirstc$3,ilastc$3
ifelse($1,`double complex',`dnl
         arrayc(ic$3)=cmplx(zero,zero)
',`dnl
         arrayc(ic$3)=zero
')dnl
      enddo

      do ir$3=0,ratio($3)-1
         do ic$3=ifirstc$3,ilastc$3
            if$3=ic$3*ratio($3)+ir$3
            arrayc(ic$3)=arrayc(ic$3)+arrayf(if$3)*lengthf
         enddo
      enddo

      do ic$3=ifirstc$3,ilastc$3
         arrayc(ic$3)=arrayc(ic$3)/lengthc
      enddo
c
      return
      end
')dnl
c
define(cart_wgtavg_outerside_body_2d,`dnl
c
c***********************************************************************
c
      lengthf=dxf($3)
      lengthc=dxc($3)

      do ic$3=ifirstc$3,ilastc$3
ifelse($1,`double complex',`dnl
         arrayc(ic$3)=cmplx(zero,zero)
',`dnl
         arrayc(ic$3)=zero
')dnl
      enddo

      do ir$3=0,ratio($3)-1
         do ic$3=ifirstc$3,ilastc$3
            if$3=ic$3*ratio($3)+ir$3
            arrayc(ic$3)=arrayc(ic$3)+arrayf(if$3)*lengthf
         enddo
      enddo

      do ic$3=ifirstc$3,ilastc$3
         arrayc(ic$3)=arrayc(ic$3)/lengthc
      enddo
c
      return
      end
')dnl
c
define(cart_wgtavg_op_cell_2d,`dnl
cart_coarsen_op_subroutine_head_2d()dnl
      $1
     &  arrayf(CELL2d(filo,fihi,0)),
     &  arrayc(CELL2d(cilo,cihi,0))
      double precision dVf,dVc
      integer ic0,ic1,if0,if1,ir0,ir1
cart_wgtavg_cell_body_2d($1)dnl
')dnl
c
define(cart_wgtavg_op_edge_2d,`dnl
cart_coarsen_op_subroutine_head_2d()dnl
      $1
     &  arrayf(EDGE2d$2(filo,fihi,0)),
     &  arrayc(EDGE2d$2(cilo,cihi,0))
      double precision lengthf, lengthc
      integer ic0,ic1,if0,if1,ir
cart_wgtavg_edge_body_2d($1,$2,$3)dnl
')dnl
define(cart_wgtavg_op_face_2d,`dnl
cart_coarsen_op_subroutine_head_2d()dnl
      $1
     &  arrayf(FACE2d$2(filo,fihi,0)),
     &  arrayc(FACE2d$2(cilo,cihi,0))
      double precision lengthf, lengthc
      integer ie$2,ic$3,if$2,if$3,ir$3
cart_wgtavg_face_body_2d($1,$2,$3)dnl
')dnl
define(cart_wgtavg_op_side_2d,`dnl
cart_coarsen_op_subroutine_head_2d()dnl
      $1
     &  arrayf(SIDE2d$2(filo,fihi,0)),
     &  arrayc(SIDE2d$2(cilo,cihi,0))
      double precision lengthf, lengthc
      integer ic0,ic1,if0,if1,ir
cart_wgtavg_side_body_2d($1,$2,$3)dnl
')dnl
define(cart_wgtavg_op_outerface_2d,`dnl
cart_coarsen_op_subroutine_head_2d()dnl
      $1
     &  arrayf(OUTERFACE2d$2(filo,fihi,0)),
     &  arrayc(OUTERFACE2d$2(cilo,cihi,0))
      double precision lengthf, lengthc
      integer ic$3,if$3,ir$3
cart_wgtavg_outerface_body_2d($1,$2,$3)dnl
')dnl
define(cart_wgtavg_op_outerside_2d,`dnl
cart_coarsen_op_subroutine_head_2d()dnl
      $1
     &  arrayf(OUTERSIDE2d$2(filo,fihi,0)),
     &  arrayc(OUTERSIDE2d$2(cilo,cihi,0))
      double precision lengthf, lengthc
      integer ic$3,if$3,ir$3
cart_wgtavg_outerside_body_2d($1,$2,$3)dnl
')dnl

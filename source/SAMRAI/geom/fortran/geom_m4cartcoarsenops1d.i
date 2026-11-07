c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for 1d Cartesian coarsen operators.
c
define(NDIM,1)dnl
include(PDAT_FORTDIR/pdat_m4arrdim1d.i)dnl
define(cart_coarsen_op_subroutine_head_1d,`dnl
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0,
     &  ratio,dxf,dxc,
     &  arrayf,arrayc)
c***********************************************************************
      implicit none
      double precision zero
      parameter (zero=0.d0)
c
      integer
     &  ifirstc0,ilastc0,
     &  filo0,fihi0,
     &  cilo0,cihi0
      integer ratio(0:NDIM-1)
      double precision
     &  dxf(0:NDIM-1),
     &  dxc(0:NDIM-1)
')dnl
c
define(cart_wgtavg_cell_body_1d,`dnl
c
c***********************************************************************
c
      dVf = dxf(0)
      dVc = dxc(0)

      do ic0=ifirstc0,ilastc0
ifelse($1,`double complex',`dnl
         arrayc(ic0) = cmplx(zero, zero)
',`dnl
         arrayc(ic0) = zero
')dnl
      enddo

      do ir0=0,ratio(0)-1
         do ic0=ifirstc0,ilastc0
            if0=ic0*ratio(0)+ir0
            arrayc(ic0)=arrayc(ic0)+arrayf(if0)*dVf
         enddo
      enddo

      do ic0=ifirstc0,ilastc0
         arrayc(ic0)=arrayc(ic0)/dVc
      enddo
c
      return
      end
')dnl
c
define(cart_injection_body_1d,`dnl
c
c***********************************************************************
c
      do $1=$2
         arrayc($1)=arrayf($1*ratio(0))
      enddo
c
      return
      end
')dnl
c
define(cart_wgtavg_op_cell_1d,`dnl
cart_coarsen_op_subroutine_head_1d()dnl
      $1
     &  arrayf(CELL1d(filo,fihi,0)),
     &  arrayc(CELL1d(cilo,cihi,0))
      double precision dVf,dVc
      integer ic0,if0,ir0
cart_wgtavg_cell_body_1d($1)dnl
')dnl
c
define(cart_wgtavg_op_edge_1d,`dnl
cart_coarsen_op_subroutine_head_1d()dnl
      $1
     &  arrayf(EDGE1d(filo,fihi,0)),
     &  arrayc(EDGE1d(cilo,cihi,0))
      double precision dVf,dVc
      integer ic0,if0,ir0
cart_wgtavg_cell_body_1d($1)dnl
')dnl
c
define(cart_wgtavg_op_face_1d,`dnl
cart_coarsen_op_subroutine_head_1d()dnl
      $1
     &  arrayf(FACE1d(filo,fihi,0)),
     &  arrayc(FACE1d(cilo,cihi,0))
      integer ie0
cart_injection_body_1d(`ie0',`ifirstc0,ilastc0+1')dnl
')dnl
c
define(cart_wgtavg_op_outerface_1d,`dnl
cart_coarsen_op_subroutine_head_1d()dnl
      $1
     &  arrayf(OUTERFACE1d(filo,fihi,0)),
     &  arrayc(OUTERFACE1d(cilo,cihi,0))
c
c***********************************************************************
c
      arrayc(1)=arrayf(1)
c
      return
      end    
')dnl
define(cart_wgtavg_op_outerside_1d,`dnl
cart_coarsen_op_subroutine_head_1d()dnl
      $1
     &  arrayf(OUTERSIDE1d(filo,fihi,0)),
     &  arrayc(OUTERSIDE1d(cilo,cihi,0))
c
c***********************************************************************
c
      arrayc(1)=arrayf(1)
c
      return
      end    
')dnl
define(cart_wgtavg_op_side_1d,`dnl
cart_coarsen_op_subroutine_head_1d()dnl
      $1
     &  arrayf(SIDE1d(filo,fihi,0)),
     &  arrayc(SIDE1d(cilo,cihi,0))
      integer ie0
cart_injection_body_1d(`ie0',`ifirstc0,ilastc0+1')dnl
')dnl
c

c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for 3d constant coarsen operators.
c
define(NDIM,3)dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl
define(con_coarsen_op_subroutine_head_3d,`dnl
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  ratio,
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
')dnl
c
define(conavg_node_body_3d,`dnl
c
c***********************************************************************
c
      do ie2=ifirstc2,ilastc2+1
         if2=ie2*ratio(2)
         do ie1=ifirstc1,ilastc1+1
            if1=ie1*ratio(1)
            do ie0=ifirstc0,ilastc0+1
               arrayc(ie0,ie1,ie2)=arrayf(ie0*ratio(0),if1,if2)
            enddo
         enddo
      enddo
c
      return
      end
')dnl
c
define(conavg_op_node_3d,`dnl
con_coarsen_op_subroutine_head_3d()dnl
      $1
     &  arrayf(NODE3d(filo,fihi,0)),
     &  arrayc(NODE3d(cilo,cihi,0))
      integer ie0,ie1,ie2,if1,if2
conavg_node_body_3d()dnl
')dnl

define(conavg_op_outernode_3d0,`dnl
con_coarsen_op_subroutine_head_3d()dnl
      $1
     &  arrayf(OUTERNODE3d0(filo,fihi,0)),
     &  arrayc(OUTERNODE3d0(cilo,cihi,0))
      integer ic1,ic2,if1,if2
c
c***********************************************************************
c
      do ic2=ifirstc2+1,ilastc2
         if2=ic2*ratio(2)
         do ic1=ifirstc1+1,ilastc1
            if1=ic1*ratio(1)
            arrayc(ic1,ic2)=arrayf(if1,if2)
         enddo
      enddo

      return
      end
')dnl
c

define(conavg_op_outernode_3d1,`dnl 
con_coarsen_op_subroutine_head_3d()dnl 
      $1
     &  arrayf(OUTERNODE3d1(filo,fihi,0)),
     &  arrayc(OUTERNODE3d1(cilo,cihi,0))
      integer ic0,ic2,if0,if2
c
c***********************************************************************
c
      do ic2=ifirstc2+1,ilastc2
         if2=ic2*ratio(2)
         do ic0=ifirstc0,ilastc0+1
            if0=ic0*ratio(0)
            arrayc(ic0,ic2)=arrayf(if0,if2)
         enddo
      enddo

      return
      end
')dnl

define(conavg_op_outernode_3d2,`dnl 
con_coarsen_op_subroutine_head_3d()dnl 
      $1
     &  arrayf(OUTERNODE3d2(filo,fihi,0)),
     &  arrayc(OUTERNODE3d2(cilo,cihi,0))
      integer ic0,ic1,if0,if1
c
c***********************************************************************
c
      do ic1=ifirstc1,ilastc1+1
         if1=ic1*ratio(1)
         do ic0=ifirstc0,ilastc0+1
            if0=ic0*ratio(0)
            arrayc(ic0,ic1)=arrayf(if0,if1)
         enddo
      enddo

      return
      end
')dnl

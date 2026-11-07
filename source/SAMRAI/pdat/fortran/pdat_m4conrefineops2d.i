c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for 2d constant refine operators.
c
define(NDIM,2)dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl
include(PDAT_FORTDIR/pdat_m4conopstuff.i)dnl
c
define(con_refine_op_subroutine_head_2d,`dnl
     &  ifirstc0,ifirstc1,ilastc0,ilastc1,
     &  ifirstf0,ifirstf1,ilastf0,ilastf1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  filo0,filo1,fihi0,fihi1,
     &  ratio,
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
')dnl
c
define(conref_cell_body_2d,`dnl
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1
coarsen_index(if1,ic1,ratio(1))dnl
         do if0=ifirstf0,ilastf0
coarsen_index(if0,ic0,ratio(0))dnl
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
')dnl
c
define(conref_edge_body_2d,`dnl
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1+$2
coarsen_index(if1,ic1,ratio(1))dnl
         do if0=ifirstf0,ilastf0+$1
coarsen_index(if0,ic0,ratio(0))dnl
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
')dnl
c
define(conref_face_body_2d,`dnl
c
c***********************************************************************
c
      do if$2=ifirstf$2,ilastf$2
coarsen_index(if$2,ic$2,ratio($2))dnl
         do if$1=ifirstf$1,ilastf$1+1
coarsen_index(if$1,ie$1,ratio($1))dnl
            arrayf(if$1,if$2)=arrayc(ie$1,ic$2)
         enddo
      enddo
c
      return
      end
')dnl
c
define(conref_side_body_2d,`dnl
c
c***********************************************************************
c
      do if1=ifirstf1,ilastf1+$1
coarsen_index(if1,ic1,ratio(1))dnl
         do if0=ifirstf0,ilastf0+$2
coarsen_index(if0,ic0,ratio(0))dnl
            arrayf(if0,if1)=arrayc(ic0,ic1)
         enddo
      enddo
c
      return
      end
')dnl
c
define(conref_outerface_body_2d,`dnl
c
c***********************************************************************
c
      do if$2=ifirstf$2,ilastf$2
coarsen_index(if$2,ic$2,ratio($2))dnl
         arrayf(if$2)=arrayc(ic$2)
      enddo
c
      return
      end
')dnl
c
define(conref_op_cell_2d,`dnl
con_refine_op_subroutine_head_2d()dnl
      $1
     &  arrayc(CELL2d(cilo,cihi,0)),
     &  arrayf(CELL2d(filo,fihi,0))
      integer ic0,ic1,if0,if1
conref_cell_body_2d()dnl
')dnl
define(conref_op_edge_2d,`dnl
con_refine_op_subroutine_head_2d()dnl
      $1
     &  arrayc(EDGE2d$2(cilo,cihi,0)),
     &  arrayf(EDGE2d$2(filo,fihi,0))
      integer ic0,ic1,if0,if1
conref_edge_body_2d($2,$3)dnl
')dnl
define(conref_op_face_2d,`dnl
con_refine_op_subroutine_head_2d()dnl
      $1
     &  arrayc(FACE2d$2(cilo,cihi,0)),
     &  arrayf(FACE2d$2(filo,fihi,0))
      integer ie$2,ic$3,if$2,if$3
conref_face_body_2d($2,$3)dnl
')dnl
define(conref_op_side_2d,`dnl
con_refine_op_subroutine_head_2d()dnl
      $1
     &  arrayc(SIDE2d$2(cilo,cihi,0)),
     &  arrayf(SIDE2d$2(filo,fihi,0))
      integer ic0,ic1,if0,if1
conref_side_body_2d($2,$3)dnl
')dnl
define(conref_op_outerface_2d,`dnl
con_refine_op_subroutine_head_2d()dnl
      $1
     &  arrayc(OUTERFACE2d$2(cilo,cihi,0)),
     &  arrayf(OUTERFACE2d$2(filo,fihi,0))
      integer ic$3,if$3
conref_outerface_body_2d($2,$3)dnl
')dnl

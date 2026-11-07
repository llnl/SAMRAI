c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for 1d constant refine operators.
c
define(NDIM,1)dnl
include(PDAT_FORTDIR/pdat_m4arrdim1d.i)dnl
include(PDAT_FORTDIR/pdat_m4conopstuff.i)dnl
c
define(con_refine_op_subroutine_head_1d,`dnl
     &  ifirstc0,ilastc0,
     &  ifirstf0,ilastf0,
     &  cilo0,cihi0,
     &  filo0,fihi0,
     &  ratio,
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
')dnl
c
define(conref_cell_body_1d,`dnl
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
coarsen_index(if0,ic0,ratio(0))dnl
         arrayf(if0)=arrayc(ic0)
      enddo
c
      return
      end
')dnl
c
define(conref_edge_body_1d,`dnl
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0
coarsen_index(if0,ie0,ratio(0))dnl
         arrayf(if0)=arrayc(ie0)
      enddo
c
      return
      end
')dnl
c
define(conref_face_body_1d,`dnl
c
c***********************************************************************
c
      do if0=ifirstf0,ilastf0+1
coarsen_face_index(if0,ie0,ratio(0))dnl
         arrayf(if0)=arrayc(ie0)
      enddo
c
      return
      end
')dnl
c
define(conref_op_cell_1d,`dnl
con_refine_op_subroutine_head_1d()dnl
      $1
     &  arrayc(CELL1d(cilo,cihi,0)),
     &  arrayf(CELL1d(filo,fihi,0))
      integer ic0,if0
conref_cell_body_1d()dnl
')dnl
define(conref_op_edge_1d,`dnl
con_refine_op_subroutine_head_1d()dnl
      $1
     &  arrayc(EDGE1d(cilo,cihi,0)),
     &  arrayf(EDGE1d(filo,fihi,0))
      integer ie0,if0
conref_edge_body_1d()dnl
')dnl
define(conref_op_face_1d,`dnl
con_refine_op_subroutine_head_1d()dnl
      $1
     &  arrayc(FACE1d(cilo,cihi,0)),
     &  arrayf(FACE1d(filo,fihi,0))
      integer ie0,if0,it
conref_face_body_1d()dnl
')dnl
define(conref_op_side_1d,`dnl
con_refine_op_subroutine_head_1d()dnl
      $1
     &  arrayc(SIDE1d(cilo,cihi,0)),
     &  arrayf(SIDE1d(filo,fihi,0))
      integer ie0,if0,it
conref_face_body_1d()dnl
')dnl
define(conref_op_outerface_1d,`dnl
con_refine_op_subroutine_head_1d()dnl
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

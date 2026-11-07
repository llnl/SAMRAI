c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for 3d constant refine operators.
c
define(NDIM,3)dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl
include(PDAT_FORTDIR/pdat_m4conopstuff.i)dnl
c
define(con_refine_op_subroutine_head_3d,`dnl
     &  ifirstc0,ifirstc1,ifirstc2,ilastc0,ilastc1,ilastc2,
     &  ifirstf0,ifirstf1,ifirstf2,ilastf0,ilastf1,ilastf2,
     &  cilo0,cilo1,cilo2,cihi0,cihi1,cihi2,
     &  filo0,filo1,filo2,fihi0,fihi1,fihi2,
     &  ratio,
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
')dnl
c
define(conref_cell_body_3d,`dnl
c
c***********************************************************************
c
      do if2=ifirstf2,ilastf2
coarsen_index(if2,ic2,ratio(2))dnl
         do if1=ifirstf1,ilastf1
coarsen_index(if1,ic1,ratio(1))dnl
            do if0=ifirstf0,ilastf0
coarsen_index(if0,ic0,ratio(0))dnl
               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)
          enddo
        enddo
      enddo
c
      return
      end
')dnl
c
define(conref_edge_body_3d,`dnl
c
c***********************************************************************
c

ifelse($1,`2',`
      do if2=ifirstf2,ilastf2
',`
      do if2=ifirstf2,ilastf2+1
')dnl
coarsen_index(if2,ic2,ratio(2))dnl

ifelse($1,`1',`
         do if1=ifirstf1,ilastf1
',`
         do if1=ifirstf1,ilastf1+1
')dnl
coarsen_index(if1,ic1,ratio(1))dnl

ifelse($1,`0',`
            do if0=ifirstf0,ilastf0
',`
            do if0=ifirstf0,ilastf0+1
')dnl
coarsen_index(if0,ic0,ratio(0))dnl

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
')dnl
c
define(conref_face_body_3d,`dnl
c
c***********************************************************************
c
      do if$3=ifirstf$3,ilastf$3
coarsen_index(if$3,ic$3,ratio($3))dnl
         do if$2=ifirstf$2,ilastf$2
coarsen_index(if$2,ic$2,ratio($2))dnl
            do if$1=ifirstf$1,ilastf$1+1
coarsen_index(if$1,ie$1,ratio($1))dnl
               arrayf(if$1,if$2,if$3)=arrayc(ie$1,ic$2,ic$3)
            enddo
         enddo
      enddo
c
      return
      end
')dnl
c
define(conref_outerface_body_3d,`dnl
c
c***********************************************************************
c
      do if$3=ifirstf$3,ilastf$3
coarsen_index(if$3,ic$3,ratio($3))dnl
         do if$2=ifirstf$2,ilastf$2
coarsen_index(if$2,ic$2,ratio($2))dnl
            arrayf(if$2,if$3)=arrayc(ic$2,ic$3)
         enddo
      enddo
c
      return
      end
')dnl
c
define(conref_side_body_3d,`dnl
c
c***********************************************************************
c

ifelse($1,`2',`
      do if2=ifirstf2,ilastf2+1
',`
      do if2=ifirstf2,ilastf2
')dnl
coarsen_index(if2,ic2,ratio(2))dnl

ifelse($1,`1',`
         do if1=ifirstf1,ilastf1+1
',`
         do if1=ifirstf1,ilastf1
')dnl
coarsen_index(if1,ic1,ratio(1))dnl

ifelse($1,`0',`
            do if0=ifirstf0,ilastf0+1
',`
            do if0=ifirstf0,ilastf0
')dnl
coarsen_index(if0,ic0,ratio(0))dnl

               arrayf(if0,if1,if2)=arrayc(ic0,ic1,ic2)

            enddo
         enddo
      enddo
c
      return
      end
')dnl
c
define(conref_op_cell_3d,`dnl
con_refine_op_subroutine_head_3d()dnl
      $1
     &  arrayc(CELL3d(cilo,cihi,0)),
     &  arrayf(CELL3d(filo,fihi,0))
      integer ic0,ic1,ic2,if0,if1,if2
conref_cell_body_3d()dnl
')dnl
define(conref_op_edge_3d,`dnl
con_refine_op_subroutine_head_3d()dnl
      $1
     &  arrayc(EDGE3d$2(cilo,cihi,0)),
     &  arrayf(EDGE3d$2(filo,fihi,0))
      integer ic0,ic1,ic2,if0,if1,if2
conref_edge_body_3d($2,$3,$4)dnl
')dnl
define(conref_op_face_3d,`dnl
con_refine_op_subroutine_head_3d()dnl
      $1
     &  arrayc(FACE3d$2(cilo,cihi,0)),
     &  arrayf(FACE3d$2(filo,fihi,0))
      integer ie$2,ic$3,ic$4,if$2,if$3,if$4
conref_face_body_3d($2,$3,$4)dnl
')dnl
define(conref_op_outerface_3d,`dnl
con_refine_op_subroutine_head_3d()dnl
      $1
     &  arrayc(OUTERFACE3d$2(cilo,cihi,0)),
     &  arrayf(OUTERFACE3d$2(filo,fihi,0))
      integer ic$3,ic$4,if$3,if$4
conref_outerface_body_3d($2,$3,$4)dnl
')dnl
define(conref_op_side_3d,`dnl
con_refine_op_subroutine_head_3d()dnl
      $1
     &  arrayc(SIDE3d$2(cilo,cihi,0)),
     &  arrayf(SIDE3d$2(filo,fihi,0))
      integer ic0,ic1,ic2,if0,if1,if2
conref_side_body_3d($2,$3,$4)dnl
')dnl

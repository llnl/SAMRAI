c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for 3d standard linear time interpolation 
c                operators.
c
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl
define(lin_time_int_subroutine_head_3d,`dnl
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ifirst1,ifirst2,ilast0,ilast1,ilast2,
     &  oilo0,oilo1,oilo2,oihi0,oihi1,oihi2,
     &  nilo0,nilo1,nilo2,nihi0,nihi1,nihi2,
     &  dilo0,dilo1,dilo2,dihi0,dihi1,dihi2
      double precision
     &  tfrac, oldfrac
')dnl
c
define(lin_time_int_body_3d,`dnl
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do $3=$6
         do $2=$5
            do $1=$4
               arraydst($1,$2,$3)=
     &                              +arrayold($1,$2,$3)*oldfrac
     &                              +arraynew($1,$2,$3)*tfrac
            enddo
         enddo
      enddo
c
      return
      end
')dnl
c
define(lin_time_int_op_cell_3d,`dnl
lin_time_int_subroutine_head_3d()dnl
      $1
     &  arrayold(CELL3d(oilo,oihi,0)),
     &  arraynew(CELL3d(nilo,nihi,0)),
     &  arraydst(CELL3d(dilo,dihi,0))
      integer ic0,ic1,ic2
lin_time_int_body_3d(`ic0',`ic1',`ic2',`ifirst0,ilast0',`ifirst1,ilast1',`ifirst2,ilast2')dnl
')dnl
c
define(lin_time_int_op_edge_3d,`dnl
lin_time_int_subroutine_head_3d()dnl
      $1
     &  arrayold(EDGE3d$2(oilo,oihi,0)),
     &  arraynew(EDGE3d$2(nilo,nihi,0)),
     &  arraydst(EDGE3d$2(dilo,dihi,0))
      integer ie0,ie1,ie2
ifelse($2,`0',`lin_time_int_body_3d(`ie0',`ie1',`ie2',`ifirst0,ilast0',`ifirst1,ilast1',`ifirst2,ilast2')dnl
',`')dnl
ifelse($2,`1',`lin_time_int_body_3d(`ie0',`ie1',`ie2',`ifirst0,ilast0',`ifirst1,ilast1',`ifirst2,ilast2')dnl
',`')dnl
ifelse($2,`2',`lin_time_int_body_3d(`ie0',`ie1',`ie2',`ifirst0,ilast0',`ifirst1,ilast1',`ifirst2,ilast2')dnl
',`')dnl
')dnl
c
define(lin_time_int_op_face_3d,`dnl
lin_time_int_subroutine_head_3d()dnl
      $1
     &  arrayold(FACE3d$2(oilo,oihi,0)),
     &  arraynew(FACE3d$2(nilo,nihi,0)),
     &  arraydst(FACE3d$2(dilo,dihi,0))
      integer ie0,ic1,ic2
lin_time_int_body_3d(`ie0',`ic1',`ic2',`ifirst0,ilast0',`ifirst1,ilast1',`ifirst2,ilast2')dnl
')dnl
c
define(lin_time_int_op_node_3d,`dnl
lin_time_int_subroutine_head_3d()dnl
      $1
     &  arrayold(NODE3d(oilo,oihi,0)),
     &  arraynew(NODE3d(nilo,nihi,0)),
     &  arraydst(NODE3d(dilo,dihi,0))
      integer ie0,ie1,ie2
lin_time_int_body_3d(`ie0',`ie1',`ie2',`ifirst0,ilast0',`ifirst1,ilast1',`ifirst2,ilast2')dnl
')dnl
c
define(lin_time_int_op_outerface_3d,`dnl
lin_time_int_subroutine_head_3d()dnl
      $1
     &  arrayold(OUTERFACE3d$2(oilo,oihi,0)),
     &  arraynew(OUTERFACE3d$2(nilo,nihi,0)),
     &  arraydst(OUTERFACE3d$2(dilo,dihi,0))
      integer ic$3,ic$4 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic$4=ifirst$4,ilast$4
         do ic$3=ifirst$3,ilast$3
            arraydst(ic$3,ic$4)=arrayold(ic$3,ic$4)*oldfrac
     &                       +arraynew(ic$3,ic$4)*tfrac
         enddo
      enddo
c
      return
      end
')dnl
c
define(lin_time_int_op_outerside_3d,`dnl
lin_time_int_subroutine_head_3d()dnl
      $1
     &  arrayold(OUTERSIDE3d$2(oilo,oihi,0)),
     &  arraynew(OUTERSIDE3d$2(nilo,nihi,0)),
     &  arraydst(OUTERSIDE3d$2(dilo,dihi,0))
      integer ic$3,ic$4 
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do ic$4=ifirst$4,ilast$4
         do ic$3=ifirst$3,ilast$3
            arraydst(ic$3,ic$4)=arrayold(ic$3,ic$4)*oldfrac
     &                       +arraynew(ic$3,ic$4)*tfrac
         enddo
      enddo
c
      return
      end
')dnl
c
define(lin_time_int_op_side_3d,`dnl
lin_time_int_subroutine_head_3d()dnl
      $1
     &  arrayold(SIDE3d$2(oilo,oihi,0)),
     &  arraynew(SIDE3d$2(nilo,nihi,0)),
     &  arraydst(SIDE3d$2(dilo,dihi,0))
      integer ie0,ie1,ie2
ifelse($2,`0',`lin_time_int_body_3d(`ie0',`ie1',`ie2',`ifirst0,ilast0',`ifirst1,ilast1',`ifirst2,ilast2')dnl
',`')dnl
ifelse($2,`1',`lin_time_int_body_3d(`ie0',`ie1',`ie2',`ifirst0,ilast0',`ifirst1,ilast1',`ifirst2,ilast2')dnl
',`')dnl
ifelse($2,`2',`lin_time_int_body_3d(`ie0',`ie1',`ie2',`ifirst0,ilast0',`ifirst1,ilast1',`ifirst2,ilast2')dnl
',`')dnl
')dnl
c

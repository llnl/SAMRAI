c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for 1d standard linear time interpolation
c                operators.
c
include(PDAT_FORTDIR/pdat_m4arrdim1d.i)dnl
define(lin_time_int_subroutine_head_1d,`dnl
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0,
     &  tfrac,
     &  arrayold,arraynew,
     &  arraydst)
c***********************************************************************
      implicit none
      double precision one
      parameter (one=1.d0)
c
      integer
     &  ifirst0,ilast0,
     &  oilo0,oihi0,
     &  nilo0,nihi0,
     &  dilo0,dihi0
      double precision
     &  tfrac, oldfrac
')dnl
c
define(lin_time_int_body_1d,`dnl
c
c***********************************************************************
c
      oldfrac=one-tfrac

      do $1=$2
         arraydst($1)=arrayold($1)*oldfrac
     &                +arraynew($1)*tfrac
      enddo
c
      return
      end
')dnl
c
define(lin_time_int_op_cell_1d,`dnl
lin_time_int_subroutine_head_1d()dnl
      $1
     &  arrayold(CELL1d(oilo,oihi,0)),
     &  arraynew(CELL1d(nilo,nihi,0)),
     &  arraydst(CELL1d(dilo,dihi,0))
      integer ic0
lin_time_int_body_1d(`ic0',`ifirst0,ilast0')dnl
')dnl
c
define(lin_time_int_op_edge_1d,`dnl
lin_time_int_subroutine_head_1d()dnl
      $1
     &  arrayold(EDGE1d(oilo,oihi,0)),
     &  arraynew(EDGE1d(nilo,nihi,0)),
     &  arraydst(EDGE1d(dilo,dihi,0))
      integer ie0
lin_time_int_body_1d(`ie0',`ifirst0,ilast0')dnl
')dnl
c
define(lin_time_int_op_face_1d,`dnl
lin_time_int_subroutine_head_1d()dnl
      $1
     &  arrayold(FACE1d(oilo,oihi,0)),
     &  arraynew(FACE1d(nilo,nihi,0)),
     &  arraydst(FACE1d(dilo,dihi,0))
      integer ie0
lin_time_int_body_1d(`ie0',`ifirst0,ilast0')dnl
')dnl
c
define(lin_time_int_op_node_1d,`dnl
lin_time_int_subroutine_head_1d()dnl
      $1
     &  arrayold(NODE1d(oilo,oihi,0)),
     &  arraynew(NODE1d(nilo,nihi,0)),
     &  arraydst(NODE1d(dilo,dihi,0))
      integer ie0
lin_time_int_body_1d(`ie0',`ifirst0,ilast0')dnl
')dnl
c
define(lin_time_int_op_outerface_1d,`dnl
lin_time_int_subroutine_head_1d()dnl
      $1
     &  arrayold(OUTERFACE1d(oilo,oihi,0)),
     &  arraynew(OUTERFACE1d(nilo,nihi,0)),
     &  arraydst(OUTERFACE1d(dilo,dihi,0))
c
c***********************************************************************
c
      oldfrac=one-tfrac

      arraydst(1)=arrayold(1)*oldfrac
     &           +arraynew(1)*tfrac
c
      return
      end
')dnl
c
define(lin_time_int_op_outerside_1d,`dnl
lin_time_int_subroutine_head_1d()dnl
      $1
     &  arrayold(OUTERSIDE1d(oilo,oihi,0)),
     &  arraynew(OUTERSIDE1d(nilo,nihi,0)),
     &  arraydst(OUTERSIDE1d(dilo,dihi,0))
c
c***********************************************************************
c
      oldfrac=one-tfrac

      arraydst(1)=arrayold(1)*oldfrac
     &           +arraynew(1)*tfrac
c
      return
      end
')dnl
c
define(lin_time_int_op_side_1d,`dnl
lin_time_int_subroutine_head_1d()dnl
      $1
     &  arrayold(SIDE1d(oilo,oihi,0)),
     &  arraynew(SIDE1d(nilo,nihi,0)),
     &  arraydst(SIDE1d(dilo,dihi,0))
      integer ie0
lin_time_int_body_1d(`ie0',`ifirst0,ilast0')dnl
')dnl

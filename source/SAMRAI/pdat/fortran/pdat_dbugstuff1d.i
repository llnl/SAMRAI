c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for 1d patchdata debugging routines.
c
include(PDAT_FORTDIR/pdat_m4arrdim1d.i)dnl
c
define(pdat_debug_subroutine_head_1d,`dnl
     &  fi0,la0,ng,ibeg0,iend0,array)
c     =============================================================
      implicit none
      integer fi0,la0,ng,ibeg0,iend0
')dnl
c
define(pdat_debug_body_1d,`dnl

      do $1=$2
         write(6,*) "array[",$1,"] = ",array($1)
      enddo

      call flush(6)
      return
      end
')dnl
define(pdat_debug_cell_1d,`dnl
pdat_debug_subroutine_head_1d()dnl
      $1
     &  array(CELL1d(fi,la,ng))
      integer ic0
c     =============================================================
pdat_debug_body_1d(`ic0',`ibeg0,iend0')dnl
')dnl
define(pdat_debug_face_1d,`dnl
pdat_debug_subroutine_head_1d()dnl
      $1
     &  array(FACE1d(fi,la,ng))
      integer ie0
c     =============================================================
pdat_debug_body_1d(`ie0',`ibeg0,iend0+1')dnl
')dnl
define(pdat_debug_node_1d,`dnl
pdat_debug_subroutine_head_1d()dnl
      $1
     &  array(NODE1d(fi,la,ng))
      integer ie0
c     =============================================================
pdat_debug_body_1d(`ie0',`ibeg0,iend0+1')dnl
')dnl
define(pdat_debug_outerface_1d,`dnl
pdat_debug_subroutine_head_1d()dnl
      $1
     &  array(OUTERFACE1d(fi,la,ng))
c     =============================================================

      write(6,*) "array[",1,"] = ",array(1)

      call flush(6)
      return
      end
')dnl
define(pdat_debug_outerside_1d,`dnl
pdat_debug_subroutine_head_1d()dnl
      $1
     &  array(OUTERSIDE1d(fi,la,ng))
c     =============================================================

      write(6,*) "array[",1,"] = ",array(1)

      call flush(6)
      return
      end
')dnl
define(pdat_debug_side_1d,`dnl
pdat_debug_subroutine_head_1d()dnl
      $1
     &  array(SIDE1d(fi,la,ng))
      integer ie0
c     =============================================================
pdat_debug_body_1d(`ie0',`ibeg0,iend0+1')dnl
')dnl
define(pdat_debug_edge_1d,`dnl
pdat_debug_subroutine_head_1d()dnl
      $1
     &  array(EDGE1d(fi,la,ng))
      integer ie0
c     =============================================================
pdat_debug_body_1d(`ie0',`ibeg0,iend0')dnl
')dnl

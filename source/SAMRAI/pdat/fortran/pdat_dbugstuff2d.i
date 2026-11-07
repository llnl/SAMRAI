c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for 2d patchdata debugging routines.
c
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl
c
define(pdat_debug_subroutine_head_2d,`dnl
     &  fi0,la0,fi1,la1,ng,
     &  ibeg0,iend0,ibeg1,iend1,array)
c     =============================================================
      implicit none
      integer fi0,la0,fi1,la1,ng,
     &   ibeg0,iend0,ibeg1,iend1
')dnl
c
define(pdat_debug_body_2d,`dnl

      do $2=$4
         do $1=$3 
            write(6,*) "array[",$1,",",$2,"] = ",array($1,$2)
         enddo
      enddo

      call flush(6)
      return
      end
')dnl
define(pdat_debug_cell_2d,`dnl
pdat_debug_subroutine_head_2d()dnl
      $1
     &  array(CELL2d(fi,la,ng))
      integer ic0,ic1
c     =============================================================
pdat_debug_body_2d(`ic0',`ic1',`ibeg0,iend0',`ibeg1,iend1')dnl
')dnl
define(pdat_debug_face_2d,`dnl
pdat_debug_subroutine_head_2d()dnl
      $1
     &  array(FACE2d$2(fi,la,ng))
      integer ie$2,ic$3
c     =============================================================
pdat_debug_body_2d(`ie$2',`ic$3',`ibeg$2,iend$2+1',`ibeg$3,iend$3')dnl
')dnl
define(pdat_debug_node_2d,`dnl
pdat_debug_subroutine_head_2d()dnl
      $1
     &  array(NODE2d(fi,la,ng))
      integer ie0,ie1
c     =============================================================
pdat_debug_body_2d(`ie0',`ie1',`ibeg0,iend0+1',`ibeg1,iend1+1')dnl
')dnl
define(pdat_debug_outerface_2d,`dnl
pdat_debug_subroutine_head_2d()dnl
      $1
     &  array(OUTERFACE2d$2(fi,la,ng))
      integer ic$3
c     =============================================================

      do ic$3=ibeg$3,iend$3
         write(6,*) "array[",ic$3,"] = ",array(ic$3)
      enddo
 
      call flush(6)
      return
      end
')dnl
define(pdat_debug_outerside_2d,`dnl
pdat_debug_subroutine_head_2d()dnl
      $1
     &  array(OUTERSIDE2d$2(fi,la,ng))
      integer ic$3
c     =============================================================

      do ic$3=ibeg$3,iend$3
         write(6,*) "array[",ic$3,"] = ",array(ic$3)
      enddo
 
      call flush(6)
      return
      end
')dnl
define(pdat_debug_side_2d,`dnl
pdat_debug_subroutine_head_2d()dnl
      $1
     &  array(SIDE2d$2(fi,la,ng))
      integer ie0,ic1
c     =============================================================
pdat_debug_body_2d(`ie0',`ic1',`ibeg0,iend0+$3',`ibeg1,iend1+$2')dnl
')dnl
define(pdat_debug_edge_2d,`dnl
pdat_debug_subroutine_head_2d()dnl
      $1
     &  array(EDGE2d$2(fi,la,ng))
      integer ie0,ic1
c     =============================================================
pdat_debug_body_2d(`ie0',`ic1',`ibeg0,iend0+$2',`ibeg1,iend1+$3')dnl
')dnl

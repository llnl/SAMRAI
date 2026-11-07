c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for Cartesian 3d boundary conditions.
c
define(NDIM,3)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl
c
c general boundary condition cases
c
define(setvalues_to_ict,`dnl
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
')dnl
define(setvalues_to_ict_reflect,`dnl
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
         arrdata(ic0,ic1,ic2,$1)= -arrdata(ic0,ic1,ic2,$1)
')dnl
define(setvalues_to_dirichlet,`dnl
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=face_values(k,$1)
         enddo
')dnl
define(setvalues_to_neumann,`dnl
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k) + $2
     &       * face_values(k,$3) * dble(ic$1-ict$1) * dx($1)
         enddo
')dnl
c
c edge boundary conditions
c
define(do_edge_ict,`dnl
             do ic2=ibeg2,iend2
               ict2 = $3
               do ic1=ibeg1,iend1
                 ict1 = $2
                 do ic0=ibeg0,iend0
                   ict0 = $1
setvalues_to_ict()dnl
                   enddo
                enddo
             enddo
')dnl
define(do_edge_ict_reflect,`dnl
             do ic2=ibeg2,iend2
               ict2 = $3
               do ic1=ibeg1,iend1
                 ict1 = $2
                 do ic0=ibeg0,iend0
                   ict0 = $1
setvalues_to_ict_reflect($4)dnl
                   enddo
                enddo
             enddo
')dnl
define(do_edge_dirichlet,`dnl
             do ic2=ibeg2,iend2
               do ic1=ibeg1,iend1
                 do ic0=ibeg0,iend0
setvalues_to_dirichlet($1)dnl
                   enddo
                enddo
             enddo
')dnl
define(do_edge_neumann,`dnl
             do ic2=ibeg2,iend2
               ict2 = $3
               do ic1=ibeg1,iend1
                 ict1 = $2
                 do ic0=ibeg0,iend0
                   ict0 = $1
setvalues_to_neumann($4,$5,$6)dnl
                   enddo
                enddo
             enddo
')dnl
define(do_bdry_edge,`dnl
           if (bdry_cond.eq.XDIRICHLET) then
do_edge_dirichlet(`face_loc0')dnl
           else if (bdry_cond.eq.YDIRICHLET) then
do_edge_dirichlet(`face_loc1')dnl
           else if (bdry_cond.eq.ZDIRICHLET) then
do_edge_dirichlet(`face_loc2')dnl
           else if (bdry_cond.eq.XNEUMANN) then
do_edge_neumann(`ipivot0',`ic1',`ic2',0,`dirsign0',`face_loc0')dnl
           else if (bdry_cond.eq.YNEUMANN) then
do_edge_neumann(`ic0',`ipivot1',`ic2',1,`dirsign1',`face_loc1')dnl
           else if (bdry_cond.eq.ZNEUMANN) then
do_edge_neumann(`ic0',`ic1',`ipivot2',2,`dirsign2',`face_loc2')dnl
           else if (bdry_cond.eq.XFLOW) then
do_edge_ict(`ipivot0',`ic1',`ic2')dnl
           else if (bdry_cond.eq.YFLOW) then
do_edge_ict(`ic0',`ipivot1',`ic2')dnl
           else if (bdry_cond.eq.ZFLOW) then
do_edge_ict(`ic0',`ic1',`ipivot2')dnl
           else if (bdry_cond.eq.XREFLECT) then
do_edge_ict_reflect(`ipivot0',`ic1',`ic2',0)dnl
           else if (bdry_cond.eq.YREFLECT) then
do_edge_ict_reflect(`ic0',`ipivot1',`ic2',1)dnl
           else if (bdry_cond.eq.ZREFLECT) then
do_edge_ict_reflect(`ic0',`ic1',`ipivot2',2)dnl
           else
             write(6,*) "INVALID EDGE bdry_cond in getcartedgebdry3d"
             write(6,*) "bdry_loc = ",bdry_loc,
     &                  "bdry_cond = ",bdry_cond
           endif
')dnl
c
c node boundary conditions
c
define(do_node_ict,`dnl
             do ic2=ibeg2,iend2
               ict2 = $3
               do ic1=ibeg1,iend1
                 ict1 = $2
                 do ic0=ibeg0,iend0
                   ict0 = $1
setvalues_to_ict()dnl
                   enddo
                enddo
             enddo
')dnl
define(do_node_ict_reflect,`dnl
             do ic2=ibeg2,iend2
               ict2 = $3
               do ic1=ibeg1,iend1
                 ict1 = $2
                 do ic0=ibeg0,iend0
                   ict0 = $1
setvalues_to_ict_reflect($4)dnl
                   enddo
                enddo
             enddo
')dnl
define(do_node_dirichlet,`dnl
             do ic2=ibeg2,iend2
               do ic1=ibeg1,iend1
                 do ic0=ibeg0,iend0
setvalues_to_dirichlet($1)dnl
                   enddo
                enddo
             enddo
')dnl
define(do_node_neumann,`dnl
             do ic2=ibeg2,iend2
               ict2 = $3
               do ic1=ibeg1,iend1
                 ict1 = $2
                 do ic0=ibeg0,iend0
                   ict0 = $1
setvalues_to_neumann($4,$5,$6)dnl
                   enddo
                enddo
             enddo
')dnl
define(do_bdry_node,`dnl
           if (bdry_cond.eq.XDIRICHLET) then
do_node_dirichlet(`face_loc0')dnl
           else if (bdry_cond.eq.YDIRICHLET) then
do_node_dirichlet(`face_loc1')dnl
           else if (bdry_cond.eq.ZDIRICHLET) then
do_node_dirichlet(`face_loc2')dnl
          else if (bdry_cond.eq.XNEUMANN) then
do_node_neumann(`ipivot0',`ic1',`ic2',0,`dirsign0',`face_loc0')dnl
          else if (bdry_cond.eq.YNEUMANN) then
do_node_neumann(`ic0',`ipivot1',`ic2',1,`dirsign1',`face_loc1')dnl
          else if (bdry_cond.eq.ZNEUMANN) then
do_node_neumann(`ic0',`ic1',`ipivot2',2,`dirsign2',`face_loc2')dnl
           else if (bdry_cond.eq.XFLOW) then
do_node_ict(`ipivot0',`ic1',`ic2')dnl
           else if (bdry_cond.eq.YFLOW) then
do_node_ict(`ic0',`ipivot1',`ic2')dnl
           else if (bdry_cond.eq.ZFLOW) then
do_node_ict(`ic0',`ic1',`ipivot2')dnl
           else if (bdry_cond.eq.XREFLECT) then
do_node_ict_reflect(`ipivot0',`ic1',`ic2',0)dnl
           else if (bdry_cond.eq.YREFLECT) then
do_node_ict_reflect(`ic0',`ipivot1',`ic2',1)dnl
           else if (bdry_cond.eq.ZREFLECT) then
do_node_ict_reflect(`ic0',`ic1',`ipivot2',2)dnl
           else
             write(6,*) "INVALID NODE bdry_cond in getcartnodebdry3d"
             write(6,*) "bdry_loc = ",bdry_loc,
     &                  "bdry_cond = ",bdry_cond
           endif
')dnl
c
c face boundary conditions
c
define(do_bdry_face,`dnl
           if (bdry_cond.eq.DIRICHLET) then
             do ic2=ibeg2,iend2
               do ic1=ibeg1,iend1
                 do ic0=ibeg0,iend0
setvalues_to_dirichlet(`face_loc')dnl
                 enddo
               enddo
             enddo
           else if (bdry_cond.eq.NEUMANN) then
             do ic2=ibeg2,iend2
               ict2 = $4
               do ic1=ibeg1,iend1
                 ict1 = $3
                 do ic0=ibeg0,iend0
                   ict0 = $2
setvalues_to_neumann($1,$5,`face_loc')dnl
                 enddo
               enddo
             enddo
           else if (bdry_cond.eq.FLOW) then
             do ic2=ibeg2,iend2
               ict2 = $4
               do ic1=ibeg1,iend1
                 ict1 = $3
                 do ic0=ibeg0,iend0
                   ict0 = $2
setvalues_to_ict()dnl
                 enddo
               enddo
             enddo
           else if (bdry_cond.eq.REFLECT) then
             do ic2=ibeg2,iend2
              ict2 = $4
               do ic1=ibeg1,iend1
                 ict1 = $3
                 do ic0=ibeg0,iend0
                   ict0 = $2
setvalues_to_ict_reflect($1)dnl
                 enddo
               enddo
             enddo
          else
             write(6,*) "INVALID FACE bdry_cond in getcartfacebdry3d"
             write(6,*) "bdry_loc = ",bdry_loc,
     &                  "bdry_cond = ",bdry_cond
          endif
')dnl
c***********************************************************************
c***********************************************************************
      subroutine getcartfacebdry3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,
     &  ngc0,ngc1,ngc2,
     &  dx,
     &  bdry_loc,
     &  bdry_cond,
     &  face_values,
     &  arrdata,
     &  arrdepth)
c***********************************************************************
      implicit none
include(FORTDIR/appu_cartbdryparams3d.i)dnl
c***********************************************************************
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      integer ngc0,ngc1,ngc2
c
      REAL dx(0:NDIM-1)
c
      integer bdry_loc,
     &        arrdepth, 
     &        bdry_cond
      REAL
     &        face_values(0:arrdepth-1,0:2*NDIM-1)
c
      REAL
     &  arrdata(CELL3dVECG(ifirst,ilast,ngc),0:arrdepth-1)
c
      integer ic2,ic1,ic0,ict0,ict1,ict2
      integer k
      integer ipivot0,ipivot1,ipivot2
      integer face_loc
      REAL    dirsign0,dirsign1,dirsign2
c
c***********************************************************************
c   bdry_loc   index     bdry_loc    index     bdry_loc   index
c       0     (-1, 0, 0)     2     (0,-1, 0)      4     (0, 0,-1)
c       1     ( 1, 0, 0)     3     (0, 1, 0)      5     (0, 0, 1)
c***********************************************************************
c***********************************************************************
      face_loc = bdry_loc
      if ((bdry_loc.eq.XLEFT).or.(bdry_loc.eq.XRIGHT)) then
         if (bdry_loc.eq.XLEFT) then
c                             x0 boundary
            ipivot0  = ifirst0
            dirsign0 = -1.d0 
         else
c                             x1 boundary
            ipivot0  = ilast0
            dirsign0 = 1.d0 
         endif
do_bdry_face(0,`ipivot0',`ic1',`ic2',`dirsign0')dnl
      else if ((bdry_loc.eq.YLEFT).or.(bdry_loc.eq.YRIGHT)) then
         if (bdry_loc.eq.YLEFT) then
c                             y0 boundary
            ipivot1  = ifirst1
            dirsign1 = -1.d0 
         else
c                             y1 boundary
            ipivot1  = ilast1
            dirsign1 = 1.d0 
         endif
do_bdry_face(1,`ic0',`ipivot1',`ic2',`dirsign1')dnl
      else if ((bdry_loc.eq.ZLEFT).or.(bdry_loc.eq.ZRIGHT)) then
         if (bdry_loc.eq.ZLEFT) then
c                             z0 boundary
            ipivot2  = ifirst2
            dirsign2 = -1.d0 
         else
c                             z1 boundary
            ipivot2  = ilast2
            dirsign2 = 1.d0 
         endif
do_bdry_face(2,`ic0',`ic1',`ipivot2',`dirsign2')dnl
      endif
c
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine getcartedgebdry3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,
     &  ngc0,ngc1,ngc2,
     &  dx,
     &  bdry_loc,
     &  bdry_cond,
     &  face_values,
     &  arrdata,
     &  arrdepth)
c***********************************************************************
      implicit none
include(FORTDIR/appu_cartbdryparams3d.i)dnl
c***********************************************************************
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      integer ngc0,ngc1,ngc2
c
      REAL dx(0:NDIM-1)
c
      integer bdry_loc,
     &        arrdepth,
     &        bdry_cond
      REAL
     &        face_values(0:arrdepth-1,0:2*NDIM-1)
c
      REAL
     &  arrdata(CELL3dVECG(ifirst,ilast,ngc),0:arrdepth-1)
c
      integer ic2,ic1,ic0,ict0,ict1,ict2
      integer k
      integer ipivot0,ipivot1,ipivot2
      integer face_loc0,face_loc1,face_loc2
      REAL    dirsign0,dirsign1,dirsign2
c
c***********************************************************************
c***********************************************************************
c    bdry_loc   index    bdry_loc   index     bdry_loc   index
c       0     (0,-1,-1)      4     (-1,0,-1)      8     (-1,-1,0)
c       1     (0, 1,-1)      5     (-1,0, 1)      9     ( 1,-1,0)
c       2     (0,-1, 1)      6     ( 1,0,-1)     10     (-1, 1,0)
c       3     (0, 1, 1)      7     ( 1,0, 1)     11     ( 1, 1,0)
c***********************************************************************
c***********************************************************************
      if ((bdry_loc.eq.Y0Z0).or.(bdry_loc.eq.Y1Z0).or.
     &    (bdry_loc.eq.Y0Z1).or.(bdry_loc.eq.Y1Z1)) then
         if (bdry_loc.eq.Y0Z0) then
            face_loc1 = YLEFT
            ipivot1   = ifirst1
            dirsign1  = -1.d0
            face_loc2 = ZLEFT
            ipivot2   = ifirst2
            dirsign2  = -1.d0
         else if (bdry_loc.eq.Y1Z0) then
            face_loc1 = YRIGHT
            ipivot1   = ilast1
            dirsign1  = 1.d0
            face_loc2 = ZLEFT
            ipivot2   = ifirst2
            dirsign2  = -1.d0
         else if (bdry_loc.eq.Y0Z1) then
            face_loc1 = YLEFT
            ipivot1   = ifirst1
            dirsign1  = -1.d0
            face_loc2 = ZRIGHT
            ipivot2   = ilast2
            dirsign2  = 1.d0
         else
            face_loc1 = YRIGHT
            ipivot1   = ilast1
            dirsign1  = 1.d0
            face_loc2 = ZRIGHT
            ipivot2   = ilast2
            dirsign2  = 1.d0
         endif
do_bdry_edge()dnl
      else if ((bdry_loc.eq.X0Z0).or.(bdry_loc.eq.X0Z1).or.
     &         (bdry_loc.eq.X1Z0).or.(bdry_loc.eq.X1Z1)) then
         if (bdry_loc.eq.X0Z0) then
            face_loc2 = ZLEFT
            ipivot2   = ifirst2
            dirsign2  = -1.d0
            face_loc0 = XLEFT
            ipivot0   = ifirst0
            dirsign0  = -1.d0
         else if (bdry_loc.eq.X0Z1) then
            face_loc2 = ZRIGHT
            ipivot2   = ilast2
            dirsign2  = 1.d0
            face_loc0 = XLEFT
            ipivot0   = ifirst0
            dirsign0  = -1.d0
         else if (bdry_loc.eq.X1Z0) then
            face_loc2 = ZLEFT
            ipivot2   = ifirst2
            dirsign2  = -1.d0
            face_loc0 = XRIGHT
            ipivot0   = ilast0
            dirsign0  = 1.d0
         else
            face_loc2 = ZRIGHT
            ipivot2   = ilast2
            dirsign2  = 1.d0
            face_loc0 = XRIGHT
            ipivot0   = ilast0
            dirsign0  = 1.d0
         endif
do_bdry_edge()dnl
      else if ((bdry_loc.eq.X0Y0).or.(bdry_loc.eq.X1Y0).or.
     &         (bdry_loc.eq.X0Y1).or.(bdry_loc.eq.X1Y1)) then
         if (bdry_loc.eq.X0Y0) then
            face_loc1 = YLEFT
            ipivot1   = ifirst1
            dirsign1  = -1.d0
            face_loc0 = XLEFT
            ipivot0   = ifirst0
            dirsign0  = -1.d0
         else if (bdry_loc.eq.X1Y0) then
            face_loc1 = YLEFT
            ipivot1   = ifirst1
            dirsign1  = -1.d0
            face_loc0 = XRIGHT
            ipivot0   = ilast0
            dirsign0  = 1.d0
         else if (bdry_loc.eq.X0Y1) then
            face_loc1 = YRIGHT
            ipivot1   = ilast1
            dirsign1  = 1.d0
            face_loc0 = XLEFT
            ipivot0   = ifirst0
            dirsign0  = -1.d0
         else
            face_loc1 = YRIGHT
            ipivot1   = ilast1
            dirsign1  = 1.d0
            face_loc0 = XRIGHT
            ipivot0   = ilast0
            dirsign0  = 1.d0
         endif
do_bdry_edge()dnl
      endif
c
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine getcartnodebdry3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,
     &  ngc0,ngc1,ngc2,
     &  dx,
     &  bdry_loc,
     &  bdry_cond,
     &  face_values,
     &  arrdata,
     &  arrdepth)
c***********************************************************************
      implicit none
include(FORTDIR/appu_cartbdryparams3d.i)dnl
c***********************************************************************
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      integer ngc0,ngc1,ngc2
c
      REAL dx(0:NDIM-1)
c
      integer bdry_loc,
     &        arrdepth,
     &        bdry_cond
      REAL
     &        face_values(0:arrdepth-1,0:2*NDIM-1)
c
      REAL
     &  arrdata(CELL3dVECG(ifirst,ilast,ngc),0:arrdepth-1)
c
      integer ic2,ic1,ic0,ict0,ict1,ict2
      integer k
      integer ipivot0,ipivot1,ipivot2
      integer face_loc0,face_loc1,face_loc2
      REAL    dirsign0,dirsign1,dirsign2
c
c***********************************************************************
c***********************************************************************
c    bdry_loc   index     bdry_loc    index
c       0     (-1,-1,-1)      4     (-1,-1, 1)
c       1     ( 1,-1,-1)      5     ( 1,-1, 1)
c       2     (-1, 1,-1)      6     (-1, 1, 1)
c       3     ( 1, 1,-1)      7     ( 1, 1, 1)
c***********************************************************************
c***********************************************************************
      if (bdry_loc.eq.X0Y0Z0) then
         face_loc0 = XLEFT
         ipivot0   = ifirst0
         dirsign0  = -1.d0
         face_loc1 = YLEFT
         ipivot1   = ifirst1
         dirsign1  = -1.d0
         face_loc2 = ZLEFT
         ipivot2   = ifirst2
         dirsign2  = -1.d0
      else if (bdry_loc.eq.X1Y0Z0) then
         face_loc0 = XRIGHT
         ipivot0   = ilast0
         dirsign0  = 1.d0
         face_loc1 = YLEFT
         ipivot1   = ifirst1
         dirsign1  = -1.d0
         face_loc2 = ZLEFT
         ipivot2   = ifirst2
         dirsign2  = -1.d0
      else if (bdry_loc.eq.X0Y1Z0) then
         face_loc0 = XLEFT
         ipivot0   = ifirst0
         dirsign0  = -1.d0
         face_loc1 = YRIGHT
         ipivot1   = ilast1
         dirsign1  = 1.d0
         face_loc2 = ZLEFT
         ipivot2   = ifirst2
         dirsign2  = -1.d0
      else if (bdry_loc.eq.X1Y1Z0) then
         face_loc0 = XRIGHT
         ipivot0   = ilast0
         dirsign0  = 1.d0
         face_loc1 = YRIGHT
         ipivot1   = ilast1
         dirsign1  = 1.d0
         face_loc2 = ZLEFT
         ipivot2   = ifirst2
         dirsign2  = -1.d0
      else if (bdry_loc.eq.X0Y0Z1) then
         face_loc0 = XLEFT
         ipivot0   = ifirst0
         dirsign0  = -1.d0
         face_loc1 = YLEFT
         ipivot1   = ifirst1
         dirsign1  = -1.d0
         face_loc2 = ZRIGHT
         ipivot2   = ilast2
         dirsign2  = 1.d0
      else if (bdry_loc.eq.X1Y0Z1) then
         face_loc0 = XRIGHT
         ipivot0   = ilast0
         dirsign0  = 1.d0
         face_loc1 = YLEFT
         ipivot1   = ifirst1
         dirsign1  = -1.d0
         face_loc2 = ZRIGHT
         ipivot2   = ilast2
         dirsign2  = 1.d0
      else if (bdry_loc.eq.X0Y1Z1) then
         face_loc0 = XLEFT
         ipivot0   = ifirst0
         dirsign0  = -1.d0
         face_loc1 = YRIGHT
         ipivot1   = ilast1
         dirsign1  = 1.d0
         face_loc2 = ZRIGHT
         ipivot2   = ilast2
         dirsign2  = 1.d0
      else if (bdry_loc.eq.X1Y1Z1) then
         face_loc0 = XRIGHT
         ipivot0   = ilast0
         dirsign0  = 1.d0
         face_loc1 = YRIGHT
         ipivot1   = ilast1
         dirsign1  = 1.d0
         face_loc2 = ZRIGHT
         ipivot2   = ilast2
         dirsign2  = 1.d0
      endif
do_bdry_node()dnl
c
      return
      end

c***********************************************************************
c***********************************************************************
      subroutine stufcartbdryloc3d(
     &  XLEFTin, XRIGHTin, YLEFTin, YRIGHTin, ZLEFTin, ZRIGHTin,
     &  Y0Z0in, Y1Z0in, Y0Z1in, Y1Z1in,
     &  X0Z0in, X0Z1in, X1Z0in, X1Z1in,
     &  X0Y0in, X1Y0in, X0Y1in, X1Y1in,
     &  X0Y0Z0in, X1Y0Z0in, X0Y1Z0in, X1Y1Z0in,
     &  X0Y0Z1in, X1Y0Z1in, X0Y1Z1in, X1Y1Z1in)
      implicit none
      integer
     &  XLEFTin, XRIGHTin, YLEFTin, YRIGHTin, ZLEFTin, ZRIGHTin,
     &  Y0Z0in, Y1Z0in, Y0Z1in, Y1Z1in,
     &  X0Z0in, X0Z1in, X1Z0in, X1Z1in,
     &  X0Y0in, X1Y0in, X0Y1in, X1Y1in,
     &  X0Y0Z0in, X1Y0Z0in, X0Y1Z0in, X1Y1Z0in,
     &  X0Y0Z1in, X1Y0Z1in, X0Y1Z1in, X1Y1Z1in
include(FORTDIR/appu_cartbdryparams3d.i)dnl

c 3d faces
        XLEFT=XLEFTin
        XRIGHT=XRIGHTin
        YLEFT=YLEFTin
        YRIGHT=YRIGHTin
        ZLEFT=ZLEFTin
        ZRIGHT=ZRIGHTin

c 3d edges
        Y0Z0=Y0Z0in
        Y1Z0=Y1Z0in
        Y0Z1=Y0Z1in
        Y1Z1=Y1Z1in
        X0Z0=X0Z0in
        X0Z1=X0Z1in
        X1Z0=X1Z0in
        X1Z1=X1Z1in
        X0Y0=X0Y0in
        X1Y0=X1Y0in
        X0Y1=X0Y1in
        X1Y1=X1Y1in

c 3d nodes
        X0Y0Z0=X0Y0Z0in
        X1Y0Z0=X1Y0Z0in
        X0Y1Z0=X0Y1Z0in
        X1Y1Z0=X1Y1Z0in
        X0Y0Z1=X0Y0Z1in
        X1Y0Z1=X1Y0Z1in
        X0Y1Z1=X0Y1Z1in
        X1Y1Z1=X1Y1Z1in

      return
      end
c***********************************************************************
c***********************************************************************
      subroutine stufcartbdrycond3d(
     &  FLOWin, XFLOWin, YFLOWin, ZFLOWin,
     &  REFLECTin, XREFLECTin, YREFLECTin, ZREFLECTin, 
     &  DIRICHLETin, XDIRICHLETin, YDIRICHLETin, ZDIRICHLETin,
     &  NEUMANNin, XNEUMANNin, YNEUMANNin, ZNEUMANNin)
      implicit none
      integer
     &  FLOWin, XFLOWin, YFLOWin, ZFLOWin,
     &  REFLECTin, XREFLECTin, YREFLECTin, ZREFLECTin,
     &  DIRICHLETin, XDIRICHLETin, YDIRICHLETin, ZDIRICHLETin,
     &  NEUMANNin, XNEUMANNin, YNEUMANNin, ZNEUMANNin
include(FORTDIR/appu_cartbdryparams3d.i)dnl

        FLOW=FLOWin
        XFLOW=XFLOWin
        YFLOW=YFLOWin
        ZFLOW=ZFLOWin

        REFLECT=REFLECTin
        XREFLECT=XREFLECTin
        YREFLECT=YREFLECTin
        ZREFLECT=ZREFLECTin

        DIRICHLET=DIRICHLETin
        XDIRICHLET=XDIRICHLETin
        YDIRICHLET=YDIRICHLETin
        ZDIRICHLET=ZDIRICHLETin

        NEUMANN=NEUMANNin
        XNEUMANN=XNEUMANNin
        YNEUMANN=YNEUMANNin
        ZNEUMANN=ZNEUMANNin

      return
      end

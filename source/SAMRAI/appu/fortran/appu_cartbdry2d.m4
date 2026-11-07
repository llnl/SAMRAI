c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for Cartesian 2d boundary conditions.
c
define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl
c
c general boundary condition cases
c
define(setvalues_to_ict,`dnl
         do k=0,arrdepth-1
           arrdata(ic0,ic1,k)=arrdata(ict0,ict1,k)
         enddo
')dnl
define(setvalues_to_ict_reflect,`dnl
         do k=0,arrdepth-1
           arrdata(ic0,ic1,k)=arrdata(ict0,ict1,k)
         enddo
         arrdata(ic0,ic1,$1)= -arrdata(ic0,ic1,$1)
')dnl
define(setvalues_to_dirichlet,`dnl
         do k=0,arrdepth-1
           arrdata(ic0,ic1,k)=edge_values(k,$1)
         enddo
')dnl
define(setvalues_to_neumann,`dnl
         do k=0,arrdepth-1
           arrdata(ic0,ic1,k)=arrdata(ict0,ict1,k) + $2
     &       * edge_values(k,$3) * dble(ic$1-ict$1) * dx($1)
         enddo
')dnl
c
c node boundary conditions
c
define(do_node_ict,`dnl
             do ic1=ibeg1,iend1
               ict1 = $2
               do ic0=ibeg0,iend0
                 ict0 = $1
setvalues_to_ict()dnl
                enddo
             enddo
')dnl
define(do_node_ict_reflect,`dnl
             do ic1=ibeg1,iend1
               ict1 = $2
               do ic0=ibeg0,iend0
                 ict0 = $1
setvalues_to_ict_reflect($3)dnl
               enddo
             enddo
')dnl
define(do_node_dirichlet,`dnl
             do ic1=ibeg1,iend1
               do ic0=ibeg0,iend0
setvalues_to_dirichlet($1)dnl
               enddo
             enddo
')dnl
define(do_node_neumann,`dnl
             do ic1=ibeg1,iend1
               ict1 = $2
                do ic0=ibeg0,iend0
                 ict0 = $1
setvalues_to_neumann($3,$4,$5)dnl
               enddo
             enddo
')dnl
define(do_bdry_node,`dnl
           if (bdry_cond.eq.XDIRICHLET) then
do_node_dirichlet(`edge_loc0')dnl
           else if (bdry_cond.eq.YDIRICHLET) then
do_node_dirichlet(`edge_loc1')dnl
           else if (bdry_cond.eq.XNEUMANN) then
do_node_neumann(`ipivot0',`ic1',0,`dirsign0',`edge_loc0')dnl
           else if (bdry_cond.eq.YNEUMANN) then
do_node_neumann(`ic0',`ipivot1',1,`dirsign1',`edge_loc1')dnl
           else if (bdry_cond.eq.XFLOW) then
do_node_ict(`ipivot0',`ic1')dnl
           else if (bdry_cond.eq.YFLOW) then
do_node_ict(`ic0',`ipivot1')dnl
           else if (bdry_cond.eq.XREFLECT) then
do_node_ict_reflect(`ipivot0',`ic1',0)dnl
           else if (bdry_cond.eq.YREFLECT) then
do_node_ict_reflect(`ic0',`ipivot1',1)dnl
           else
             write(6,*) "INVALID NODE bdry_cond in getcartnodebdry2d"
             write(6,*) "bdry_loc = ",bdry_loc,
     &                  "bdry_cond = ",bdry_cond
           endif
')dnl
c
c edge boundary conditions
c
define(do_bdry_edge,`dnl
           if (bdry_cond.eq.DIRICHLET) then
             do ic1=ibeg1,iend1
               do ic0=ibeg0,iend0
setvalues_to_dirichlet(`edge_loc')dnl
               enddo
             enddo
           else if (bdry_cond.eq.NEUMANN) then
             do ic1=ibeg1,iend1
               ict1 = $3
               do ic0=ibeg0,iend0
                 ict0 = $2
setvalues_to_neumann($1,$4,`edge_loc')dnl
               enddo
             enddo
           else if (bdry_cond.eq.FLOW) then
             do ic1=ibeg1,iend1
               ict1 = $3
               do ic0=ibeg0,iend0
                 ict0 = $2
setvalues_to_ict()dnl
               enddo
             enddo
           else if (bdry_cond.eq.REFLECT) then
             do ic1=ibeg1,iend1
               ict1 = $3
               do ic0=ibeg0,iend0
                 ict0 = $2
setvalues_to_ict_reflect($1)dnl
               enddo
             enddo
          else
             write(6,*) "INVALID EDGE bdry_cond in getcartedgebdry2d"
             write(6,*) "bdry_loc = ",bdry_loc,
     &                  "bdry_cond = ",bdry_cond
          endif
')dnl
c***********************************************************************
c***********************************************************************
      subroutine getcartedgebdry2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  ibeg0,iend0,ibeg1,iend1,
     &  ngc0,ngc1,
     &  dx,
     &  bdry_loc,
     &  bdry_cond,
     &  edge_values,
     &  arrdata,
     &  arrdepth)
c***********************************************************************
      implicit none
include(FORTDIR/appu_cartbdryparams2d.i)dnl
c***********************************************************************
      integer ifirst0,ilast0,ifirst1,ilast1
      integer ibeg0,iend0,ibeg1,iend1
      integer ngc0,ngc1
c
      REAL dx(0:NDIM-1)
c
      integer bdry_loc,
     &        arrdepth, 
     &        bdry_cond
      REAL
     &        edge_values(0:arrdepth-1,0:2*NDIM-1)
c
      REAL
     &  arrdata(CELL2dVECG(ifirst,ilast,ngc),0:arrdepth-1)
c
      integer ic1,ic0,ict0,ict1
      integer k
      integer ipivot0,ipivot1
      integer edge_loc
      REAL    dirsign0,dirsign1
c
c***********************************************************************
c   bdry_loc   index     bdry_loc  index   
c       0     (-1, 0)       2     (0,-1)  
c       1     ( 1, 0)       3     (0, 1)  
c***********************************************************************
c***********************************************************************
      edge_loc = bdry_loc
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
do_bdry_edge(0,`ipivot0',`ic1',`dirsign0')dnl
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
do_bdry_edge(1,`ic0',`ipivot1',`dirsign1')dnl
      endif
c
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine getcartnodebdry2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  ibeg0,iend0,ibeg1,iend1,
     &  ngc0,ngc1,
     &  dx,
     &  bdry_loc,
     &  bdry_cond,
     &  edge_values,
     &  arrdata,
     &  arrdepth)
c***********************************************************************
      implicit none
include(FORTDIR/appu_cartbdryparams2d.i)dnl
c***********************************************************************
      integer ifirst0,ilast0,ifirst1,ilast1
      integer ibeg0,iend0,ibeg1,iend1
      integer ngc0,ngc1
c
      REAL dx(0:NDIM-1)
c
      integer bdry_loc,
     &        arrdepth,
     &        bdry_cond
      REAL
     &        edge_values(0:arrdepth-1,0:2*NDIM-1)
c
      REAL
     &  arrdata(CELL2dVECG(ifirst,ilast,ngc),0:arrdepth-1)
c
      integer ic1,ic0,ict0,ict1
      integer k
      integer ipivot0,ipivot1
      integer edge_loc0,edge_loc1
      REAL    dirsign0,dirsign1
c
c***********************************************************************
c***********************************************************************
c    bdry_loc   index 
c       0     (-1,-1)
c       1     ( 1,-1)
c       2     (-1, 1)
c       3     ( 1, 1)
c***********************************************************************
c***********************************************************************
      if (bdry_loc.eq.X0Y0) then
         edge_loc0 = XLEFT
         ipivot0   = ifirst0
         dirsign0  = -1.d0
         edge_loc1 = YLEFT
         ipivot1   = ifirst1
         dirsign1  = -1.d0
      else if (bdry_loc.eq.X1Y0) then
         edge_loc0 = XRIGHT
         ipivot0   = ilast0
         dirsign0  = 1.d0
         edge_loc1 = YLEFT
         ipivot1   = ifirst1
         dirsign1  = -1.d0
      else if (bdry_loc.eq.X0Y1) then
         edge_loc0 = XLEFT
         ipivot0   = ifirst0
         dirsign0  = -1.d0
         edge_loc1 = YRIGHT
         ipivot1   = ilast1
         dirsign1  = 1.d0
      else if (bdry_loc.eq.X1Y1) then
         edge_loc0 = XRIGHT
         ipivot0   = ilast0
         dirsign0  = 1.d0
         edge_loc1 = YRIGHT
         ipivot1   = ilast1
         dirsign1  = 1.d0
      endif
do_bdry_node()dnl
c
      return
      end

c***********************************************************************
c***********************************************************************
      subroutine stufcartbdryloc2d(
     &  XLEFTin, XRIGHTin, YLEFTin, YRIGHTin,
     &  X0Y0in, X1Y0in, X0Y1in, X1Y1in)
      implicit none
      integer
     &  XLEFTin, XRIGHTin, YLEFTin, YRIGHTin,
     &  X0Y0in, X1Y0in, X0Y1in, X1Y1in
include(FORTDIR/appu_cartbdryparams2d.i)dnl

c 2d edges
        XLEFT=XLEFTin
        XRIGHT=XRIGHTin
        YLEFT=YLEFTin
        YRIGHT=YRIGHTin

c 2d nodes
        X0Y0=X0Y0in
        X1Y0=X1Y0in
        X0Y1=X0Y1in
        X1Y1=X1Y1in

      return
      end
c***********************************************************************
c***********************************************************************
      subroutine stufcartbdrycond2d(
     &  FLOWin, XFLOWin, YFLOWin, 
     &  REFLECTin, XREFLECTin, YREFLECTin,
     &  DIRICHLETin, XDIRICHLETin, YDIRICHLETin,
     &  NEUMANNin, XNEUMANNin, YNEUMANNin)
      implicit none
      integer
     &  FLOWin, XFLOWin, YFLOWin, 
     &  REFLECTin, XREFLECTin, YREFLECTin,
     &  DIRICHLETin, XDIRICHLETin, YDIRICHLETin,
     &  NEUMANNin, XNEUMANNin, YNEUMANNin
include(FORTDIR/appu_cartbdryparams2d.i)dnl

        FLOW=FLOWin
        XFLOW=XFLOWin
        YFLOW=YFLOWin

        REFLECT=REFLECTin
        XREFLECT=XREFLECTin
        YREFLECT=YREFLECTin

        DIRICHLET=DIRICHLETin
        XDIRICHLET=XDIRICHLETin
        YDIRICHLET=YDIRICHLETin

        NEUMANN=NEUMANNin
        XNEUMANN=XNEUMANNin
        YNEUMANN=YNEUMANNin

      return
      end

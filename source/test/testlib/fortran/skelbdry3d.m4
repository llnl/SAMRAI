c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   FORTRAN routines for general 3d boundary condition cases,
c                node boundary conditions, and edge boundary conditions
c
c***********************************************************************
c***********************************************************************
      subroutine getskelfacebdry3d(
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
c
c  3d boundary constant common blocks
c
      common/skelbdrylocparams3d/
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,ZLEFT,ZRIGHT, 
     &  Y0Z0,Y1Z0,Y0Z1,Y1Z1,
     &  X0Z0,X0Z1,X1Z0,X1Z1,
     &  X0Y0,X1Y0,X0Y1,X1Y1,
     &  X0Y0Z0,X1Y0Z0,X0Y1Z0,X1Y1Z0,
     &  X0Y0Z1,X1Y0Z1,X0Y1Z1,X1Y1Z1
      integer
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,ZLEFT,ZRIGHT,
     &  Y0Z0,Y1Z0,Y0Z1,Y1Z1,
     &  X0Z0,X0Z1,X1Z0,X1Z1,
     &  X0Y0,X1Y0,X0Y1,X1Y1,
     &  X0Y0Z0,X1Y0Z0,X0Y1Z0,X1Y1Z0,
     &  X0Y0Z1,X1Y0Z1,X0Y1Z1,X1Y1Z1  
c
c
      common/skelbdrycondparams3d/
     &  FLOW,XFLOW,YFLOW,ZFLOW,
     &  REFLECT,XREFLECT,YREFLECT,ZREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,ZDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN,ZNEUMANN
      integer
     &  FLOW,XFLOW,YFLOW,ZFLOW,
     &  REFLECT,XREFLECT,YREFLECT,ZREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,ZDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN,ZNEUMANN
c***********************************************************************
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      integer ngc0,ngc1,ngc2
c
      double precision dx(0:3-1)
c
      integer bdry_loc,
     &        arrdepth, 
     &        bdry_cond
      double precision
     &        face_values(0:arrdepth-1,0:2*3-1)
c
      double precision
     &  arrdata(ifirst0-ngc0:ilast0+ngc0,
     &          ifirst1-ngc1:ilast1+ngc1,
     &          ifirst2-ngc2:ilast2+ngc2,0:arrdepth-1)
c
      integer ic2,ic1,ic0,ict0,ict1,ict2
      integer k
      integer ipivot0,ipivot1,ipivot2
      integer face_loc
      double precision    dirsign0,dirsign1,dirsign2
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
           if (bdry_cond.eq.DIRICHLET) then
             do ic2=ibeg2,iend2
               do ic1=ibeg1,iend1
                 do ic0=ibeg0,iend0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=face_values(k,face_loc)
         enddo
                 enddo
               enddo
             enddo
           else if (bdry_cond.eq.NEUMANN) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ipivot0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k) + dirsign0
     &       * face_values(k,face_loc) * dble(ic0-ict0) * dx(0)
         enddo
                 enddo
               enddo
             enddo
           else if (bdry_cond.eq.FLOW) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ipivot0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
                 enddo
               enddo
             enddo
           else if (bdry_cond.eq.REFLECT) then
             do ic2=ibeg2,iend2
              ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ipivot0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
         arrdata(ic0,ic1,ic2,0)= -arrdata(ic0,ic1,ic2,0)
                 enddo
               enddo
             enddo
          else
             write(6,*) "INVALID FACE bdry_cond in getskelfacebdry3d"
             write(6,*) "bdry_loc = ",bdry_loc,
     &                  "bdry_cond = ",bdry_cond
          endif
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
           if (bdry_cond.eq.DIRICHLET) then
             do ic2=ibeg2,iend2
               do ic1=ibeg1,iend1
                 do ic0=ibeg0,iend0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=face_values(k,face_loc)
         enddo
                 enddo
               enddo
             enddo
           else if (bdry_cond.eq.NEUMANN) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ipivot1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k) + dirsign1
     &       * face_values(k,face_loc) * dble(ic1-ict1) * dx(1)
         enddo
                 enddo
               enddo
             enddo
           else if (bdry_cond.eq.FLOW) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ipivot1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
                 enddo
               enddo
             enddo
           else if (bdry_cond.eq.REFLECT) then
             do ic2=ibeg2,iend2
              ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ipivot1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
         arrdata(ic0,ic1,ic2,1)= -arrdata(ic0,ic1,ic2,1)
                 enddo
               enddo
             enddo
          else
             write(6,*) "INVALID FACE bdry_cond in getskelfacebdry3d"
             write(6,*) "bdry_loc = ",bdry_loc,
     &                  "bdry_cond = ",bdry_cond
          endif
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
           if (bdry_cond.eq.DIRICHLET) then
             do ic2=ibeg2,iend2
               do ic1=ibeg1,iend1
                 do ic0=ibeg0,iend0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=face_values(k,face_loc)
         enddo
                 enddo
               enddo
             enddo
           else if (bdry_cond.eq.NEUMANN) then
             do ic2=ibeg2,iend2
               ict2 = ipivot2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k) + dirsign2
     &       * face_values(k,face_loc) * dble(ic2-ict2) * dx(2)
         enddo
                 enddo
               enddo
             enddo
           else if (bdry_cond.eq.FLOW) then
             do ic2=ibeg2,iend2
               ict2 = ipivot2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
                 enddo
               enddo
             enddo
           else if (bdry_cond.eq.REFLECT) then
             do ic2=ibeg2,iend2
              ict2 = ipivot2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
         arrdata(ic0,ic1,ic2,2)= -arrdata(ic0,ic1,ic2,2)
                 enddo
               enddo
             enddo
          else
             write(6,*) "INVALID FACE bdry_cond in getskelfacebdry3d"
             write(6,*) "bdry_loc = ",bdry_loc,
     &                  "bdry_cond = ",bdry_cond
          endif
      endif
c
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine getskeledgebdry3d(
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
c
c  3d boundary constant common blocks
c
      common/skelbdrylocparams3d/
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,ZLEFT,ZRIGHT, 
     &  Y0Z0,Y1Z0,Y0Z1,Y1Z1,
     &  X0Z0,X0Z1,X1Z0,X1Z1,
     &  X0Y0,X1Y0,X0Y1,X1Y1,
     &  X0Y0Z0,X1Y0Z0,X0Y1Z0,X1Y1Z0,
     &  X0Y0Z1,X1Y0Z1,X0Y1Z1,X1Y1Z1
      integer
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,ZLEFT,ZRIGHT,
     &  Y0Z0,Y1Z0,Y0Z1,Y1Z1,
     &  X0Z0,X0Z1,X1Z0,X1Z1,
     &  X0Y0,X1Y0,X0Y1,X1Y1,
     &  X0Y0Z0,X1Y0Z0,X0Y1Z0,X1Y1Z0,
     &  X0Y0Z1,X1Y0Z1,X0Y1Z1,X1Y1Z1  
c
c
      common/skelbdrycondparams3d/
     &  FLOW,XFLOW,YFLOW,ZFLOW,
     &  REFLECT,XREFLECT,YREFLECT,ZREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,ZDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN,ZNEUMANN
      integer
     &  FLOW,XFLOW,YFLOW,ZFLOW,
     &  REFLECT,XREFLECT,YREFLECT,ZREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,ZDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN,ZNEUMANN
c***********************************************************************
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      integer ngc0,ngc1,ngc2
c
      double precision dx(0:3-1)
c
      integer bdry_loc,
     &        arrdepth,
     &        bdry_cond
      double precision
     &        face_values(0:arrdepth-1,0:2*3-1)
c
      double precision
     &  arrdata(ifirst0-ngc0:ilast0+ngc0,
     &          ifirst1-ngc1:ilast1+ngc1,
     &          ifirst2-ngc2:ilast2+ngc2,0:arrdepth-1)
c
      integer ic2,ic1,ic0,ict0,ict1,ict2
      integer k
      integer ipivot0,ipivot1,ipivot2
      integer face_loc0,face_loc1,face_loc2
      double precision    dirsign0,dirsign1,dirsign2
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
           if (bdry_cond.eq.XDIRICHLET) then
             do ic2=ibeg2,iend2
               do ic1=ibeg1,iend1
                 do ic0=ibeg0,iend0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=face_values(k,face_loc0)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.YDIRICHLET) then
             do ic2=ibeg2,iend2
               do ic1=ibeg1,iend1
                 do ic0=ibeg0,iend0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=face_values(k,face_loc1)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.ZDIRICHLET) then
             do ic2=ibeg2,iend2
               do ic1=ibeg1,iend1
                 do ic0=ibeg0,iend0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=face_values(k,face_loc2)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.XNEUMANN) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ipivot0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k) + dirsign0
     &       * face_values(k,face_loc0) * dble(ic0-ict0) * dx(0)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.YNEUMANN) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ipivot1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k) + dirsign1
     &       * face_values(k,face_loc1) * dble(ic1-ict1) * dx(1)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.ZNEUMANN) then
             do ic2=ibeg2,iend2
               ict2 = ipivot2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k) + dirsign2
     &       * face_values(k,face_loc2) * dble(ic2-ict2) * dx(2)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.XFLOW) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ipivot0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.YFLOW) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ipivot1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.ZFLOW) then
             do ic2=ibeg2,iend2
               ict2 = ipivot2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.XREFLECT) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ipivot0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
         arrdata(ic0,ic1,ic2,0)= -arrdata(ic0,ic1,ic2,0)
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.YREFLECT) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ipivot1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
         arrdata(ic0,ic1,ic2,1)= -arrdata(ic0,ic1,ic2,1)
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.ZREFLECT) then
             do ic2=ibeg2,iend2
               ict2 = ipivot2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
         arrdata(ic0,ic1,ic2,2)= -arrdata(ic0,ic1,ic2,2)
                   enddo
                enddo
             enddo
           else
             write(6,*) "INVALID EDGE bdry_cond in getskeledgebdry3d"
             write(6,*) "bdry_loc = ",bdry_loc,
     &                  "bdry_cond = ",bdry_cond
           endif
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
           if (bdry_cond.eq.XDIRICHLET) then
             do ic2=ibeg2,iend2
               do ic1=ibeg1,iend1
                 do ic0=ibeg0,iend0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=face_values(k,face_loc0)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.YDIRICHLET) then
             do ic2=ibeg2,iend2
               do ic1=ibeg1,iend1
                 do ic0=ibeg0,iend0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=face_values(k,face_loc1)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.ZDIRICHLET) then
             do ic2=ibeg2,iend2
               do ic1=ibeg1,iend1
                 do ic0=ibeg0,iend0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=face_values(k,face_loc2)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.XNEUMANN) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ipivot0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k) + dirsign0
     &       * face_values(k,face_loc0) * dble(ic0-ict0) * dx(0)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.YNEUMANN) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ipivot1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k) + dirsign1
     &       * face_values(k,face_loc1) * dble(ic1-ict1) * dx(1)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.ZNEUMANN) then
             do ic2=ibeg2,iend2
               ict2 = ipivot2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k) + dirsign2
     &       * face_values(k,face_loc2) * dble(ic2-ict2) * dx(2)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.XFLOW) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ipivot0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.YFLOW) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ipivot1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.ZFLOW) then
             do ic2=ibeg2,iend2
               ict2 = ipivot2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.XREFLECT) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ipivot0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
         arrdata(ic0,ic1,ic2,0)= -arrdata(ic0,ic1,ic2,0)
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.YREFLECT) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ipivot1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
         arrdata(ic0,ic1,ic2,1)= -arrdata(ic0,ic1,ic2,1)
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.ZREFLECT) then
             do ic2=ibeg2,iend2
               ict2 = ipivot2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
         arrdata(ic0,ic1,ic2,2)= -arrdata(ic0,ic1,ic2,2)
                   enddo
                enddo
             enddo
           else
             write(6,*) "INVALID EDGE bdry_cond in getskeledgebdry3d"
             write(6,*) "bdry_loc = ",bdry_loc,
     &                  "bdry_cond = ",bdry_cond
           endif
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
           if (bdry_cond.eq.XDIRICHLET) then
             do ic2=ibeg2,iend2
               do ic1=ibeg1,iend1
                 do ic0=ibeg0,iend0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=face_values(k,face_loc0)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.YDIRICHLET) then
             do ic2=ibeg2,iend2
               do ic1=ibeg1,iend1
                 do ic0=ibeg0,iend0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=face_values(k,face_loc1)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.ZDIRICHLET) then
             do ic2=ibeg2,iend2
               do ic1=ibeg1,iend1
                 do ic0=ibeg0,iend0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=face_values(k,face_loc2)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.XNEUMANN) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ipivot0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k) + dirsign0
     &       * face_values(k,face_loc0) * dble(ic0-ict0) * dx(0)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.YNEUMANN) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ipivot1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k) + dirsign1
     &       * face_values(k,face_loc1) * dble(ic1-ict1) * dx(1)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.ZNEUMANN) then
             do ic2=ibeg2,iend2
               ict2 = ipivot2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k) + dirsign2
     &       * face_values(k,face_loc2) * dble(ic2-ict2) * dx(2)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.XFLOW) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ipivot0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.YFLOW) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ipivot1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.ZFLOW) then
             do ic2=ibeg2,iend2
               ict2 = ipivot2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.XREFLECT) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ipivot0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
         arrdata(ic0,ic1,ic2,0)= -arrdata(ic0,ic1,ic2,0)
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.YREFLECT) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ipivot1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
         arrdata(ic0,ic1,ic2,1)= -arrdata(ic0,ic1,ic2,1)
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.ZREFLECT) then
             do ic2=ibeg2,iend2
               ict2 = ipivot2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
         arrdata(ic0,ic1,ic2,2)= -arrdata(ic0,ic1,ic2,2)
                   enddo
                enddo
             enddo
           else
             write(6,*) "INVALID EDGE bdry_cond in getskeledgebdry3d"
             write(6,*) "bdry_loc = ",bdry_loc,
     &                  "bdry_cond = ",bdry_cond
           endif
      endif
c
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine getskelnodebdry3d(
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
c
c  3d boundary constant common blocks
c
      common/skelbdrylocparams3d/
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,ZLEFT,ZRIGHT, 
     &  Y0Z0,Y1Z0,Y0Z1,Y1Z1,
     &  X0Z0,X0Z1,X1Z0,X1Z1,
     &  X0Y0,X1Y0,X0Y1,X1Y1,
     &  X0Y0Z0,X1Y0Z0,X0Y1Z0,X1Y1Z0,
     &  X0Y0Z1,X1Y0Z1,X0Y1Z1,X1Y1Z1
      integer
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,ZLEFT,ZRIGHT,
     &  Y0Z0,Y1Z0,Y0Z1,Y1Z1,
     &  X0Z0,X0Z1,X1Z0,X1Z1,
     &  X0Y0,X1Y0,X0Y1,X1Y1,
     &  X0Y0Z0,X1Y0Z0,X0Y1Z0,X1Y1Z0,
     &  X0Y0Z1,X1Y0Z1,X0Y1Z1,X1Y1Z1  
c
c
      common/skelbdrycondparams3d/
     &  FLOW,XFLOW,YFLOW,ZFLOW,
     &  REFLECT,XREFLECT,YREFLECT,ZREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,ZDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN,ZNEUMANN
      integer
     &  FLOW,XFLOW,YFLOW,ZFLOW,
     &  REFLECT,XREFLECT,YREFLECT,ZREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,ZDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN,ZNEUMANN
c***********************************************************************
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      integer ngc0,ngc1,ngc2
c
      double precision dx(0:3-1)
c
      integer bdry_loc,
     &        arrdepth,
     &        bdry_cond
      double precision
     &        face_values(0:arrdepth-1,0:2*3-1)
c
      double precision
     &  arrdata(ifirst0-ngc0:ilast0+ngc0,
     &          ifirst1-ngc1:ilast1+ngc1,
     &          ifirst2-ngc2:ilast2+ngc2,0:arrdepth-1)
c
      integer ic2,ic1,ic0,ict0,ict1,ict2
      integer k
      integer ipivot0,ipivot1,ipivot2
      integer face_loc0,face_loc1,face_loc2
      double precision    dirsign0,dirsign1,dirsign2
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
           if (bdry_cond.eq.XDIRICHLET) then
             do ic2=ibeg2,iend2
               do ic1=ibeg1,iend1
                 do ic0=ibeg0,iend0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=face_values(k,face_loc0)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.YDIRICHLET) then
             do ic2=ibeg2,iend2
               do ic1=ibeg1,iend1
                 do ic0=ibeg0,iend0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=face_values(k,face_loc1)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.ZDIRICHLET) then
             do ic2=ibeg2,iend2
               do ic1=ibeg1,iend1
                 do ic0=ibeg0,iend0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=face_values(k,face_loc2)
         enddo
                   enddo
                enddo
             enddo
          else if (bdry_cond.eq.XNEUMANN) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ipivot0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k) + dirsign0
     &       * face_values(k,face_loc0) * dble(ic0-ict0) * dx(0)
         enddo
                   enddo
                enddo
             enddo
          else if (bdry_cond.eq.YNEUMANN) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ipivot1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k) + dirsign1
     &       * face_values(k,face_loc1) * dble(ic1-ict1) * dx(1)
         enddo
                   enddo
                enddo
             enddo
          else if (bdry_cond.eq.ZNEUMANN) then
             do ic2=ibeg2,iend2
               ict2 = ipivot2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k) + dirsign2
     &       * face_values(k,face_loc2) * dble(ic2-ict2) * dx(2)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.XFLOW) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ipivot0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.YFLOW) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ipivot1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.ZFLOW) then
             do ic2=ibeg2,iend2
               ict2 = ipivot2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.XREFLECT) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ipivot0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
         arrdata(ic0,ic1,ic2,0)= -arrdata(ic0,ic1,ic2,0)
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.YREFLECT) then
             do ic2=ibeg2,iend2
               ict2 = ic2
               do ic1=ibeg1,iend1
                 ict1 = ipivot1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
         arrdata(ic0,ic1,ic2,1)= -arrdata(ic0,ic1,ic2,1)
                   enddo
                enddo
             enddo
           else if (bdry_cond.eq.ZREFLECT) then
             do ic2=ibeg2,iend2
               ict2 = ipivot2
               do ic1=ibeg1,iend1
                 ict1 = ic1
                 do ic0=ibeg0,iend0
                   ict0 = ic0
         do k=0,arrdepth-1
           arrdata(ic0,ic1,ic2,k)=arrdata(ict0,ict1,ict2,k)
         enddo
         arrdata(ic0,ic1,ic2,2)= -arrdata(ic0,ic1,ic2,2)
                   enddo
                enddo
             enddo
           else
             write(6,*) "INVALID NODE bdry_cond in getskelnodebdry3d"
             write(6,*) "bdry_loc = ",bdry_loc,
     &                  "bdry_cond = ",bdry_cond
           endif
c
      return
      end

c***********************************************************************
c***********************************************************************
      subroutine stufskelbdryloc3d(
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
c
c  3d boundary constant common blocks
c
      common/skelbdrylocparams3d/
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,ZLEFT,ZRIGHT, 
     &  Y0Z0,Y1Z0,Y0Z1,Y1Z1,
     &  X0Z0,X0Z1,X1Z0,X1Z1,
     &  X0Y0,X1Y0,X0Y1,X1Y1,
     &  X0Y0Z0,X1Y0Z0,X0Y1Z0,X1Y1Z0,
     &  X0Y0Z1,X1Y0Z1,X0Y1Z1,X1Y1Z1
      integer
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,ZLEFT,ZRIGHT,
     &  Y0Z0,Y1Z0,Y0Z1,Y1Z1,
     &  X0Z0,X0Z1,X1Z0,X1Z1,
     &  X0Y0,X1Y0,X0Y1,X1Y1,
     &  X0Y0Z0,X1Y0Z0,X0Y1Z0,X1Y1Z0,
     &  X0Y0Z1,X1Y0Z1,X0Y1Z1,X1Y1Z1  
c
c
      common/skelbdrycondparams3d/
     &  FLOW,XFLOW,YFLOW,ZFLOW,
     &  REFLECT,XREFLECT,YREFLECT,ZREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,ZDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN,ZNEUMANN
      integer
     &  FLOW,XFLOW,YFLOW,ZFLOW,
     &  REFLECT,XREFLECT,YREFLECT,ZREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,ZDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN,ZNEUMANN

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
      subroutine stufskelbdrycond3d(
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
c
c  3d boundary constant common blocks
c
      common/skelbdrylocparams3d/
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,ZLEFT,ZRIGHT, 
     &  Y0Z0,Y1Z0,Y0Z1,Y1Z1,
     &  X0Z0,X0Z1,X1Z0,X1Z1,
     &  X0Y0,X1Y0,X0Y1,X1Y1,
     &  X0Y0Z0,X1Y0Z0,X0Y1Z0,X1Y1Z0,
     &  X0Y0Z1,X1Y0Z1,X0Y1Z1,X1Y1Z1
      integer
     &  XLEFT,XRIGHT,YLEFT,YRIGHT,ZLEFT,ZRIGHT,
     &  Y0Z0,Y1Z0,Y0Z1,Y1Z1,
     &  X0Z0,X0Z1,X1Z0,X1Z1,
     &  X0Y0,X1Y0,X0Y1,X1Y1,
     &  X0Y0Z0,X1Y0Z0,X0Y1Z0,X1Y1Z0,
     &  X0Y0Z1,X1Y0Z1,X0Y1Z1,X1Y1Z1  
c
c
      common/skelbdrycondparams3d/
     &  FLOW,XFLOW,YFLOW,ZFLOW,
     &  REFLECT,XREFLECT,YREFLECT,ZREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,ZDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN,ZNEUMANN
      integer
     &  FLOW,XFLOW,YFLOW,ZFLOW,
     &  REFLECT,XREFLECT,YREFLECT,ZREFLECT,
     &  DIRICHLET,XDIRICHLET,YDIRICHLET,ZDIRICHLET,
     &  NEUMANN,XNEUMANN,YNEUMANN,ZNEUMANN

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

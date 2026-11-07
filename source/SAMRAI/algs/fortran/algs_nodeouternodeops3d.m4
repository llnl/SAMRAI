c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for summing outernode data with other node or
c                outernode data.
c
define(NDIM,3)dnl
include(PDAT_FORTDIR/pdat_m4arrdim3d.i)dnl

define(do_cfsum_along_patch_face,`dnl
      do id = 1, depth
         do ic$2=cilo$2_loop,cihi$2_loop
            if$2=ic$2*ratio($2)
            do ic$1=cilo$1_loop,cihi$1_loop
               if$1=ic$1*ratio($1)
               fnode(if0,if1,if2,id) = fnode(if0,if1,if2,id) +
     &                                 couternode$3(ic$1,ic$2,id)
            enddo
         enddo
      enddo
')dnl

define(interp_interior_loop_begin,`dnl
         do ic$2=cilo$2_loop,cihi$2_loop
            iflower$2=ic$2*ratio($2)
            ifupper$2=iflower$2+ratio($2)
            do ir$2=1,ratio($2)-1
               x$2 = dble(ir$2)/dratio$2
               do ic$1=cilo$1_loop,cihi$1_loop
                  iflower$1=ic$1*ratio($1)
                  ifupper$1=iflower$1+ratio($1)
                  do ir$1=1,ratio($1)-1
                     x$1 = dble(ir$1)/dratio$1
')dnl
define(interp_interior_loop_end,`dnl
                  enddo
               enddo
            enddo
         enddo
')dnl

define(do_interp_interior0,`dnl
interp_interior_loop_begin(1,2)dnl
                     fnode(if0,iflower1+ir1,iflower2+ir2,id) =
     &               (fnode(if0,iflower1,iflower2,id) * (one-x1) +
     &                fnode(if0,ifupper1,iflower2,id) * x1) * (one-x2) +
     &               (fnode(if0,iflower1,ifupper2,id) * (one-x1) +
     &                fnode(if0,ifupper1,ifupper2,id) * x1) * x2
interp_interior_loop_end()dnl
')dnl

define(do_interp_interior1,`dnl
interp_interior_loop_begin(0,2)dnl
                     fnode(iflower0+ir0,if1,iflower2+ir2,id) =
     &               (fnode(iflower0,if1,iflower2,id) * (one-x0) +
     &                fnode(ifupper0,if1,iflower2,id) * x0) * (one-x2) +
     &               (fnode(iflower0,if1,ifupper2,id) * (one-x0) +
     &                fnode(ifupper0,if1,ifupper2,id) * x0) * x2
interp_interior_loop_end()dnl
')dnl

define(do_interp_interior2,`dnl
interp_interior_loop_begin(0,1)dnl
                     fnode(iflower0+ir0,iflower1+ir1,if2,id) =
     &               (fnode(iflower0,iflower1,if2,id) * (one-x0) +
     &                fnode(ifupper0,iflower1,if2,id) * x0) * (one-x1) +
     &               (fnode(iflower0,ifupper1,if2,id) * (one-x0) +
     &                fnode(ifupper0,ifupper1,if2,id) * x0) * x1
interp_interior_loop_end()dnl
')dnl

define(do_interp_lines0,`dnl
         do ic$1=cilo$1_loop,cihi$1_loop+1
            if$1=ic$1*ratio($1)
            do ic0=cilo0_loop,cihi0_loop
               iflower0=ic0*ratio(0)
               ifupper0=iflower0+ratio(0)
               do ir0=1,ratio(0)-1
                  x0 = dble(ir0)/dratio0
                  fnode(iflower0+ir0,if1,if2,id) =
     &               fnode(ifupper0,if1,if2,id) * x0
     &           +   fnode(iflower0,if1,if2,id) * (one-x0)
               enddo
            enddo
         enddo
')dnl

define(do_interp_lines1,`dnl
         do ic1=cilo1_loop,cihi1_loop
            iflower1=ic1*ratio(1)
            ifupper1=iflower1+ratio(1)
            do ir1=1,ratio(1)-1
               x1 = dble(ir1)/dratio1
               do ic$1=cilo$1_loop,cihi$1_loop+1
                  if$1=ic$1*ratio($1)
                  fnode(if0,iflower1+ir1,if2,id) =
     &               fnode(if0,ifupper1,if2,id) * x1
     &           +   fnode(if0,iflower1,if2,id) * (one-x1)
               enddo
            enddo
         enddo
')dnl

define(do_interp_lines2,`dnl
         do ic2=cilo2_loop,cihi2_loop
            iflower2=ic2*ratio(2)
            ifupper2=iflower2+ratio(2)
            do ir2=1,ratio(2)-1
               x2 = dble(ir2)/dratio2
               do ic$1=cilo$1_loop,cihi$1_loop+1
                  if$1=ic$1*ratio($1)
                  fnode(if0,if1,iflower2+ir2,id) =
     &               fnode(if0,if1,ifupper2,id) * x2
     &           +   fnode(if0,if1,iflower2,id) * (one-x2)
               enddo
            enddo
         enddo
')dnl

c
c***********************************************************************
c Sum node data on fine patch with outernode data on coarse patch
c and store back in node data on fine patch.
c i.e. fine_node_data += coarse_outernode_data
c***********************************************************************
c
      subroutine nodeouternodesum3d(
     &  filo0,filo1,filo2,
     &  fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,
     &  cihi0,cihi1,cihi2,
     &  ratio,
     &  depth,
     &  ngc0,ngc1,ngc2,
     &  fnode,
     &  couternodelower0,couternodeupper0,
     &  couternodelower1,couternodeupper1,
     &  couternodelower2,couternodeupper2)
c***********************************************************************
      implicit none

      integer
     &  filo0,filo1,filo2,
     &  fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,
     &  cihi0,cihi1,cihi2
      integer ratio(0:NDIM-1),depth
      integer ngc0,ngc1,ngc2
      double precision
     &  fnode(NODE3dVECG(filo,fihi,ngc),depth),
     &  couternodelower0(OUTERNODE3d0(cilo,cihi,0),depth),
     &  couternodeupper0(OUTERNODE3d0(cilo,cihi,0),depth),
     &  couternodelower1(OUTERNODE3d1(cilo,cihi,0),depth),
     &  couternodeupper1(OUTERNODE3d1(cilo,cihi,0),depth),
     &  couternodelower2(OUTERNODE3d2(cilo,cihi,0),depth),
     &  couternodeupper2(OUTERNODE3d2(cilo,cihi,0),depth)
      integer ic0,if0,ic1,if1,ic2,if2,id
      integer cilo0_loop,cilo1_loop,cilo2_loop,
     &        cihi0_loop,cihi1_loop,cihi2_loop
c
c***********************************************************************
c 

c sum along XLOWER face
      cilo1_loop=cilo1+1
      cihi1_loop=cihi1
      cilo2_loop=cilo2+1
      cihi2_loop=cihi2

      if0=cilo0*ratio(0)
do_cfsum_along_patch_face(1,2,`lower0')dnl

c sum along XUPPER face
      cilo1_loop=cilo1+1
      cihi1_loop=cihi1
      cilo2_loop=cilo2+1
      cihi2_loop=cihi2

      if0=(cihi0+1)*ratio(0)
do_cfsum_along_patch_face(1,2,`upper0')dnl

c sum along YLOWER face
      cilo0_loop=cilo0
      cihi0_loop=cihi0+1
      cilo2_loop=cilo2+1
      cihi2_loop=cihi2

      if1=cilo1*ratio(1)
do_cfsum_along_patch_face(0,2,`lower1')dnl

c sum along YUPPER face
      cilo0_loop=cilo0
      cihi0_loop=cihi0+1
      cilo2_loop=cilo2+1
      cihi2_loop=cihi2

      if1=(cihi1+1)*ratio(1)
do_cfsum_along_patch_face(0,2,`upper1')dnl

c sum along ZLOWER face
      cilo0_loop=cilo0
      cihi0_loop=cihi0+1
      cilo1_loop=cilo1
      cihi1_loop=cihi1+1

      if2=cilo2*ratio(2)
do_cfsum_along_patch_face(0,1,`lower2')dnl

c sum along ZUPPER face
      cilo0_loop=cilo0
      cihi0_loop=cihi0+1
      cilo1_loop=cilo1
      cihi1_loop=cihi1+1

      if2=(cihi2+1)*ratio(2)
do_cfsum_along_patch_face(0,1,`upper2')dnl

c
      return
      end
c

c
c***********************************************************************
c Fill hanging nodes on fine patch (i.e. those nodes that do not overlap
c a coarse level node) along patch boundary by interpolation from
c appropriate neighboring coarse nodes.
c***********************************************************************
c
      subroutine nodehangnodeinterp3d(
     &  filo0,filo1,filo2,
     &  fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,
     &  cihi0,cihi1,cihi2,
     &  bboxilo0,bboxilo1,bboxilo2,
     &  bboxihi0,bboxihi1,bboxihi2,
     &  bboxloc,
     &  ratio,
     &  depth,
     &  ngc0,ngc1,ngc2,
     &  fnode)
c***********************************************************************
      implicit none

      integer XLOWER,XUPPER,YLOWER,YUPPER,ZLOWER,ZUPPER
      parameter (XLOWER=0)
      parameter (XUPPER=1)
      parameter (YLOWER=2)
      parameter (YUPPER=3)
      parameter (ZLOWER=4)
      parameter (ZUPPER=5)

      double precision one
      parameter (one=1.0d0)

c***********************************************************************
      integer
     &  filo0,filo1,filo2,
     &  fihi0,fihi1,fihi2,
     &  cilo0,cilo1,cilo2,
     &  cihi0,cihi1,cihi2
      integer
     &  bboxilo0,bboxilo1,bboxilo2,
     &  bboxihi0,bboxihi1,bboxihi2,
     &  bboxloc
      integer ratio(0:NDIM-1),depth
      integer ngc0,ngc1,ngc2
      double precision
     &  fnode(NODE3dVECG(filo,fihi,ngc),depth)
      integer ic0,if0,ir0,ic1,if1,ir1,ic2,if2,ir2,id
      integer iflower0,iflower1,iflower2,ifupper0,ifupper1,ifupper2
      integer cilo0_loop,cilo1_loop,cilo2_loop,
     &        cihi0_loop,cihi1_loop,cihi2_loop
      double precision dratio0,dratio1,dratio2,x0,x1,x2
c
c***********************************************************************
c

      dratio0 = dble(ratio(0))
      dratio1 = dble(ratio(1))
      dratio2 = dble(ratio(2))

      cilo0_loop=max(cilo0, bboxilo0)
      cihi0_loop=min(cihi0, bboxihi0)
      cilo1_loop=max(cilo1, bboxilo1)
      cihi1_loop=min(cihi1, bboxihi1)
      cilo2_loop=max(cilo2, bboxilo2)
      cihi2_loop=min(cihi2, bboxihi2)

      if (bboxloc.eq.XLOWER) then
      
         if0=cilo0*ratio(0)

         do id = 1, depth

c interpolate nodes along lines in y-direction
do_interp_lines1(2)dnl

c interpolate nodes along lines in z-direction
do_interp_lines2(1)dnl

c interpolate face nodes in coarse cell interiors
do_interp_interior0()dnl

         enddo

      else if (bboxloc.eq.XUPPER) then

         if0=(cihi0+1)*ratio(0)

         do id = 1, depth

c interpolate nodes along lines in y-direction
do_interp_lines1(2)dnl

c interpolate nodes along lines in z-direction
do_interp_lines2(1)dnl

c interpolate face nodes in coarse cell interiors
do_interp_interior0()dnl

         enddo

      else if (bboxloc.eq.YLOWER) then

         if1=cilo1*ratio(1)

         do id = 1, depth

c interpolate nodes along lines in x-direction
do_interp_lines0(2)dnl

c interpolate nodes along lines in z-direction
do_interp_lines2(0)dnl

c interpolate face nodes in coarse cell interiors
do_interp_interior1()dnl

         enddo

      else if (bboxloc.eq.YUPPER) then

         if1=(cihi1+1)*ratio(1)
   
         do id = 1, depth

c interpolate nodes along lines in x-direction
do_interp_lines0(2)dnl

c interpolate nodes along lines in z-direction
do_interp_lines2(0)dnl

c interpolate face nodes in coarse cell interiors
do_interp_interior1()dnl

         enddo

      else if (bboxloc.eq.ZLOWER) then

         if2=cilo2*ratio(2)

         do id = 1, depth

c interpolate nodes along lines in x-direction
do_interp_lines0(1)dnl

c interpolate nodes along lines in y-direction
do_interp_lines1(0)dnl

c interpolate face nodes in coarse cell interiors
do_interp_interior2()dnl

         enddo

      else if (bboxloc.eq.ZUPPER) then

         if2=(cihi2+1)*ratio(2)

         do id = 1, depth

c interpolate nodes along lines in x-direction
do_interp_lines0(1)dnl

c interpolate nodes along lines in y-direction
do_interp_lines1(0)dnl

c interpolate face nodes in coarse cell interiors
do_interp_interior2()dnl

         enddo

      endif
c
      return
      end
c


c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for summing outernode data with other node or
c                outernode data.
c
define(NDIM,2)dnl
include(PDAT_FORTDIR/pdat_m4arrdim2d.i)dnl

define(do_cfsum_along_patch_side,`dnl
      do id = 1, depth
         do ic$1=cilo$1_loop,cihi$1_loop
            if$1=ic$1*ratio($1)
            fnode(if0,if1,id) = fnode(if0,if1,id) + 
     &                          couternode$2(ic$1,id)
         enddo
      enddo
')dnl

define(do_hanging_node_interpolation,`dnl
      do id = 1, depth
         do ic$1=cilo$1_loop,cihi$1_loop
            if$1=ic$1*ratio($1)
            do ir$1=1,ratio($1)-1
               x$1 = dble(ir$1)/dratio$1
               fnode($2,id) =
     &            fnode(if0,if1,id) * (one - x$1) +
     &            fnode($3,id) * x$1
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
      subroutine nodeouternodesum2d(
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  ratio,
     &  depth,
     &  ngc0,ngc1,
     &  fnode,
     &  couternodelower0,couternodeupper0,
     &  couternodelower1,couternodeupper1)
c***********************************************************************
      implicit none

      integer
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1
      integer ratio(0:NDIM-1), depth
      integer ngc0,ngc1
      double precision
     &  fnode(NODE2dVECG(filo,fihi,ngc),depth),
     &  couternodelower0(OUTERNODE2d0(cilo,cihi,0),depth),
     &  couternodeupper0(OUTERNODE2d0(cilo,cihi,0),depth),
     &  couternodelower1(OUTERNODE2d1(cilo,cihi,0),depth),
     &  couternodeupper1(OUTERNODE2d1(cilo,cihi,0),depth)
      integer ic0,if0,ic1,if1,id
      integer cilo0_loop,cilo1_loop,cihi0_loop,cihi1_loop
c
c***********************************************************************
c

c sum along XLOWER side
      cilo1_loop=cilo1+1
      cihi1_loop=cihi1

      if0=cilo0*ratio(0)
do_cfsum_along_patch_side(1,`lower0')dnl

c sum along XUPPER side
      cilo1_loop=cilo1+1
      cihi1_loop=cihi1

      if0=(cihi0+1)*ratio(0)
do_cfsum_along_patch_side(1,`upper0')dnl

c sum along YLOWER side
      cilo0_loop=cilo0
      cihi0_loop=cihi0+1

      if1=cilo1*ratio(1)
do_cfsum_along_patch_side(0,`lower1')dnl

c sum along YUPPER side
      cilo0_loop=cilo0
      cihi0_loop=cihi0+1

      if1=(cihi1+1)*ratio(1)
do_cfsum_along_patch_side(0,`upper1')dnl

c
      return
      end

c
c***********************************************************************
c Fill hanging nodes on fine patch (i.e. those nodes that do not overlap
c a coarse level node) along patch boundary by interpolation from 
c appropriate neighboring coarse nodes
c***********************************************************************
c
      subroutine nodehangnodeinterp2d(
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1,
     &  bboxilo0,bboxilo1,bboxihi0,bboxihi1,
     &  bboxloc,
     &  ratio,
     &  depth,
     &  ngc0,ngc1,
     &  fnode)
c***********************************************************************
      implicit none

      integer XLOWER,XUPPER,YLOWER,YUPPER
      parameter (XLOWER=0)
      parameter (XUPPER=1)
      parameter (YLOWER=2)
      parameter (YUPPER=3)

      double precision one
      parameter (one=1.0d0)
c***********************************************************************
c
      integer
     &  filo0,filo1,fihi0,fihi1,
     &  cilo0,cilo1,cihi0,cihi1
      integer bboxilo0,bboxilo1,bboxihi0,bboxihi1,bboxloc
      integer ratio(0:NDIM-1),depth
      integer ngc0,ngc1
      double precision
     &  fnode(NODE2dVECG(filo,fihi,ngc),depth)
      integer ic0,if0,ir0,ic1,if1,ir1,id
      integer cilo0_loop,cilo1_loop,cihi0_loop,cihi1_loop
      double precision dratio0,dratio1,x0,x1
c
c***********************************************************************
c
      dratio0 = dble(ratio(0))
      dratio1 = dble(ratio(1))

      cilo0_loop=max(cilo0, bboxilo0)
      cihi0_loop=min(cihi0, bboxihi0)
      cilo1_loop=max(cilo1, bboxilo1)
      cihi1_loop=min(cihi1, bboxihi1)
c
c***********************************************************************
c 
      if (bboxloc.eq.XLOWER) then

         if0=cilo0*ratio(0)
do_hanging_node_interpolation(1,`if0,if1+ir1',`if0,if1+ratio(1)')dnl

      else if (bboxloc.eq.XUPPER) then
 
         if0=(cihi0+1)*ratio(0)
do_hanging_node_interpolation(1,`if0,if1+ir1',`if0,if1+ratio(1)')dnl

      else if (bboxloc.eq.YLOWER) then
      
         if1=cilo1*ratio(1)
do_hanging_node_interpolation(0,`if0+ir0,if1',`if0+ratio(0),if1')dnl

      else if (bboxloc.eq.YUPPER) then

         if1=(cihi1+1)*ratio(1)
do_hanging_node_interpolation(0,`if0+ir0,if1',`if0+ratio(0),if1')dnl
  
      endif
c
      return
      end

c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   FORTRAN routines for spatial coarsening of 2d patch data
c                on a regular Cartesian mesh.
c
include(PDAT_FORTDIR/pdat_m4concoarsenops2d.i)dnl
c
c***********************************************************************
c Constant coarsening for 2d node-centered double data
c***********************************************************************
c
      subroutine conavgnodedoub2d(
conavg_op_node_2d(`double precision')dnl
c
c***********************************************************************
c Constant coarsening for 2d node-centered float data 
c***********************************************************************
c
      subroutine conavgnodeflot2d(
conavg_op_node_2d(`real')dnl
c
c***********************************************************************
c Constant coarsening for 2d node-centered complex data
c***********************************************************************
c
      subroutine conavgnodecplx2d(
conavg_op_node_2d(`double complex')dnl
c
c***********************************************************************
c Constant coarsening for 2d node-centered integer data
c***********************************************************************
c
      subroutine conavgnodeintg2d(
conavg_op_node_2d(`integer')dnl
c
c

c***********************************************************************
c Constant coarsening for 2d outernode-centered double data
c***********************************************************************
c
      subroutine conavgouternodedoub2d0(
conavg_op_outernode_2d0(`double precision')dnl

      subroutine conavgouternodedoub2d1(
conavg_op_outernode_2d1(`double precision')dnl

c***********************************************************************
c Constant coarsening for 2d outernode-centered float data
c***********************************************************************
c
      subroutine conavgouternodeflot2d0(
conavg_op_outernode_2d0(`real')dnl

      subroutine conavgouternodeflot2d1(
conavg_op_outernode_2d1(`real')dnl

c
c***********************************************************************
c Constant coarsening for 2d outernode-centered complex data
c***********************************************************************
c
      subroutine conavgouternodecplx2d0(
conavg_op_outernode_2d0(`complex')dnl

      subroutine conavgouternodecplx2d1(
conavg_op_outernode_2d1(`complex')dnl

c
c***********************************************************************
c Constant coarsening for 2d outernode-centered integer data
c***********************************************************************
c
      subroutine conavgouternodeint2d0(
conavg_op_outernode_2d0(`integer')dnl

      subroutine conavgouternodeint2d1(
conavg_op_outernode_2d1(`integer')dnl

c


c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   FORTRAN routines for standard linear time interpolation 
c                of 3d patch data types.
c
include(PDAT_FORTDIR/pdat_m4lintimeintops3d.i)dnl
c
c***********************************************************************
c Linear time interpolation for 3d cell-centered double data
c***********************************************************************
c
      subroutine lintimeintcelldoub3d(
lin_time_int_op_cell_3d(`double precision')dnl
c
c***********************************************************************
c Linear time interpolation for 3d cell-centered float data
c***********************************************************************
c
      subroutine lintimeintcellfloat3d(
lin_time_int_op_cell_3d(`real')dnl
c
c***********************************************************************
c Linear time interpolation for 3d cell-centered complex data
c***********************************************************************
c
      subroutine lintimeintcellcmplx3d(
lin_time_int_op_cell_3d(`double complex')dnl
c
c***********************************************************************
c Linear time interpolation for 3d edge-centered double data
c***********************************************************************
c
      subroutine lintimeintedgedoub3d0(
lin_time_int_op_edge_3d(`double precision',0,1,2)dnl
c
      subroutine lintimeintedgedoub3d1(
lin_time_int_op_edge_3d(`double precision',1,2,0)dnl
c
      subroutine lintimeintedgedoub3d2(
lin_time_int_op_edge_3d(`double precision',2,0,1)dnl
c
c***********************************************************************
c Linear time interpolation for 3d edge-centered float data
c***********************************************************************
c
      subroutine lintimeintedgefloat3d0(
lin_time_int_op_edge_3d(`real',0,1,2)dnl
c
      subroutine lintimeintedgefloat3d1(
lin_time_int_op_edge_3d(`real',1,2,0)dnl
c
      subroutine lintimeintedgefloat3d2(
lin_time_int_op_edge_3d(`real',2,0,1)dnl
c
c***********************************************************************
c Linear time interpolation for 3d edge-centered complex data
c***********************************************************************
c
      subroutine lintimeintedgecmplx3d0(
lin_time_int_op_edge_3d(`double complex',0,1,2)dnl
c
      subroutine lintimeintedgecmplx3d1(
lin_time_int_op_edge_3d(`double complex',1,2,0)dnl
c
      subroutine lintimeintedgecmplx3d2(
lin_time_int_op_edge_3d(`double complex',2,0,1)dnl
c
c***********************************************************************
c Linear time interpolation for 3d face-centered double data
c***********************************************************************
c
      subroutine lintimeintfacedoub3d0(
lin_time_int_op_face_3d(`double precision',0,1,2)dnl
c
      subroutine lintimeintfacedoub3d1(
lin_time_int_op_face_3d(`double precision',1,2,0)dnl
c
      subroutine lintimeintfacedoub3d2(
lin_time_int_op_face_3d(`double precision',2,0,1)dnl
c
c***********************************************************************
c Linear time interpolation for 3d face-centered float data
c***********************************************************************
c
      subroutine lintimeintfacefloat3d0(
lin_time_int_op_face_3d(`real',0,1,2)dnl
c
      subroutine lintimeintfacefloat3d1(
lin_time_int_op_face_3d(`real',1,2,0)dnl
c
      subroutine lintimeintfacefloat3d2(
lin_time_int_op_face_3d(`real',2,0,1)dnl
c
c***********************************************************************
c Linear time interpolation for 3d face-centered complex data
c***********************************************************************
c
      subroutine lintimeintfacecmplx3d0(
lin_time_int_op_face_3d(`double complex',0,1,2)dnl
c
      subroutine lintimeintfacecmplx3d1(
lin_time_int_op_face_3d(`double complex',1,2,0)dnl
c
      subroutine lintimeintfacecmplx3d2(
lin_time_int_op_face_3d(`double complex',2,0,1)dnl
c
c***********************************************************************
c Linear time interpolation for 3d node-centered double data
c***********************************************************************
c
      subroutine lintimeintnodedoub3d(
lin_time_int_op_node_3d(`double precision')dnl
c
c***********************************************************************
c Linear time interpolation for 3d node-centered float data
c***********************************************************************
c
      subroutine lintimeintnodefloat3d(
lin_time_int_op_node_3d(`real')dnl
c
c***********************************************************************
c Linear time interpolation for 3d node-centered complex data
c***********************************************************************
c
      subroutine lintimeintnodecmplx3d(
lin_time_int_op_node_3d(`double complex')dnl
c
c***********************************************************************
c Linear time interpolation for 3d outerface double data
c***********************************************************************
c
      subroutine lintimeintoutfacedoub3d0(
lin_time_int_op_outerface_3d(`double precision',0,1,2)dnl
c
      subroutine lintimeintoutfacedoub3d1(
lin_time_int_op_outerface_3d(`double precision',1,2,0)dnl
c
      subroutine lintimeintoutfacedoub3d2(
lin_time_int_op_outerface_3d(`double precision',2,0,1)dnl
c
c***********************************************************************
c Linear time interpolation for 3d outerface float data
c***********************************************************************
c
      subroutine lintimeintoutfacefloat3d0(
lin_time_int_op_outerface_3d(`real',0,1,2)dnl
c
      subroutine lintimeintoutfacefloat3d1(
lin_time_int_op_outerface_3d(`real',1,2,0)dnl
c
      subroutine lintimeintoutfacefloat3d2(
lin_time_int_op_outerface_3d(`real',2,0,1)dnl
c
c***********************************************************************
c Linear time interpolation for 3d outerface complex data
c***********************************************************************
c
      subroutine lintimeintoutfacecmplx3d0(
lin_time_int_op_outerface_3d(`double complex',0,1,2)dnl
c
      subroutine lintimeintoutfacecmplx3d1(
lin_time_int_op_outerface_3d(`double complex',1,2,0)dnl
c
      subroutine lintimeintoutfacecmplx3d2(
lin_time_int_op_outerface_3d(`double complex',2,0,1)dnl
c
c***********************************************************************
c Linear time interpolation for 3d outerside double data
c***********************************************************************
c
      subroutine lintimeintoutsidedoub3d0(
lin_time_int_op_outerside_3d(`double precision',0,1,2)dnl
c
      subroutine lintimeintoutsidedoub3d1(
lin_time_int_op_outerside_3d(`double precision',1,0,2)dnl
c
      subroutine lintimeintoutsidedoub3d2(
lin_time_int_op_outerside_3d(`double precision',2,0,1)dnl
c
c***********************************************************************
c Linear time interpolation for 3d outerside float data
c***********************************************************************
c
      subroutine lintimeintoutsidefloat3d0(
lin_time_int_op_outerside_3d(`real',0,1,2)dnl
c
      subroutine lintimeintoutsidefloat3d1(
lin_time_int_op_outerside_3d(`real',1,0,2)dnl
c
      subroutine lintimeintoutsidefloat3d2(
lin_time_int_op_outerside_3d(`real',2,0,1)dnl
c
c***********************************************************************
c Linear time interpolation for 3d outerside complex data
c***********************************************************************
c
      subroutine lintimeintoutsidecmplx3d0(
lin_time_int_op_outerside_3d(`double complex',0,1,2)dnl
c
      subroutine lintimeintoutsidecmplx3d1(
lin_time_int_op_outerside_3d(`double complex',1,0,2)dnl
c
      subroutine lintimeintoutsidecmplx3d2(
lin_time_int_op_outerside_3d(`double complex',2,0,1)dnl
c
c***********************************************************************
c Linear time interpolation for 3d side-centered double data
c***********************************************************************
c
      subroutine lintimeintsidedoub3d0(
lin_time_int_op_side_3d(`double precision',0,1,2)dnl
c
      subroutine lintimeintsidedoub3d1(
lin_time_int_op_side_3d(`double precision',1,2,0)dnl
c
      subroutine lintimeintsidedoub3d2(
lin_time_int_op_side_3d(`double precision',2,0,1)dnl
c
c***********************************************************************
c Linear time interpolation for 3d side-centered float data
c***********************************************************************
c
      subroutine lintimeintsidefloat3d0(
lin_time_int_op_side_3d(`real',0,1,2)dnl
c
      subroutine lintimeintsidefloat3d1(
lin_time_int_op_side_3d(`real',1,2,0)dnl
c
      subroutine lintimeintsidefloat3d2(
lin_time_int_op_side_3d(`real',2,0,1)dnl
c
c***********************************************************************
c Linear time interpolation for 3d side-centered complex data
c***********************************************************************
c
      subroutine lintimeintsidecmplx3d0(
lin_time_int_op_side_3d(`double complex',0,1,2)dnl
c
      subroutine lintimeintsidecmplx3d1(
lin_time_int_op_side_3d(`double complex',1,2,0)dnl
c
      subroutine lintimeintsidecmplx3d2(
lin_time_int_op_side_3d(`double complex',2,0,1)dnl
c

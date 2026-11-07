c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   FORTRAN routines for standard linear time interpolation 
c                of 2d patch data types.
c
include(PDAT_FORTDIR/pdat_m4lintimeintops2d.i)dnl
c
c***********************************************************************
c Linear time interpolation for 2d cell-centered double data
c***********************************************************************
c
      subroutine lintimeintcelldoub2d(
lin_time_int_op_cell_2d(`double precision')dnl
c
c***********************************************************************
c Linear time interpolation for 2d cell-centered float data
c***********************************************************************
c
      subroutine lintimeintcellfloat2d(
lin_time_int_op_cell_2d(`real')dnl
c
c***********************************************************************
c Linear time interpolation for 2d cell-centered complex data
c***********************************************************************
c
      subroutine lintimeintcellcmplx2d(
lin_time_int_op_cell_2d(`double complex')dnl
c
c***********************************************************************
c Linear time interpolation for 2d edge-centered double data
c***********************************************************************
c
      subroutine lintimeintedgedoub2d0(
lin_time_int_op_edge_2d(`double precision',0,1)dnl
c
      subroutine lintimeintedgedoub2d1(
lin_time_int_op_edge_2d(`double precision',1,0)dnl
c
c***********************************************************************
c Linear time interpolation for 2d edge-centered float data
c***********************************************************************
c
      subroutine lintimeintedgefloat2d0(
lin_time_int_op_edge_2d(`real',0,1)dnl
c
      subroutine lintimeintedgefloat2d1(
lin_time_int_op_edge_2d(`real',1,0)dnl
c
c***********************************************************************
c Linear time interpolation for 2d edge-centered complex data
c***********************************************************************
c
      subroutine lintimeintedgecmplx2d0(
lin_time_int_op_edge_2d(`double complex',0,1)dnl
c
      subroutine lintimeintedgecmplx2d1(
lin_time_int_op_edge_2d(`double complex',1,0)dnl
c
c***********************************************************************
c Linear time interpolation for 2d face-centered double data
c***********************************************************************
c
      subroutine lintimeintfacedoub2d0(
lin_time_int_op_face_2d(`double precision',0,1)dnl
c
      subroutine lintimeintfacedoub2d1(
lin_time_int_op_face_2d(`double precision',1,0)dnl
c
c***********************************************************************
c Linear time interpolation for 2d face-centered float data
c***********************************************************************
c
      subroutine lintimeintfacefloat2d0(
lin_time_int_op_face_2d(`real',0,1)dnl
c
      subroutine lintimeintfacefloat2d1(
lin_time_int_op_face_2d(`real',1,0)dnl
c
c***********************************************************************
c Linear time interpolation for 2d face-centered complex data
c***********************************************************************
c
      subroutine lintimeintfacecmplx2d0(
lin_time_int_op_face_2d(`double complex',0,1)dnl
c
      subroutine lintimeintfacecmplx2d1(
lin_time_int_op_face_2d(`double complex',1,0)dnl
c
c***********************************************************************
c Linear time interpolation for 2d node-centered double data
c***********************************************************************
c
      subroutine lintimeintnodedoub2d(
lin_time_int_op_node_2d(`double precision')dnl
c
c***********************************************************************
c Linear time interpolation for 2d node-centered float data
c***********************************************************************
c
      subroutine lintimeintnodefloat2d(
lin_time_int_op_node_2d(`real')dnl
c
c***********************************************************************
c Linear time interpolation for 2d node-centered complex data
c***********************************************************************
c
      subroutine lintimeintnodecmplx2d(
lin_time_int_op_node_2d(`double complex')dnl
c
c***********************************************************************
c Linear time interpolation for 2d outerface double data
c***********************************************************************
c
      subroutine lintimeintoutfacedoub2d0(
lin_time_int_op_outerface_2d(`double precision',0,1)dnl
c
      subroutine lintimeintoutfacedoub2d1(
lin_time_int_op_outerface_2d(`double precision',1,0)dnl
c
c***********************************************************************
c Linear time interpolation for 2d outerface float data
c***********************************************************************
c
      subroutine lintimeintoutfacefloat2d0(
lin_time_int_op_outerface_2d(`real',0,1)dnl
c
      subroutine lintimeintoutfacefloat2d1(
lin_time_int_op_outerface_2d(`real',1,0)dnl
c
c***********************************************************************
c Linear time interpolation for 2d outerface complex data
c***********************************************************************
c
      subroutine lintimeintoutfacecmplx2d0(
lin_time_int_op_outerface_2d(`double complex',0,1)dnl
c
      subroutine lintimeintoutfacecmplx2d1(
lin_time_int_op_outerface_2d(`double complex',1,0)dnl
c
c***********************************************************************
c Linear time interpolation for 2d outerside double data
c***********************************************************************
c
      subroutine lintimeintoutsidedoub2d0(
lin_time_int_op_outerside_2d(`double precision',0,1)dnl
c
      subroutine lintimeintoutsidedoub2d1(
lin_time_int_op_outerside_2d(`double precision',1,0)dnl
c
c***********************************************************************
c Linear time interpolation for 2d outerside float data
c***********************************************************************
c
      subroutine lintimeintoutsidefloat2d0(
lin_time_int_op_outerside_2d(`real',0,1)dnl
c
      subroutine lintimeintoutsidefloat2d1(
lin_time_int_op_outerside_2d(`real',1,0)dnl
c
c***********************************************************************
c Linear time interpolation for 2d outerside complex data
c***********************************************************************
c
      subroutine lintimeintoutsidecmplx2d0(
lin_time_int_op_outerside_2d(`double complex',0,1)dnl
c
      subroutine lintimeintoutsidecmplx2d1(
lin_time_int_op_outerside_2d(`double complex',1,0)dnl
c
c***********************************************************************
c Linear time interpolation for 2d side-centered double data
c***********************************************************************
c
      subroutine lintimeintsidedoub2d0(
lin_time_int_op_side_2d(`double precision',0,1)dnl
c
      subroutine lintimeintsidedoub2d1(
lin_time_int_op_side_2d(`double precision',1,0)dnl
c
c***********************************************************************
c Linear time interpolation for 2d side-centered float data
c***********************************************************************
c
      subroutine lintimeintsidefloat2d0(
lin_time_int_op_side_2d(`real',0,1)dnl
c
      subroutine lintimeintsidefloat2d1(
lin_time_int_op_side_2d(`real',1,0)dnl
c
c***********************************************************************
c Linear time interpolation for 2d side-centered complex data
c***********************************************************************
c
      subroutine lintimeintsidecmplx2d0(
lin_time_int_op_side_2d(`double complex',0,1)dnl
c
      subroutine lintimeintsidecmplx2d1(
lin_time_int_op_side_2d(`double complex',1,0)dnl
c

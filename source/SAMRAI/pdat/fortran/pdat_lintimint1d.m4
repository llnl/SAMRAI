c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   FORTRAN routines for standard linear time interpolation
c                of 1d patch data types.
c
include(PDAT_FORTDIR/pdat_m4lintimeintops1d.i)dnl
c
c***********************************************************************
c Linear time interpolation for 1d cell-centered double data
c***********************************************************************
c
      subroutine lintimeintcelldoub1d(
lin_time_int_op_cell_1d(`double precision')dnl
c
c***********************************************************************
c Linear time interpolation for 1d cell-centered float data
c***********************************************************************
c
      subroutine lintimeintcellfloat1d(
lin_time_int_op_cell_1d(`real')dnl
c
c***********************************************************************
c Linear time interpolation for 1d cell-centered complex data
c***********************************************************************
c
      subroutine lintimeintcellcmplx1d(
lin_time_int_op_cell_1d(`double complex')dnl
c
c***********************************************************************
c Linear time interpolation for 1d edge-centered double data
c***********************************************************************
c
      subroutine lintimeintedgedoub1d(
lin_time_int_op_edge_1d(`double precision')dnl
c
c***********************************************************************
c Linear time interpolation for 1d edge-centered float data
c***********************************************************************
c
      subroutine lintimeintedgefloat1d(
lin_time_int_op_edge_1d(`real',)dnl
c
c***********************************************************************
c Linear time interpolation for 1d edge-centered complex data
c***********************************************************************
c
      subroutine lintimeintedgecmplx1d(
lin_time_int_op_edge_1d(`double complex')dnl
c
c***********************************************************************
c Linear time interpolation for 1d face-centered double data
c***********************************************************************
c
      subroutine lintimeintfacedoub1d(
lin_time_int_op_face_1d(`double precision')dnl
c
c***********************************************************************
c Linear time interpolation for 1d face-centered float data
c***********************************************************************
c
      subroutine lintimeintfacefloat1d(
lin_time_int_op_face_1d(`real',)dnl
c
c***********************************************************************
c Linear time interpolation for 1d face-centered complex data
c***********************************************************************
c
      subroutine lintimeintfacecmplx1d(
lin_time_int_op_face_1d(`double complex')dnl
c
c***********************************************************************
c Linear time interpolation for 1d node-centered double data
c***********************************************************************
c
      subroutine lintimeintnodedoub1d(
lin_time_int_op_node_1d(`double precision')dnl
c
c***********************************************************************
c Linear time interpolation for 1d node-centered float data
c***********************************************************************
c
      subroutine lintimeintnodefloat1d(
lin_time_int_op_node_1d(`real')dnl
c
c***********************************************************************
c Linear time interpolation for 1d node-centered complex data
c***********************************************************************
c
      subroutine lintimeintnodecmplx1d(
lin_time_int_op_node_1d(`double complex')dnl
c
c***********************************************************************
c Linear time interpolation for 1d outerface double data
c***********************************************************************
c
      subroutine lintimeintoutfacedoub1d(
lin_time_int_op_outerface_1d(`double precision')dnl
c
c***********************************************************************
c Linear time interpolation for 1d outerface float data
c***********************************************************************
c
      subroutine lintimeintoutfacefloat1d(
lin_time_int_op_outerface_1d(`real')dnl
c
c***********************************************************************
c Linear time interpolation for 1d outerface complex data
c***********************************************************************
c
      subroutine lintimeintoutfacecmplx1d(
lin_time_int_op_outerface_1d(`double complex')dnl
c
c***********************************************************************
c Linear time interpolation for 1d outerside double data
c***********************************************************************
c
      subroutine lintimeintoutsidedoub1d(
lin_time_int_op_outerside_1d(`double precision')dnl
c
c***********************************************************************
c Linear time interpolation for 1d outerside float data
c***********************************************************************
c
      subroutine lintimeintoutsidefloat1d(
lin_time_int_op_outerside_1d(`real')dnl
c
c***********************************************************************
c Linear time interpolation for 1d outerside complex data
c***********************************************************************
c
      subroutine lintimeintoutsidecmplx1d(
lin_time_int_op_outerside_1d(`double complex')dnl
c
c***********************************************************************
c Linear time interpolation for 1d side-centered double data
c***********************************************************************
c
      subroutine lintimeintsidedoub1d(
lin_time_int_op_side_1d(`double precision')dnl
c
c***********************************************************************
c Linear time interpolation for 1d side-centered float data
c***********************************************************************
c
      subroutine lintimeintsidefloat1d(
lin_time_int_op_side_1d(`real',)dnl
c
c***********************************************************************
c Linear time interpolation for 1d side-centered complex data
c***********************************************************************
c
      subroutine lintimeintsidecmplx1d(
lin_time_int_op_side_1d(`double complex')dnl
c

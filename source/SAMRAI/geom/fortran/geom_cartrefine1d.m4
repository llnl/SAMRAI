c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   FORTRAN routines for spatial refining of 1d patch data
c                on a regular Cartesian mesh.
c
include(FORTDIR/geom_m4cartrefineops1d.i)dnl
c
c***********************************************************************
c Linear interpolation for 1d cell-centered double data
c***********************************************************************
c
      subroutine cartlinrefcelldoub1d(
cart_linref_op_cell_1d(`double precision')dnl
c
c***********************************************************************
c Linear interpolation for 1d cell-centered float data
c***********************************************************************
c
      subroutine cartlinrefcellflot1d(
cart_linref_op_cell_1d(`real')dnl
c
c***********************************************************************
c Linear interpolation for 1d cell-centered complex data
c***********************************************************************
c
      subroutine cartlinrefcellcplx1d(
cart_linref_op_cell_1d(`double complex')dnl
c
c***********************************************************************
c Conservative linear interpolation for 1d cell-centered double data
c***********************************************************************
c
      subroutine cartclinrefcelldoub1d(
cart_clinref_op_cell_1d(`double precision')dnl
c
c***********************************************************************
c Conservative linear interpolation for 1d cell-centered float data
c***********************************************************************
c
      subroutine cartclinrefcellflot1d(
cart_clinref_op_cell_1d(`real')dnl
c
c***********************************************************************
c Conservative linear interpolation for 1d cell-centered complex data
c***********************************************************************
c
      subroutine cartclinrefcellcplx1d(
cart_clinref_op_cell_1d(`double complex')dnl
c
c***********************************************************************
c Conservative linear interpolation for 1d edge-centered double data
c***********************************************************************
c
      subroutine cartclinrefedgedoub1d(
cart_clinref_op_edge_1d(`double precision')dnl
c
c***********************************************************************
c Conservative linear interpolation for 1d edge-centered float data
c***********************************************************************
c
      subroutine cartclinrefedgeflot1d(
cart_clinref_op_edge_1d(`real')dnl
c
c***********************************************************************
c Conservative linear interpolation for 1d face-centered double data
c***********************************************************************
c
      subroutine cartclinreffacedoub1d(
cart_clinref_op_face_1d(`double precision')dnl
c
c***********************************************************************
c Conservative linear interpolation for 1d face-centered float data
c***********************************************************************
c
      subroutine cartclinreffaceflot1d(
cart_clinref_op_face_1d(`real')dnl
c
c***********************************************************************
c Conservative linear interpolation for 1d face-centered complex data
c***********************************************************************
c
c      subroutine cartclinreffacecplx1d(
ccart_clinref_op_face_1d(`double complex')dnl
c
c***********************************************************************
c Linear interpolation for 1d node-centered double data
c***********************************************************************
c
       subroutine cartlinrefnodedoub1d(
cart_linref_op_node_1d(`double precision')dnl
c
c***********************************************************************
c Linear interpolation for 1d node-centered float data
c***********************************************************************
c
       subroutine cartlinrefnodeflot1d(
cart_linref_op_node_1d(`real')dnl
c
c***********************************************************************
c Linear interpolation for 1d node-centered complex data
c***********************************************************************
c
       subroutine cartlinrefnodecplx1d(
cart_linref_op_node_1d(`double complex')dnl
c
c***********************************************************************
c Conservative linear interpolation for 1d outerface double data
c***********************************************************************
c
      subroutine cartclinrefoutfacedoub1d(
cart_clinref_op_outerface_1d(`double precision')dnl
c
c***********************************************************************
c Conservative linear interpolation for 1d outerface float data
c***********************************************************************
c
      subroutine cartclinrefoutfaceflot1d(
cart_clinref_op_face_1d(`real')dnl
c
c***********************************************************************
c Conservative linear interpolation for 1d outerface complex data
c***********************************************************************
c
c      subroutine cartclinrefoutfacecplx1d(
ccart_clinref_op_face_1d(`double complex')dnl
c
c***********************************************************************
c Conservative linear interpolation for 1d side-centered double data
c***********************************************************************
c
      subroutine cartclinrefsidedoub1d(
cart_clinref_op_side_1d(`double precision')dnl
c
c***********************************************************************
c Conservative linear interpolation for 1d side-centered float data
c***********************************************************************
c
      subroutine cartclinrefsideflot1d(
cart_clinref_op_side_1d(`real')dnl
c
c***********************************************************************
c Conservative linear interpolation for 1d side-centered complex data
c***********************************************************************
c
c      subroutine cartclinrefsidecplx1d(
ccart_clinref_op_side_1d(`double complex')dnl
c

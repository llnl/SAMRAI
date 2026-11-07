c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   FORTRAN routines for spatial refining of 3d patch data
c                on a regular Cartesian mesh.
c
include(FORTDIR/geom_m4cartrefineops3d.i)dnl
c
c***********************************************************************
c Linear interpolation for 3d cell-centered double data
c***********************************************************************
c
      subroutine cartlinrefcelldoub3d(
cart_linref_op_cell_3d(`double precision')dnl
c
c***********************************************************************
c Linear interpolation for 3d cell-centered float data
c***********************************************************************
c
      subroutine cartlinrefcellflot3d(
cart_linref_op_cell_3d(`real')dnl
c
c***********************************************************************
c Linear interpolation for 3d cell-centered complex data
c***********************************************************************
c
      subroutine cartlinrefcellcplx3d(
cart_linref_op_cell_3d(`double complex')dnl
c
c***********************************************************************
c Conservative linear interpolation for 3d cell-centered double data
c***********************************************************************
c
      subroutine cartclinrefcelldoub3d(
cart_clinref_op_cell_3d(`double precision')dnl
c
c***********************************************************************
c Conservative linear interpolation for 3d cell-centered float data
c***********************************************************************
c
      subroutine cartclinrefcellflot3d(
cart_clinref_op_cell_3d(`real')dnl
c
c***********************************************************************
c Conservative linear interpolation for 3d cell-centered complex data
c***********************************************************************
c
      subroutine cartclinrefcellcplx3d(
cart_clinref_op_cell_3d(`double complex')dnl
c
c***********************************************************************
c Conservative linear interpolation for 3d edge-centered double data
c***********************************************************************
c
      subroutine cartclinrefedgedoub3d0(
cart_clinref_op_edge_3d(`double precision',0,1,2)dnl
c
      subroutine cartclinrefedgedoub3d1(
cart_clinref_op_edge_3d(`double precision',1,2,0)dnl
c
      subroutine cartclinrefedgedoub3d2(
cart_clinref_op_edge_3d(`double precision',2,0,1)dnl
c
c***********************************************************************
c Conservative linear interpolation for 3d edge-centered float data
c***********************************************************************
c
      subroutine cartclinrefedgeflot3d0(
cart_clinref_op_edge_3d(`real',0,1,2)dnl
c
      subroutine cartclinrefedgeflot3d1(
cart_clinref_op_edge_3d(`real',1,2,0)dnl
c
      subroutine cartclinrefedgeflot3d2(
cart_clinref_op_edge_3d(`real',2,0,1)dnl
c
c***********************************************************************
c Conservative linear interpolation for 3d face-centered double data
c***********************************************************************
c
      subroutine cartclinreffacedoub3d0(
cart_clinref_op_face_3d(`double precision',0,1,2)dnl
c
      subroutine cartclinreffacedoub3d1(
cart_clinref_op_face_3d(`double precision',1,2,0)dnl
c
      subroutine cartclinreffacedoub3d2(
cart_clinref_op_face_3d(`double precision',2,0,1)dnl
c
c***********************************************************************
c Conservative linear interpolation for 3d face-centered float data
c***********************************************************************
c
      subroutine cartclinreffaceflot3d0(
cart_clinref_op_face_3d(`real',0,1,2)dnl
c
      subroutine cartclinreffaceflot3d1(
cart_clinref_op_face_3d(`real',1,2,0)dnl
c
      subroutine cartclinreffaceflot3d2(
cart_clinref_op_face_3d(`real',2,0,1)dnl
c
c***********************************************************************
c Conservative linear interpolation for 3d face-centered complex data
c***********************************************************************
c
c     subroutine cartclinreffacecplx3d0(
ccart_clinref_op_face_3d(`double complex',0,1,2)dnl
c
c      subroutine cartclinreffacecplx3d1(
ccart_clinref_op_face_3d(`double complex',1,2,0)dnl
c
c      subroutine cartclinreffacecplx3d2(
ccart_clinref_op_face_3d(`double complex',2,0,1)dnl
c
c***********************************************************************
c Linear interpolation for 3d node-centered double data
c***********************************************************************
c
       subroutine cartlinrefnodedoub3d(
cart_linref_op_node_3d(`double precision')dnl
c
c***********************************************************************
c Linear interpolation for 3d node-centered float data
c***********************************************************************
c
       subroutine cartlinrefnodeflot3d(
cart_linref_op_node_3d(`real')dnl
c
c***********************************************************************
c Linear interpolation for 3d node-centered complex data
c***********************************************************************
c
       subroutine cartlinrefnodecplx3d(
cart_linref_op_node_3d(`double complex')dnl
c
c***********************************************************************
c Conservative linear interpolation for 3d outerface double data
c***********************************************************************
c
      subroutine cartclinrefoutfacedoub3d0(
cart_clinref_op_outerface_3d(`double precision',0,1,2)dnl
c
      subroutine cartclinrefoutfacedoub3d1(
cart_clinref_op_outerface_3d(`double precision',1,2,0)dnl
c
      subroutine cartclinrefoutfacedoub3d2(
cart_clinref_op_outerface_3d(`double precision',2,0,1)dnl
c
c***********************************************************************
c Conservative linear interpolation for 3d outerface float data
c***********************************************************************
c
      subroutine cartclinrefoutfaceflot3d0(
cart_clinref_op_outerface_3d(`real',0,1,2)dnl
c
      subroutine cartclinrefoutfaceflot3d1(
cart_clinref_op_outerface_3d(`real',1,2,0)dnl
c
      subroutine cartclinrefoutfaceflot3d2(
cart_clinref_op_outerface_3d(`real',2,0,1)dnl
c
c***********************************************************************
c Conservative linear interpolation for 3d outerface complex data
c***********************************************************************
c
c     subroutine cartclinrefoutfacecplx3d0(
ccart_clinref_op_outerface_3d(`double complex',0,1,2)dnl
c
c      subroutine cartclinrefoutfacecplx3d1(
ccart_clinref_op_outerface_3d(`double complex',1,2,0)dnl
c
c      subroutine cartclinrefoutfacecplx3d2(
ccart_clinref_op_outerface_3d(`double complex',2,0,1)dnl
c
c***********************************************************************
c Conservative linear interpolation for 3d side-centered double data
c***********************************************************************
c
      subroutine cartclinrefsidedoub3d0(
cart_clinref_op_side_3d(`double precision',0,1,2)dnl
c
      subroutine cartclinrefsidedoub3d1(
cart_clinref_op_side_3d(`double precision',1,2,0)dnl
c
      subroutine cartclinrefsidedoub3d2(
cart_clinref_op_side_3d(`double precision',2,0,1)dnl
c
c***********************************************************************
c Conservative linear interpolation for 3d side-centered float data
c***********************************************************************
c
      subroutine cartclinrefsideflot3d0(
cart_clinref_op_side_3d(`real',0,1,2)dnl
c
      subroutine cartclinrefsideflot3d1(
cart_clinref_op_side_3d(`real',1,2,0)dnl
c
      subroutine cartclinrefsideflot3d2(
cart_clinref_op_side_3d(`real',2,0,1)dnl
c
c***********************************************************************
c Conservative linear interpolation for 3d side-centered complex data
c***********************************************************************
c
c     subroutine cartclinrefsidecplx3d0(
ccart_clinref_op_side_3d(`double complex',0,1,2)dnl
c
c      subroutine cartclinrefsidecplx3d1(
ccart_clinref_op_side_3d(`double complex',1,2,0)dnl
c
c      subroutine cartclinrefsidecplx3d2(
ccart_clinref_op_side_3d(`double complex',2,0,1)dnl
c

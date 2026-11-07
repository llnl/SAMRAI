c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   FORTRAN routines for spatial refining of 2d patch data
c                on a regular Cartesian mesh.
c
include(FORTDIR/geom_m4cartrefineops2d.i)dnl
c
c***********************************************************************
c Linear interpolation for 2d cell-centered double data
c***********************************************************************
c
      subroutine cartlinrefcelldoub2d(
cart_linref_op_cell_2d(`double precision')dnl
c
c***********************************************************************
c Linear interpolation for 2d cell-centered float data
c***********************************************************************
c
      subroutine cartlinrefcellflot2d(
cart_linref_op_cell_2d(`real')dnl
c
c***********************************************************************
c Linear interpolation for 2d cell-centered complex data
c***********************************************************************
c
      subroutine cartlinrefcellcplx2d(
cart_linref_op_cell_2d(`double complex')dnl
c
c***********************************************************************
c Conservative linear interpolation for 2d cell-centered double data
c***********************************************************************
c
      subroutine cartclinrefcelldoub2d(
cart_clinref_op_cell_2d(`double precision')dnl
c
c***********************************************************************
c Conservative linear interpolation for 2d cell-centered float data
c***********************************************************************
c
      subroutine cartclinrefcellflot2d(
cart_clinref_op_cell_2d(`real')dnl
c
c***********************************************************************
c Conservative linear interpolation for 2d cell-centered complex data
c***********************************************************************
c
      subroutine cartclinrefcellcplx2d(
cart_clinref_op_cell_2d(`double complex')dnl
c
c***********************************************************************
c Conservative linear interpolation for 2d edge-centered double data
c***********************************************************************
c
      subroutine cartclinrefedgedoub2d0(
cart_clinref_op_edge_2d(`double precision',0,1)dnl
c
      subroutine cartclinrefedgedoub2d1(
cart_clinref_op_edge_2d(`double precision',1,0)dnl
c
c***********************************************************************
c Conservative linear interpolation for 2d edge-centered float data
c***********************************************************************
c
      subroutine cartclinrefedgeflot2d0(
cart_clinref_op_edge_2d(`real',0,1)dnl
c
      subroutine cartclinrefedgeflot2d1(
cart_clinref_op_edge_2d(`real',1,0)dnl
c
c***********************************************************************
c Conservative linear interpolation for 2d face-centered double data
c***********************************************************************
c
      subroutine cartclinreffacedoub2d0(
cart_clinref_op_face_2d(`double precision',0,1)dnl
c
      subroutine cartclinreffacedoub2d1(
cart_clinref_op_face_2d(`double precision',1,0)dnl
c
c***********************************************************************
c Conservative linear interpolation for 2d face-centered float data
c***********************************************************************
c
      subroutine cartclinreffaceflot2d0(
cart_clinref_op_face_2d(`real',0,1)dnl
c
      subroutine cartclinreffaceflot2d1(
cart_clinref_op_face_2d(`real',1,0)dnl
c
c***********************************************************************
c Conservative linear interpolation for 2d face-centered complex data
c***********************************************************************
c
c      subroutine cartclinreffacecplx2d0(
ccart_clinref_op_face_2d(`double complex',0,1)dnl
cc
c      subroutine cartclinreffacecplx2d1(
ccart_clinref_op_face_2d(`double complex',1,0)dnl
c
c***********************************************************************
c Linear interpolation for 2d node-centered double data
c***********************************************************************
c
       subroutine cartlinrefnodedoub2d(
cart_linref_op_node_2d(`double precision')dnl
c
c***********************************************************************
c Linear interpolation for 2d node-centered float data
c***********************************************************************
c
       subroutine cartlinrefnodeflot2d(
cart_linref_op_node_2d(`real')dnl
c
c***********************************************************************
c Linear interpolation for 2d node-centered complex data
c***********************************************************************
c
       subroutine cartlinrefnodecplx2d(
cart_linref_op_node_2d(`double complex')dnl
c
c***********************************************************************
c Conservative linear interpolation for 2d outerface double data
c***********************************************************************
c
      subroutine cartclinrefoutfacedoub2d0(
cart_clinref_op_outerface_2d(`double precision',0,1)dnl
c
      subroutine cartclinrefoutfacedoub2d1(
cart_clinref_op_outerface_2d(`double precision',1,0)dnl
c
c***********************************************************************
c Conservative linear interpolation for 2d outerface float data
c***********************************************************************
c
      subroutine cartclinrefoutfaceflot2d0(
cart_clinref_op_outerface_2d(`real',0,1)dnl
c
      subroutine cartclinrefoutfaceflot2d1(
cart_clinref_op_outerface_2d(`real',1,0)dnl
c
c***********************************************************************
c Conservative linear interpolation for 2d outerface complex data
c***********************************************************************
c
c      subroutine cartclinrefoutfacecplx2d0(
ccart_clinref_op_outerface_2d(`double complex',0,1)dnl
cc
c      subroutine cartclinrefoutfacecplx2d1(
ccart_clinref_op_outerface_2d(`double complex',1,0)dnl
c
c***********************************************************************
c Conservative linear interpolation for 2d side-centered double data
c***********************************************************************
c
      subroutine cartclinrefsidedoub2d0(
cart_clinref_op_side_2d(`double precision',0,1)dnl
c
      subroutine cartclinrefsidedoub2d1(
cart_clinref_op_side_2d(`double precision',1,0)dnl
c
c***********************************************************************
c Conservative linear interpolation for 2d side-centered float data
c***********************************************************************
c
      subroutine cartclinrefsideflot2d0(
cart_clinref_op_side_2d(`real',0,1)dnl
c
      subroutine cartclinrefsideflot2d1(
cart_clinref_op_side_2d(`real',1,0)dnl
c
c***********************************************************************
c Conservative linear interpolation for 2d side-centered complex data
c***********************************************************************
c
c      subroutine cartclinrefsidecplx2d0(
ccart_clinref_op_side_2d(`double complex',0,1)dnl
cc
c      subroutine cartclinrefsidecplx2d1(
ccart_clinref_op_side_2d(`double complex',1,0)dnl
c

c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   FORTRAN routines for spatial refining of 2d patch data
c                on a regular Cartesian mesh.
c
include(PDAT_FORTDIR/pdat_m4conrefineops2d.i)dnl
c
c***********************************************************************
c Constant interpolation for 2d cell-centered double data
c***********************************************************************
c
      subroutine conrefcelldoub2d(
conref_op_cell_2d(`double precision')dnl
c
c***********************************************************************
c Constant interpolation for 2d cell-centered float data
c***********************************************************************
c
      subroutine conrefcellflot2d(
conref_op_cell_2d(`real')dnl
c
c***********************************************************************
c Constant interpolation for 2d cell-centered complex data
c***********************************************************************
c
      subroutine conrefcellcplx2d(
conref_op_cell_2d(`double complex')dnl
c
c***********************************************************************
c Constant interpolation for 2d cell-centered integer data
c***********************************************************************
c
      subroutine conrefcellintg2d(
conref_op_cell_2d(`integer')dnl
c
c***********************************************************************
c Constant interpolation for 2d edge-centered double data
c***********************************************************************
c
      subroutine conrefedgedoub2d0(
conref_op_edge_2d(`double precision',0,1)dnl
c
      subroutine conrefedgedoub2d1(
conref_op_edge_2d(`double precision',1,0)dnl
c
c***********************************************************************
c Constant interpolation for 2d edge-centered float data
c***********************************************************************
c
      subroutine conrefedgeflot2d0(
conref_op_edge_2d(`real',0,1)dnl
c
      subroutine conrefedgeflot2d1(
conref_op_edge_2d(`real',1,0)dnl
c
c***********************************************************************
c Constant interpolation for 2d edge-centered complex data
c***********************************************************************

      subroutine conrefedgecplx2d0(
conref_op_edge_2d(`double complex',0,1)dnl
c
      subroutine conrefedgecplx2d1(
conref_op_edge_2d(`double complex',1,0)dnl
c
c***********************************************************************
c Constant interpolation for 2d edge-centered integer data
c***********************************************************************
c
      subroutine conrefedgeintg2d0(
conref_op_edge_2d(`integer',0,1)dnl
c
      subroutine conrefedgeintg2d1(
conref_op_edge_2d(`integer',1,0)dnl
c
c***********************************************************************
c Constant interpolation for 2d face-centered double data
c***********************************************************************
c
      subroutine conreffacedoub2d0(
conref_op_face_2d(`double precision',0,1)dnl
c
      subroutine conreffacedoub2d1(
conref_op_face_2d(`double precision',1,0)dnl
c
c***********************************************************************
c Constant interpolation for 2d face-centered float data
c***********************************************************************
c
      subroutine conreffaceflot2d0(
conref_op_face_2d(`real',0,1)dnl
c
      subroutine conreffaceflot2d1(
conref_op_face_2d(`real',1,0)dnl
c
c***********************************************************************
c Constant interpolation for 2d face-centered complex data
c***********************************************************************

      subroutine conreffacecplx2d0(
conref_op_face_2d(`double complex',0,1)dnl
c
      subroutine conreffacecplx2d1(
conref_op_face_2d(`double complex',1,0)dnl
c
c***********************************************************************
c Constant interpolation for 2d face-centered integer data
c***********************************************************************
c
      subroutine conreffaceintg2d0(
conref_op_face_2d(`integer',0,1)dnl
c
      subroutine conreffaceintg2d1(
conref_op_face_2d(`integer',1,0)dnl
c
c***********************************************************************
c Constant interpolation for 2d outerface double data
c***********************************************************************
c
      subroutine conrefoutfacedoub2d0(
conref_op_outerface_2d(`double precision',0,1)dnl
c
      subroutine conrefoutfacedoub2d1(
conref_op_outerface_2d(`double precision',1,0)dnl
c
c***********************************************************************
c Constant interpolation for 2d outerface float data
c***********************************************************************
c
      subroutine conrefoutfaceflot2d0(
conref_op_outerface_2d(`real',0,1)dnl
c
      subroutine conrefoutfaceflot2d1(
conref_op_outerface_2d(`real',1,0)dnl
c
c***********************************************************************
c Constant interpolation for 2d outerface complex data
c***********************************************************************

      subroutine conrefoutfacecplx2d0(
conref_op_outerface_2d(`double complex',0,1)dnl
c
      subroutine conrefoutfacecplx2d1(
conref_op_outerface_2d(`double complex',1,0)dnl
c
c***********************************************************************
c Constant interpolation for 2d outerface integer data
c***********************************************************************
c
      subroutine conrefoutfaceintg2d0(
conref_op_outerface_2d(`integer',0,1)dnl
c
      subroutine conrefoutfaceintg2d1(
conref_op_outerface_2d(`integer',1,0)dnl
c
c***********************************************************************
c Constant interpolation for 2d side-centered double data
c***********************************************************************
c
      subroutine conrefsidedoub2d0(
conref_op_side_2d(`double precision',0,1)dnl
c
      subroutine conrefsidedoub2d1(
conref_op_side_2d(`double precision',1,0)dnl
c
c***********************************************************************
c Constant interpolation for 2d side-centered float data
c***********************************************************************
c
      subroutine conrefsideflot2d0(
conref_op_side_2d(`real',0,1)dnl
c
      subroutine conrefsideflot2d1(
conref_op_side_2d(`real',1,0)dnl
c
c***********************************************************************
c Constant interpolation for 2d side-centered complex data
c***********************************************************************

      subroutine conrefsidecplx2d0(
conref_op_side_2d(`double complex',0,1)dnl
c
      subroutine conrefsidecplx2d1(
conref_op_side_2d(`double complex',1,0)dnl
c
c***********************************************************************
c Constant interpolation for 2d side-centered integer data
c***********************************************************************
c
      subroutine conrefsideintg2d0(
conref_op_side_2d(`integer',0,1)dnl
c
      subroutine conrefsideintg2d1(
conref_op_side_2d(`integer',1,0)dnl
c

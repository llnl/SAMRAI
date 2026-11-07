c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   FORTRAN routines for spatial coarsening of 4d patch data
c                on a regular Cartesian mesh.
c
include(FORTDIR/geom_m4cartcoarsenops4d.i)dnl
c
c***********************************************************************
c Weighted averaging for 4d cell-centered double data
c***********************************************************************
c
      subroutine cartwgtavgcelldoub4d(
cart_wgtavg_op_cell_4d(`double precision')dnl
c
c***********************************************************************
c Weighted averaging for 4d cell-centered float data
c***********************************************************************
c
      subroutine cartwgtavgcellflot4d(
cart_wgtavg_op_cell_4d(`real')dnl
c
c***********************************************************************
c Weighted averaging for 4d cell-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgcellcplx4d(
cart_wgtavg_op_cell_4d(`double complex')dnl
c
c***********************************************************************
c Weighted averaging for 4d face-centered double data
c***********************************************************************
c
      subroutine cartwgtavgfacedoub4d0(
cart_wgtavg_op_face_4d(`double precision',0,1,2,3)dnl
c
      subroutine cartwgtavgfacedoub4d1(
cart_wgtavg_op_face_4d(`double precision',1,2,3,0)dnl
c
      subroutine cartwgtavgfacedoub4d2(
cart_wgtavg_op_face_4d(`double precision',2,3,0,1)dnl
c
      subroutine cartwgtavgfacedoub4d3(
cart_wgtavg_op_face_4d(`double precision',3,0,1,2)dnl
c
c***********************************************************************
c Weighted averaging for 4d face-centered float data
c***********************************************************************
c
      subroutine cartwgtavgfaceflot4d0(
cart_wgtavg_op_face_4d(`real',0,1,2,3)dnl
c
      subroutine cartwgtavgfaceflot4d1(
cart_wgtavg_op_face_4d(`real',1,2,3,0)dnl
c
      subroutine cartwgtavgfaceflot4d2(
cart_wgtavg_op_face_4d(`real',2,3,0,1)dnl
c
      subroutine cartwgtavgfaceflot4d3(
cart_wgtavg_op_face_4d(`real',3,0,1,2)dnl
c
c***********************************************************************
c Weighted averaging for 4d face-centered complex data
c***********************************************************************
c
      subroutine cartwgtavgfacecplx4d0(
cart_wgtavg_op_face_4d(`double complex',0,1,2,3)dnl
c
      subroutine cartwgtavgfacecplx4d1(
cart_wgtavg_op_face_4d(`double complex',1,2,3,0)dnl
c
      subroutine cartwgtavgfacecplx4d2(
cart_wgtavg_op_face_4d(`double complex',2,3,0,1)dnl
c
      subroutine cartwgtavgfacecplx4d3(
cart_wgtavg_op_face_4d(`double complex',3,0,1,2)dnl

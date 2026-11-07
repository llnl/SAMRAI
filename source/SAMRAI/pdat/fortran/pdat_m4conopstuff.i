c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for constant patchdata transfer routines.
c
define(coarsen_index,`dnl
         if ($1.lt.0) then
            $2=($1+1)/$3-1
         else
            $2=$1/$3
         endif
')dnl
define(coarsen_face_index,`dnl
         it=2*$1+$3
         if (it.le.0) then
            $2=it/(2*$3)-1
         else
            $2=(it-1)/(2*$3)
         endif
')dnl

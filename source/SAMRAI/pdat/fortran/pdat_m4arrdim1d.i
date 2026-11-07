c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for dimensioning 1d arrays in FORTRAN routines.
c
define(SAMRAICELL1d,`$1$4-$3:$2$4+$3')dnl
define(SAMRAICELL1d0G,`$1$3:$2$3')dnl
define(SAMRAICELL1dVECG,`$1$4-$3$4:$2$4+$3$4')dnl
define(SAMRAIEDGE1d,`$1$4-$3:$2$4+$3')dnl
define(SAMRAIEDGE1d0G,`$1$3:$2$3')dnl
define(SAMRAIEDGE1dVECG,`$1$4-$3$4:$2$4+$3$4')dnl
define(SAMRAIFACE1d,`$1$4-$3:$2$4+1+$3')dnl
define(SAMRAIFACE1d0G,`$1$3:$2$3+1')dnl
define(SAMRAIFACE1dVECG,`$1$4-$3$4:$2$4+1+$3$4')dnl
define(SAMRAINODE1d,`$1$4-$3:$2$4+1+$3')dnl
define(SAMRAINODE1d0G,`$1$3:$2$3+1')dnl
define(SAMRAINODE1dVECG,`$1$4-$3$4:$2$4+1+$3$4')dnl
define(SAMRAIOUTERFACE1d,`1')dnl
define(SAMRAIOUTERSIDE1d,`1')dnl
define(SAMRAIOUTERNODE1d,`1')dnl
define(SAMRAISIDE1d,`$1$4-$3:$2$4+1+$3')dnl
define(SAMRAISIDE1d0G,`$1$3:$2$3+1')dnl
define(SAMRAISIDE1dVECG,`$1$4-$3$4:$2$4+1+$3$4')dnl
define(CELL1d,`ifelse($3,`0',`SAMRAICELL1d0G($1,$2,0)',`SAMRAICELL1d($1,$2,$3,0)')')dnl 
define(EDGE1d,`ifelse($3,`0',`SAMRAIEDGE1d0G($1,$2,0)',`SAMRAIEDGE1d($1,$2,$3,0)')')dnl 
define(FACE1d,`ifelse($3,`0',`SAMRAIFACE1d0G($1,$2,0)',`SAMRAIFACE1d($1,$2,$3,0)')')dnl 
define(NODE1d,`ifelse($3,`0',`SAMRAINODE1d0G($1,$2,0)',`SAMRAINODE1d($1,$2,$3,0)')')dnl 
define(OUTERFACE1d,`SAMRAIOUTERFACE1d')dnl 
define(OUTERSIDE1d,`SAMRAIOUTERSIDE1d')dnl 
define(OUTERNODE1d,`SAMRAIOUTERNODE1d')dnl 
define(SIDE1d,`ifelse($3,`0',`SAMRAISIDE1d0G($1,$2,0)',`SAMRAISIDE1d($1,$2,$3,0)')')dnl
define(CELL1dVECG,`SAMRAICELL1dVECG($1,$2,$3,0)')dnl 
define(EDGE1dVECG,`SAMRAIEDGE1dVECG($1,$2,$3,0)')dnl 
define(FACE1dVECG,`SAMRAIFACE1dVECG($1,$2,$3,0)')dnl 
define(NODE1dVECG,`SAMRAINODE1dVECG($1,$2,$3,0)')dnl
define(SIDE1dVECG,`SAMRAIFACE1dVECG($1,$2,$3,0)')dnl 

c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for dimensioning 3d arrays in FORTRAN routines.
c
define(SAMRAICELL3d,`$1$4-$3:$2$4+$3,
     &          $1$5-$3:$2$5+$3,
     &          $1$6-$3:$2$6+$3')dnl
define(SAMRAICELL3d0G,`$1$3:$2$3,
     &          $1$4:$2$4,
     &          $1$5:$2$5')dnl
define(SAMRAICELL3dVECG,`$1$4-$3$4:$2$4+$3$4,
     &          $1$5-$3$5:$2$5+$3$5,
     &          $1$6-$3$6:$2$6+$3$6')dnl
define(SAMRAIEDGE3d0,`$1`0'-$3:$2`0'+$3,
     &          $1`1'-$3:$2`1'+1+$3,
     &          $1`2'-$3:$2`2'+1+$3')dnl
define(SAMRAIEDGE3d1,`$1`0'-$3:$2`0'+1+$3,
     &          $1`1'-$3:$2`1'+$3,
     &          $1`2'-$3:$2`2'+1+$3')dnl
define(SAMRAIEDGE3d2,`$1`0'-$3:$2`0'+1+$3,
     &          $1`1'-$3:$2`1'+1+$3,
     &          $1`2'-$3:$2`2'+$3')dnl
define(SAMRAIEDGE3d0G0,`$1`0':$2`0',
     &          $1`1':$2`1'+1,
     &          $1`2':$2`2'+1')dnl
define(SAMRAIEDGE3d0G1,`$1`0':$2`0'+1,
     &          $1`1':$2`1',
     &          $1`2':$2`2'+1')dnl
define(SAMRAIEDGE3d0G2,`$1`0':$2`0'+1,
     &          $1`1':$2`1'+1,
     &          $1`2':$2`2'')dnl
define(SAMRAIEDGE3d0VECG,`$1`0'-$3`0':$2`0'+$3`0',
     &          $1`1'-$3`1':$2`1'+1+$3`1',
     &          $1`2'-$3`2':$2`2'+1+$3`2'')dnl
define(SAMRAIEDGE3d1VECG,`$1`0'-$3`0':$2`0'+1+$3`0',
     &          $1`1'-$3`1':$2`1'+$3`1',
     &          $1`2'-$3`2':$2`2'+1+$3`2'')dnl
define(SAMRAIEDGE3d2VECG,`$1`0'-$3`0':$2`0'+1+$3`0',
     &          $1`1'-$3`1':$2`1'+1+$3`1',
     &          $1`2'-$3`2':$2`2'+$3`2'')dnl
define(SAMRAIFACE3d,`$1$4-$3:$2$4+1+$3,
     &          $1$5-$3:$2$5+$3,
     &          $1$6-$3:$2$6+$3')dnl
define(SAMRAIFACE3d0G,`$1$3:$2$3+1,
     &          $1$4:$2$4,
     &          $1$5:$2$5')dnl
define(SAMRAIFACE3dVECG,`$1$4-$3$4:$2$4+1+$3$4,
     &          $1$5-$3$5:$2$5+$3$5,
     &          $1$6-$3$6:$2$6+$3$6')dnl
define(SAMRAINODE3d,`$1$4-$3:$2$4+1+$3,
     &          $1$5-$3:$2$5+1+$3,
     &          $1$6-$3:$2$6+1+$3')dnl
define(SAMRAINODE3d0G,`$1$3:$2$3+1,
     &          $1$4:$2$4+1,
     &          $1$5:$2$5+1')dnl
define(SAMRAINODE3dVECG,`$1$4-$3$4:$2$4+1+$3$4,
     &          $1$5-$3$5:$2$5+1+$3$5,
     &          $1$6-$3$6:$2$6+1+$3$6')dnl
define(SAMRAIOUTERFACE3d,`$1$4-$3:$2$4+$3,
     &          $1$5-$3:$2$5+$3')dnl
define(SAMRAIOUTERFACE3d0G,`$1$3:$2$3,
     &          $1$4:$2$4')dnl
define(SAMRAIOUTERSIDE3d,`$1$4-$3:$2$4+$3,
     &          $1$5-$3:$2$5+$3')dnl
define(SAMRAIOUTERSIDE3d0G,`$1$3:$2$3,
     &          $1$4:$2$4')dnl
define(SAMRAIOUTERNODE3d0G,`$1$3+1:$2$3,
     &          $1$4+1:$2$4')dnl
define(SAMRAIOUTERNODE3d1G,`$1$3:$2$3+1,
     &          $1$4+1:$2$4')dnl
define(SAMRAIOUTERNODE3d2G,`$1$3:$2$3+1,
     &          $1$4:$2$4+1')dnl
define(SAMRAISIDE3d0,`$1`0'-$3:$2`0'+1+$3,
     &          $1`1'-$3:$2`1'+$3,
     &          $1`2'-$3:$2`2'+$3')dnl
define(SAMRAISIDE3d1,`$1`0'-$3:$2`0'+$3,
     &          $1`1'-$3:$2`1'+1+$3,
     &          $1`2'-$3:$2`2'+$3')dnl
define(SAMRAISIDE3d2,`$1`0'-$3:$2`0'+$3,
     &          $1`1'-$3:$2`1'+$3,
     &          $1`2'-$3:$2`2'+1+$3')dnl
define(SAMRAISIDE3d0G0,`$1`0':$2`0'+1,
     &          $1`1':$2`1',
     &          $1`2':$2`2'')dnl
define(SAMRAISIDE3d0G1,`$1`0':$2`0',
     &          $1`1':$2`1'+1,
     &          $1`2':$2`2'')dnl
define(SAMRAISIDE3d0G2,`$1`0':$2`0',
     &          $1`1':$2`1',
     &          $1`2':$2`2'+1')dnl
define(SAMRAISIDE3d0VECG,`$1`0'-$3`0':$2`0'+1+$3`0',
     &          $1`1'-$3`1':$2`1'+$3`1',
     &          $1`2'-$3`2':$2`2'+$3`2'')dnl
define(SAMRAISIDE3d1VECG,`$1`0'-$3`0':$2`0'+$3`0',
     &          $1`1'-$3`1':$2`1'+1+$3`1',
     &          $1`2'-$3`2':$2`2'+$3`2'')dnl
define(SAMRAISIDE3d2VECG,`$1`0'-$3`0':$2`0'+$3`0',
     &          $1`1'-$3`1':$2`1'+$3`1',
     &          $1`2'-$3`2':$2`2'+1+$3`2'')dnl
define(CELL3d,`ifelse($3,`0',`SAMRAICELL3d0G($1,$2,0,1,2)',`SAMRAICELL3d($1,$2,$3,0,1,2)')')dnl 
define(EDGE3d0,`ifelse($3,`0',`SAMRAIEDGE3d0G0($1,$2)',`SAMRAIEDGE3d0($1,$2,$3)')')dnl 
define(EDGE3d1,`ifelse($3,`0',`SAMRAIEDGE3d0G1($1,$2)',`SAMRAIEDGE3d1($1,$2,$3)')')dnl 
define(EDGE3d2,`ifelse($3,`0',`SAMRAIEDGE3d0G2($1,$2)',`SAMRAIEDGE3d2($1,$2,$3)')')dnl 
define(FACE3d0,`ifelse($3,`0',`SAMRAIFACE3d0G($1,$2,0,1,2)',`SAMRAIFACE3d($1,$2,$3,0,1,2)')')dnl 
define(FACE3d1,`ifelse($3,`0',`SAMRAIFACE3d0G($1,$2,1,2,0)',`SAMRAIFACE3d($1,$2,$3,1,2,0)')')dnl 
define(FACE3d2,`ifelse($3,`0',`SAMRAIFACE3d0G($1,$2,2,0,1)',`SAMRAIFACE3d($1,$2,$3,2,0,1)')')dnl
define(NODE3d,`ifelse($3,`0',`SAMRAINODE3d0G($1,$2,0,1,2)',`SAMRAINODE3d($1,$2,$3,0,1,2)')')dnl 
define(OUTERFACE3d0,`ifelse($3,`0',`SAMRAIOUTERFACE3d0G($1,$2,1,2)',`SAMRAIOUTERFACE3d($1,$2,$3,1,2)')')dnl 
define(OUTERFACE3d1,`ifelse($3,`0',`SAMRAIOUTERFACE3d0G($1,$2,2,0)',`SAMRAIOUTERFACE3d($1,$2,$3,2,0)')')dnl 
define(OUTERFACE3d2,`ifelse($3,`0',`SAMRAIOUTERFACE3d0G($1,$2,0,1)',`SAMRAIOUTERFACE3d($1,$2,$3,0,1)')')dnl
define(OUTERSIDE3d0,`ifelse($3,`0',`SAMRAIOUTERSIDE3d0G($1,$2,1,2)',`SAMRAIOUTERSIDE3d($1,$2,$3,1,2)')')dnl 
define(OUTERSIDE3d1,`ifelse($3,`0',`SAMRAIOUTERSIDE3d0G($1,$2,0,2)',`SAMRAIOUTERSIDE3d($1,$2,$3,0,2)')')dnl 
define(OUTERSIDE3d2,`ifelse($3,`0',`SAMRAIOUTERSIDE3d0G($1,$2,0,1)',`SAMRAIOUTERSIDE3d($1,$2,$3,0,1)')')dnl
define(OUTERNODE3d0,`SAMRAIOUTERNODE3d0G($1,$2,1,2)')dnl
define(OUTERNODE3d1,`SAMRAIOUTERNODE3d1G($1,$2,0,2)')dnl
define(OUTERNODE3d2,`SAMRAIOUTERNODE3d2G($1,$2,0,1)')dnl
define(SIDE3d0,`ifelse($3,`0',`SAMRAISIDE3d0G0($1,$2)',`SAMRAISIDE3d0($1,$2,$3)')')dnl 
define(SIDE3d1,`ifelse($3,`0',`SAMRAISIDE3d0G1($1,$2)',`SAMRAISIDE3d1($1,$2,$3)')')dnl 
define(SIDE3d2,`ifelse($3,`0',`SAMRAISIDE3d0G2($1,$2)',`SAMRAISIDE3d2($1,$2,$3)')')dnl 
define(CELL3dVECG,`SAMRAICELL3dVECG($1,$2,$3,0,1,2)')dnl 
define(EDGE3d0VECG,`SAMRAIEDGE3d0VECG($1,$2,$3)')dnl
define(EDGE3d1VECG,`SAMRAIEDGE3d1VECG($1,$2,$3)')dnl
define(EDGE3d2VECG,`SAMRAIEDGE3d2VECG($1,$2,$3)')dnl 
define(FACE3d0VECG,`SAMRAIFACE3dVECG($1,$2,$3,0,1,2)')dnl 
define(FACE3d1VECG,`SAMRAIFACE3dVECG($1,$2,$3,1,2,0)')dnl 
define(FACE3d2VECG,`SAMRAIFACE3dVECG($1,$2,$3,2,0,1)')dnl 
define(NODE3dVECG,`SAMRAINODE3dVECG($1,$2,$3,0,1,2)')dnl
define(SIDE3d0VECG,`SAMRAISIDE3d0VECG($1,$2,$3)')dnl
define(SIDE3d1VECG,`SAMRAISIDE3d1VECG($1,$2,$3)')dnl
define(SIDE3d2VECG,`SAMRAISIDE3d2VECG($1,$2,$3)')dnl 

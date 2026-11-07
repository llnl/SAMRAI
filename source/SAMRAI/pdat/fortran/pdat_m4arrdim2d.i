c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for dimensioning 2d arrays in FORTRAN routines.
c
define(SAMRAICELL2d,`$1$4-$3:$2$4+$3,
     &          $1$5-$3:$2$5+$3')dnl
define(SAMRAICELL2d0G,`$1$3:$2$3,
     &          $1$4:$2$4')dnl
define(SAMRAICELL2dVECG,`$1$4-$3$4:$2$4+$3$4,
     &          $1$5-$3$5:$2$5+$3$5')dnl
define(SAMRAIEDGE2d0,`$1`0'-$3:$2`0'+$3,
     &          $1`1'-$3:$2`1'+1+$3')dnl
define(SAMRAIEDGE2d1,`$1`0'-$3:$2`0'+1+$3,
     &          $1`1'-$3:$2`1'+$3')dnl
define(SAMRAIEDGE2d0G0,`$1`0':$2`0',
     &          $1`1':$2`1'+1')dnl
define(SAMRAIEDGE2d0G1,`$1`0':$2`0'+1,
     &          $1`1':$2`1'')dnl
define(SAMRAIEDGE2d0VECG,`$1`0'-$3`0':$2`0'+$3`0',
     &          $1`1'-$3`1':$2`1'+1+$3`1'')dnl
define(SAMRAIEDGE2d1VECG,`$1`0'-$3`0':$2`0'+1+$3`0',
     &          $1`1'-$3`1':$2`1'+$3`1'')dnl
define(SAMRAIFACE2d,`$1$4-$3:$2$4+1+$3,
     &          $1$5-$3:$2$5+$3')dnl
define(SAMRAIFACE2d0G,`$1$3:$2$3+1,
     &          $1$4:$2$4')dnl
define(SAMRAIFACE2dVECG,`$1$4-$3$4:$2$4+1+$3$4,
     &          $1$5-$3$5:$2$5+$3$5')dnl
define(SAMRAINODE2d,`$1$4-$3:$2$4+1+$3,
     &          $1$5-$3:$2$5+1+$3')dnl
define(SAMRAINODE2d0G,`$1$3:$2$3+1,
     &          $1$4:$2$4+1')dnl
define(SAMRAINODE2dVECG,`$1$4-$3$4:$2$4+1+$3$4,
     &          $1$5-$3$5:$2$5+1+$3$5')dnl
define(SAMRAIOUTERFACE2d,`$1$4-$3:$2$4+$3')dnl
define(SAMRAIOUTERFACE2d0G,`$1$3:$2$3')dnl
define(SAMRAIOUTERSIDE2d,`$1$4-$3:$2$4+$3')dnl
define(SAMRAIOUTERSIDE2d0G,`$1$3:$2$3')dnl
define(SAMRAIOUTERNODE2d0G,`$1$3+1:$2$3')dnl
define(SAMRAIOUTERNODE2d1G,`$1$3:$2$3+1')dnl
define(SAMRAISIDE2d0,`$1`0'-$3:$2`0'+1+$3,
     &          $1`1'-$3:$2`1'+$3')dnl
define(SAMRAISIDE2d1,`$1`0'-$3:$2`0'+$3,
     &          $1`1'-$3:$2`1'+1+$3')dnl
define(SAMRAISIDE2d0G0,`$1`0':$2`0'+1,
     &          $1`1':$2`1'')dnl
define(SAMRAISIDE2d0G1,`$1`0':$2`0',
     &          $1`1':$2`1'+1')dnl
define(SAMRAISIDE2d0VECG,`$1`0'-$3`0':$2`0'+1+$3`0',
     &          $1`1'-$3`1':$2`1'+$3`1'')dnl
define(SAMRAISIDE2d1VECG,`$1`0'-$3`0':$2`0'+$3`0',
     &          $1`1'-$3`1':$2`1'+1+$3`1'')dnl
define(CELL2d,`ifelse($3,`0',`SAMRAICELL2d0G($1,$2,0,1)',`SAMRAICELL2d($1,$2,$3,
0,1)')')dnl
define(EDGE2d0,`ifelse($3,`0',`SAMRAIEDGE2d0G0($1,$2)',`SAMRAIEDGE2d0($1,$2,$3)')')dnl 
define(EDGE2d1,`ifelse($3,`0',`SAMRAIEDGE2d0G1($1,$2)',`SAMRAIEDGE2d1($1,$2,$3)')')dnl 
define(FACE2d0,`ifelse($3,`0',`SAMRAIFACE2d0G($1,$2,0,1)',`SAMRAIFACE2d($1,$2,$3,0,1)')')dnl 
define(FACE2d1,`ifelse($3,`0',`SAMRAIFACE2d0G($1,$2,1,0)',`SAMRAIFACE2d($1,$2,$3,1,0)')')dnl 
define(NODE2d,`ifelse($3,`0',`SAMRAINODE2d0G($1,$2,0,1)',`SAMRAINODE2d($1,$2,$3,0,1)')')dnl 
define(OUTERFACE2d0,`ifelse($3,`0',`SAMRAIOUTERFACE2d0G($1,$2,1)',`SAMRAIOUTERFACE2d($1,$2,$3,1)')')dnl 
define(OUTERFACE2d1,`ifelse($3,`0',`SAMRAIOUTERFACE2d0G($1,$2,0)',`SAMRAIOUTERFACE2d($1,$2,$3,0)')')dnl
define(OUTERSIDE2d0,`ifelse($3,`0',`SAMRAIOUTERSIDE2d0G($1,$2,1)',`SAMRAIOUTERSIDE2d($1,$2,$3,1)')')dnl 
define(OUTERSIDE2d1,`ifelse($3,`0',`SAMRAIOUTERSIDE2d0G($1,$2,0)',`SAMRAIOUTERSIDE2d($1,$2,$3,0)')')dnl
define(OUTERNODE2d0,`SAMRAIOUTERNODE2d0G($1,$2,1)')dnl 
define(OUTERNODE2d1,`SAMRAIOUTERNODE2d1G($1,$2,0)')dnl 
define(SIDE2d0,`ifelse($3,`0',`SAMRAISIDE2d0G0($1,$2)',`SAMRAISIDE2d0($1,$2,$3)')')dnl 
define(SIDE2d1,`ifelse($3,`0',`SAMRAISIDE2d0G1($1,$2)',`SAMRAISIDE2d1($1,$2,$3)')')dnl 
define(CELL2dVECG,`SAMRAICELL2dVECG($1,$2,$3,0,1)')dnl
define(EDGE2d0VECG,`SAMRAIEDGE2d0VECG($1,$2,$3)')dnl
define(EDGE2d1VECG,`SAMRAIEDGE2d1VECG($1,$2,$3)')dnl
define(FACE2d0VECG,`SAMRAIFACE2dVECG($1,$2,$3,0,1)')dnl 
define(FACE2d1VECG,`SAMRAIFACE2dVECG($1,$2,$3,1,0)')dnl 
define(NODE2dVECG,`SAMRAINODE2dVECG($1,$2,$3,0,1)')dnl
define(SIDE2d0VECG,`SAMRAISIDE2d0VECG($1,$2,$3)')dnl
define(SIDE2d1VECG,`SAMRAISIDE2d1VECG($1,$2,$3)')dnl

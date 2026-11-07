c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for dimensioning 4d arrays in FORTRAN routines.
c
define(SAMRAICELL4d,`$1$4-$3:$2$4+$3,
     &          $1$5-$3:$2$5+$3,
     &          $1$6-$3:$2$6+$3,
     &          $1$7-$3:$2$7+$3')dnl
define(SAMRAICELL4d0G,`$1$3:$2$3,
     &          $1$4:$2$4,
     &          $1$5:$2$5,
     &          $1$6:$2$6')dnl
define(SAMRAICELL4dVECG,`$1$4-$3$4:$2$4+$3$4,
     &          $1$5-$3$5:$2$5+$3$5,
     &          $1$6-$3$6:$2$6+$3$6,
     &          $1$7-$3$7:$2$7+$3$7')dnl
define(SAMRAIEDGE4d0,`$1`0'-$3:$2`0'+$3,
     &          $1`1'-$3:$2`1'+1+$3,
     &          $1`2'-$3:$2`2'+1+$3,
     &          $1`3'-$3:$2`3'+1+$3')dnl
define(SAMRAIEDGE4d1,`$1`0'-$3:$2`0'+1+$3,
     &          $1`1'-$3:$2`1'+$3,
     &          $1`2'-$3:$2`2'+1+$3,
     &          $1`3'-$3:$2`3'+1+$3')dnl
define(SAMRAIEDGE4d2,`$1`0'-$3:$2`0'+1+$3,
     &          $1`1'-$3:$2`1'+1+$3,
     &          $1`2'-$3:$2`2'+$3,
     &          $1`3'-$3:$2`3'+1+$3')dnl
define(SAMRAIEDGE4d3,`$1`0'-$3:$2`0'+1+$3,
     &          $1`1'-$3:$2`1'+1+$3,
     &          $1`2'-$3:$2`2'+1+$3,
     &          $1`3'-$3:$2`3'+$3')dnl
define(SAMRAIEDGE4d0G0,`$1`0':$2`0',
     &          $1`1':$2`1'+1,
     &          $1`2':$2`2'+1,
     &          $1`3':$2`3'+1')dnl
define(SAMRAIEDGE4d0G1,`$1`0':$2`0'+1,
     &          $1`1':$2`1',
     &          $1`2':$2`2'+1,
     &          $1`3':$2`3'+1')dnl
define(SAMRAIEDGE4d0G2,`$1`0':$2`0'+1,
     &          $1`1':$2`1'+1,
     &          $1`2':$2`2',
     &          $1`3':$2`3'+1')dnl
define(SAMRAIEDGE4d0G3,`$1`0':$2`0'+1,
     &          $1`1':$2`1'+1,
     &          $1`2':$2`2'+1,
     &          $1`3':$2`3'')dnl
define(SAMRAIEDGE4d0VECG,`$1`0'-$3`0':$2`0'+$3`0',
     &          $1`1'-$3`1':$2`1'+1+$3`1',
     &          $1`2'-$3`2':$2`2'+1+$3`2',
     &          $1`3'-$3`3':$2`3'+1+$3`3'')dnl
define(SAMRAIEDGE4d1VECG,`$1`0'-$3`0':$2`0'+1+$3`0',
     &          $1`1'-$3`1':$2`1'+$3`1',
     &          $1`2'-$3`2':$2`2'+1+$3`2',
     &          $1`3'-$3`3':$2`3'+1+$3`3'')dnl
define(SAMRAIEDGE4d2VECG,`$1`0'-$3`0':$2`0'+1+$3`0',
     &          $1`1'-$3`1':$2`1'+1+$3`1',
     &          $1`2'-$3`2':$2`2'+$3`2',
     &          $1`3'-$3`3':$2`3'+1+$3`3'')dnl
define(SAMRAIEDGE4d3VECG,`$1`0'-$3`0':$2`0'+1+$3`0',
     &          $1`1'-$3`1':$2`1'+1+$3`1',
     &          $1`2'-$3`2':$2`2'+1+$3`2',
     &          $1`3'-$3`3':$2`3'+$3`3'')dnl
define(SAMRAIFACE4d,`$1$4-$3:$2$4+1+$3,
     &          $1$5-$3:$2$5+$3,
     &          $1$6-$3:$2$6+$3,
     &          $1$7-$3:$2$7+$3')dnl
define(SAMRAIFACE4d0G,`$1$3:$2$3+1,
     &          $1$4:$2$4,
     &          $1$5:$2$5,
     &          $1$6:$2$6')dnl
define(SAMRAIFACE4dVECG,`$1$4-$3$4:$2$4+1+$3$4,
     &          $1$5-$3$5:$2$5+$3$5,
     &          $1$6-$3$6:$2$6+$3$6,
     &          $1$7-$3$7:$2$7+$3$7')dnl
define(SAMRAINODE4d,`$1$4-$3:$2$4+1+$3,
     &          $1$5-$3:$2$5+1+$3,
     &          $1$6-$3:$2$6+1+$3,
     &          $1$7-$3:$2$7+1+$3')dnl
define(SAMRAINODE4d0G,`$1$3:$2$3+1,
     &          $1$4:$2$4+1,
     &          $1$5:$2$5+1,
     &          $1$6:$2$6+1')dnl
define(SAMRAINODE4dVECG,`$1$4-$3$4:$2$4+1+$3$4,
     &          $1$5-$3$5:$2$5+1+$3$5,
     &          $1$6-$3$6:$2$6+1+$3$6,
     &          $1$7-$3$7:$2$7+1+$3$7')dnl
define(SAMRAIOUTERFACE4d,`$1$4-$3:$2$4+$3,
     &          $1$5-$3:$2$5+$3,
     &          $1$6-$3:$2$6+$3')dnl
define(SAMRAIOUTERFACE4d0G,`$1$3:$2$3,
     &          $1$4:$2$4,
     &          $1$5:$2$5')dnl
define(SAMRAIOUTERSIDE4d,`$1$4-$3:$2$4+$3,
     &          $1$5-$3:$2$5+$3,
     &          $1$6-$3:$2$6+$3')dnl
define(SAMRAIOUTERSIDE4d0G,`$1$3:$2$3,
     &          $1$4:$2$4,
     &          $1$5:$2$5')dnl
define(SAMRAIOUTERNODE4d0G,`$1$3+1:$2$3,
     &          $1$4+1:$2$4,
     &          $1$5+1:$2$5')dnl
define(SAMRAIOUTERNODE4d1G,`$1$3:$2$3+1,
     &          $1$4+1:$2$4,
     &          $1$5+1:$2$5')dnl
define(SAMRAIOUTERNODE4d2G,`$1$3:$2$3+1,
     &          $1$4:$2$4+1,
     &          $1$5:$2$5+1')dnl
define(SAMRAIOUTERNODE4d3G,`$1$3:$2$3+1,
     &          $1$4:$2$4+1,
     &          $1$5:$2$5+1')dnl
define(SAMRAISIDE4d0,`$1`0'-$3:$2`0'+1+$3,
     &          $1`1'-$3:$2`1'+$3,
     &          $1`2'-$3:$2`2'+$3,
     &          $1`3'-$3:$2`3'+$3')dnl
define(SAMRAISIDE4d1,`$1`0'-$3:$2`0'+$3,
     &          $1`1'-$3:$2`1'+1+$3,
     &          $1`2'-$3:$2`2'+$3,
     &          $1`3'-$3:$2`3'+$3')dnl
define(SAMRAISIDE4d2,`$1`0'-$3:$2`0'+$3,
     &          $1`1'-$3:$2`1'+$3,
     &          $1`2'-$3:$2`2'+1+$3,
     &          $1`3'-$3:$2`3'+$3')dnl
define(SAMRAISIDE4d3,`$1`0'-$3:$2`0'+$3,
     &          $1`1'-$3:$2`1'+$3,
     &          $1`2'-$3:$2`2'+$3,
     &          $1`3'-$3:$2`3'+1+$3')dnl
define(SAMRAISIDE4d0G0,`$1`0':$2`0'+1,
     &          $1`1':$2`1',
     &          $1`2':$2`2',
     &          $1`3':$2`3'')dnl
define(SAMRAISIDE4d0G1,`$1`0':$2`0',
     &          $1`1':$2`1'+1,
     &          $1`2':$2`2',
     &          $1`3':$2`3'')dnl
define(SAMRAISIDE4d0G2,`$1`0':$2`0',
     &          $1`1':$2`1',
     &          $1`2':$2`2'+1,
     &          $1`3':$2`3'')dnl
define(SAMRAISIDE4d0G3,`$1`0':$2`0',
     &          $1`1':$2`1',
     &          $1`2':$2`2',
     &          $1`3':$2`3'+1')dnl
define(SAMRAISIDE4d0VECG,`$1`0'-$3`0':$2`0'+1+$3`0',
     &          $1`1'-$3`1':$2`1'+$3`1',
     &          $1`2'-$3`2':$2`2'+$3`2',
     &          $1`3'-$3`3':$2`3'+$3`3'')dnl
define(SAMRAISIDE4d1VECG,`$1`0'-$3`0':$2`0'+$3`0',
     &          $1`1'-$3`1':$2`1'+1+$3`1',
     &          $1`2'-$3`2':$2`2'+$3`2',
     &          $1`3'-$3`3':$2`3'+$3`3'')dnl
define(SAMRAISIDE4d2VECG,`$1`0'-$3`0':$2`0'+$3`0',
     &          $1`1'-$3`1':$2`1'+$3`1',
     &          $1`2'-$3`2':$2`2'+1+$3`2',
     &          $1`3'-$3`3':$2`3'+$3`3'')dnl
define(SAMRAISIDE4d3VECG,`$1`0'-$3`0':$2`0'+$3`0',
     &          $1`1'-$3`1':$2`1'+$3`1',
     &          $1`2'-$3`2':$2`2'+$3`2',
     &          $1`3'-$3`3':$2`3'+1+$3`3'')dnl
define(CELL4d,`ifelse($3,`0',`SAMRAICELL4d0G($1,$2,0,1,2,3)',`SAMRAICELL4d($1,$2,$3,0,1,2,3)')')dnl 
define(EDGE4d0,`ifelse($3,`0',`SAMRAIEDGE4d0G0($1,$2)',`SAMRAIEDGE4d0($1,$2,$3)')')dnl 
define(EDGE4d1,`ifelse($3,`0',`SAMRAIEDGE4d0G1($1,$2)',`SAMRAIEDGE4d1($1,$2,$3)')')dnl 
define(EDGE4d2,`ifelse($3,`0',`SAMRAIEDGE4d0G2($1,$2)',`SAMRAIEDGE4d2($1,$2,$3)')')dnl 
define(EDGE4d3,`ifelse($3,`0',`SAMRAIEDGE4d0G3($1,$2)',`SAMRAIEDGE4d3($1,$2,$3)')')dnl 
define(FACE4d0,`ifelse($3,`0',`SAMRAIFACE4d0G($1,$2,0,1,2,3)',`SAMRAIFACE4d($1,$2,$3,0,1,2,3)')')dnl 
define(FACE4d1,`ifelse($3,`0',`SAMRAIFACE4d0G($1,$2,1,2,3,0)',`SAMRAIFACE4d($1,$2,$3,1,2,3,0)')')dnl 
define(FACE4d2,`ifelse($3,`0',`SAMRAIFACE4d0G($1,$2,2,3,0,1)',`SAMRAIFACE4d($1,$2,$3,2,3,0,1)')')dnl
define(FACE4d3,`ifelse($3,`0',`SAMRAIFACE4d0G($1,$2,3,0,1,2)',`SAMRAIFACE4d($1,$2,$3,3,0,1,2)')')dnl
define(NODE4d,`ifelse($3,`0',`SAMRAINODE4d0G($1,$2,0,1,2,3)',`SAMRAINODE4d($1,$2,$3,0,1,2,3)')')dnl 
define(OUTERFACE4d0,`ifelse($3,`0',`SAMRAIOUTERFACE4d0G($1,$2,1,2,3)',`SAMRAIOUTERFACE4d($1,$2,$3,1,2,3)')')dnl 
define(OUTERFACE4d1,`ifelse($3,`0',`SAMRAIOUTERFACE4d0G($1,$2,2,3,0)',`SAMRAIOUTERFACE4d($1,$2,$3,2,3,0)')')dnl 
define(OUTERFACE4d2,`ifelse($3,`0',`SAMRAIOUTERFACE4d0G($1,$2,3,0,1)',`SAMRAIOUTERFACE4d($1,$2,$3,3,0,1)')')dnl
define(OUTERFACE4d3,`ifelse($3,`0',`SAMRAIOUTERFACE4d0G($1,$2,0,1,2)',`SAMRAIOUTERFACE4d($1,$2,$3,0,1,2)')')dnl
define(OUTERSIDE4d0,`ifelse($3,`0',`SAMRAIOUTERSIDE4d0G($1,$2,1,2,3)',`SAMRAIOUTERSIDE4d($1,$2,$3,1,2,3)')')dnl 
define(OUTERSIDE4d1,`ifelse($3,`0',`SAMRAIOUTERSIDE4d0G($1,$2,0,2,3)',`SAMRAIOUTERSIDE4d($1,$2,$3,0,2,3)')')dnl 
define(OUTERSIDE4d2,`ifelse($3,`0',`SAMRAIOUTERSIDE4d0G($1,$2,0,1,3)',`SAMRAIOUTERSIDE4d($1,$2,$3,0,1,3)')')dnl
define(OUTERSIDE4d3,`ifelse($3,`0',`SAMRAIOUTERSIDE4d0G($1,$2,0,1,2)',`SAMRAIOUTERSIDE4d($1,$2,$3,0,1,2)')')dnl
define(OUTERNODE4d0,`SAMRAIOUTERNODE4d0G($1,$2,1,2,3)')dnl
define(OUTERNODE4d1,`SAMRAIOUTERNODE4d1G($1,$2,0,2,3)')dnl
define(OUTERNODE4d2,`SAMRAIOUTERNODE4d2G($1,$2,0,1,3)')dnl
define(OUTERNODE4d3,`SAMRAIOUTERNODE4d3G($1,$2,0,1,2)')dnl
define(SIDE4d0,`ifelse($3,`0',`SAMRAISIDE4d0G0($1,$2)',`SAMRAISIDE4d0($1,$2,$3)')')dnl 
define(SIDE4d1,`ifelse($3,`0',`SAMRAISIDE4d0G1($1,$2)',`SAMRAISIDE4d1($1,$2,$3)')')dnl 
define(SIDE4d2,`ifelse($3,`0',`SAMRAISIDE4d0G2($1,$2)',`SAMRAISIDE4d2($1,$2,$3)')')dnl 
define(SIDE4d3,`ifelse($3,`0',`SAMRAISIDE4d0G3($1,$2)',`SAMRAISIDE4d3($1,$2,$3)')')dnl 
define(CELL4dVECG,`SAMRAICELL4dVECG($1,$2,$3,0,1,2,3)')dnl 
define(EDGE4d0VECG,`SAMRAIEDGE4d0VECG($1,$2,$3)')dnl
define(EDGE4d1VECG,`SAMRAIEDGE4d1VECG($1,$2,$3)')dnl
define(EDGE4d2VECG,`SAMRAIEDGE4d2VECG($1,$2,$3)')dnl 
define(EDGE4d3VECG,`SAMRAIEDGE4d3VECG($1,$2,$3)')dnl 
define(FACE4d0VECG,`SAMRAIFACE4dVECG($1,$2,$3,0,1,2,3)')dnl 
define(FACE4d1VECG,`SAMRAIFACE4dVECG($1,$2,$3,1,2,3,0)')dnl 
define(FACE4d2VECG,`SAMRAIFACE4dVECG($1,$2,$3,2,3,0,1)')dnl 
define(FACE4d3VECG,`SAMRAIFACE4dVECG($1,$2,$3,3,0,1,2)')dnl 
define(NODE4dVECG,`SAMRAINODE4dVECG($1,$2,$3,0,1,2,3)')dnl
define(SIDE4d0VECG,`SAMRAISIDE4d0VECG($1,$2,$3)')dnl
define(SIDE4d1VECG,`SAMRAISIDE4d1VECG($1,$2,$3)')dnl
define(SIDE4d2VECG,`SAMRAISIDE4d2VECG($1,$2,$3)')dnl 
define(SIDE4d3VECG,`SAMRAISIDE4d3VECG($1,$2,$3)')dnl 

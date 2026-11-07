c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for 2d flux calculation.
c
define(riemann_solve,`dnl
c     write(6,*) "checking onedr sol in riemann solve  "
c     write(6,*) "         dt= ",dt
c     write(6,*) "  calculating flux$1, 1+extra_cell= ",$3
c     write(6,*) "  ic$2=",ifirst$2-1-$3,ilast$2+1+$3
c     write(6,*) "  ie$1=",ifirst$1-1-$3,ilast$1+1+1+$3
        do ic$2=ifirst$2-$3,ilast$2+$3
          do ie$1=ifirst$1-$3,ilast$1+1+$3
 
            if (advecspeed($1).ge.zero) then
               riemst= trlft$1(ie$1,ic$2)
             else
               riemst= trrgt$1(ie$1,ic$2)
             endif
 
            flux$1(ie$1,ic$2)= dt*riemst*advecspeed($1)
c           write(6,*) "   flux$1(",ie$1,ic$2,")= ",flux$1(ie$1,ic$2,1),
c    &                   flux$1(ie$1,ic$2,2),
c    &          flux$1(ie$1,ic$2,3),flux$1(ie$1,ic$2,4)
          enddo
        enddo
')dnl

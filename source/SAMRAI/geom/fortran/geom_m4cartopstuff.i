c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   m4 include file for cartesian geometry transfer routines.
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

define(node_refloop_pre,`dnl
      do ic$1=ifirstc$1,ilastc$1
         if$1=ic$1*ratio($1)
         do ir$1=0,ratio($1)
            ie$1=if$1+ir$1
            if ((ie$1.ge.filo$1).and.(ie$1.le.(fihi$1+1))) then
')dnl

define(node_refloop_post,`dnl
           endif
         end do
      end do
')dnl

define(coarse_fine_cell_deltas,`dnl
      do ir$1=0,ratio($1)-1
         deltax(ir$1,$1)=(dble(ir$1)+half)*dxf($1)-dxc($1)*half
      enddo
')dnl
define(coarse_fine_face_deltas,`dnl
      do ir$1=0,ratio($1)-1
         deltax(ir$1,$1)=dble(ir$1)*dxf($1)
      enddo
')dnl
define(muscl_limited_cell_slopes,`dnl
      do ie$1=ifirstc$1,ilastc$1+1
         diff$1(ie$1)=arrayc($2)
     &               -arrayc($3)
      enddo
      do ic$1=ifirstc$1,ilastc$1
         coef2=half*(diff$1(ic$1+1)+diff$1(ic$1))
         bound=two*min(abs(diff$1(ic$1+1)),abs(diff$1(ic$1)))
         if (diff$1(ic$1)*diff$1(ic$1+1).gt.zero) then
            slope$1($4)=sign(min(abs(coef2),bound),coef2)
     &                  /dxc($1)
         else
            slope$1($4)=zero
         endif
      enddo
')dnl
define(muscl_limited_cell_slopes_complex,`dnl
      do ie$1=ifirstc$1,ilastc$1+1
         diff$1(ie$1)=arrayc($2)
     &               -arrayc($3)
      enddo
      do ic$1=ifirstc$1,ilastc$1
         diff0real=dble(diff$1(ic$1))
         diff0imag=imag(diff$1(ic$1))
         diff1real=dble(diff$1(ic$1+1))
         diff1imag=imag(diff$1(ic$1+1))
         coef2real=half*(diff0real+diff1real)
         coef2imag=half*(diff0imag+diff1imag)
         boundreal=two*min(abs(diff1real),abs(diff0real))
         boundimag=two*min(abs(diff1imag),abs(diff0imag))
         if (diff0real*diff1real.gt.zero) then
            slopereal=sign(min(abs(coef2real),boundreal),coef2real)
     &                /dxc($1)
         else
            slopereal=zero
         endif
         if (diff0imag*diff1imag.gt.zero) then
            slopeimag=sign(min(abs(coef2imag),boundimag),coef2imag) 
     &                /dxc($1)
         else
            slopeimag=zero
         endif
         slope$1($4) = slopereal + cmplx(zero,one)*slopeimag
      enddo
')dnl
define(muscl_limited_face_slopes,`dnl
      do ic$1=ifirstc$1-1,ilastc$1+1
         diff$1(ic$1)=arrayc($2)
     &                -arrayc($3)
      enddo
      do ie$1=ifirstc$1,ilastc$1+1
         coef2=half*(diff$1(ie$1-1)+diff$1(ie$1))
         bound=two*min(abs(diff$1(ie$1-1)),abs(diff$1(ie$1)))
         if (diff$1(ie$1)*diff$1(ie$1-1).gt.zero) then
           slope$1($4)=sign(min(abs(coef2),bound),coef2)
     &                 /dxc($1) 
         else
            slope$1($4)=zero
         endif
      enddo
')dnl

c
c This file is part of the SAMRAI distribution.  For full copyright
c information, see COPYRIGHT and LICENSE.
c
c Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
c Description:   F77 routines for Cartesian 2d Robin boundary conditions.
c
c***********************************************************************
c***********************************************************************
      subroutine settype1cells2d(
     &  data, 
     &  difirst, dilast, djfirst, djlast,
     &  a, b, g,
     &  ifirst, ilast, jfirst, jlast,
     &  ibeg, iend, jbeg, jend,
     &  face, ghos, inte, location,
     &  h, zerog )
c***********************************************************************
      implicit none
      integer difirst, dilast, djfirst, djlast
      double precision data(difirst:dilast,djfirst:djlast)
      integer ifirst, ilast, jfirst, jlast
      double precision a(ifirst:ilast,jfirst:jlast)
      double precision b(ifirst:ilast,jfirst:jlast)
      double precision g(ifirst:ilast,jfirst:jlast)
      integer ibeg, iend
      integer jbeg, jend
      integer face, ghos, inte, location, zerog
      double precision h
      integer i, j
      if ( zerog .eq. 1 ) then
c        Assume g value of zero
         if ( location/2 .eq. 0 ) then
            do i=ibeg,iend
               do j=jbeg,jend
                  data(ghos,j) =
     &            ( data(inte,j)*( b(face,j) - 0.5*h*a(face,j) ) )
     &                 / ( b(face,j)+0.5*h*a(face,j) )
               enddo
            enddo
         elseif ( location/2 .eq. 1 ) then
            do i=ibeg,iend
               do j=jbeg,jend
                  data(i,ghos) =
     &            ( data(i,inte)*( b(i,face) - 0.5*h*a(i,face) ) )
     &                 / ( b(i,face)+0.5*h*a(i,face) )
               enddo
            enddo
         endif
      else
c        Assume finite g
         if ( location/2 .eq. 0 ) then
            do i=ibeg,iend
               do j=jbeg,jend
                  data(ghos,j) = ( h*g(face,j)
     &                 + data(inte,j)*( b(face,j) - 0.5*h*a(face,j) ) )
     &                 / ( b(face,j)+0.5*h*a(face,j) )
               enddo
            enddo
         elseif ( location/2 .eq. 1 ) then
            do i=ibeg,iend
               do j=jbeg,jend
                  data(i,ghos) = ( h*g(i,face)
     &                 + data(i,inte)*( b(i,face) - 0.5*h*a(i,face) ) )
     &                 / ( b(i,face)+0.5*h*a(i,face) )
               enddo
            enddo
         endif
      endif
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine settype2cells2d(
     &  data, 
     &  difirst, dilast, djfirst, djlast,
     &  lower, upper, location )
c***********************************************************************
      implicit none
      integer difirst, dilast, djfirst, djlast
      double precision data(difirst:dilast,djfirst:djlast)
      integer lower(0:1), upper(0:1), location
      integer i, j
      if ( location .eq. 0 ) then
c        min i min j node
         i = lower(0)
         j = lower(1)
         data(i,j) = -data(i+1,j+1) + ( data(i+1,j) + data(i,j+1) )
      elseif ( location .eq. 1 ) then
c        max i min j node
         i = upper(0)
         j = lower(1)
         data(i,j) = -data(i-1,j+1) + ( data(i-1,j) + data(i,j+1) )
      elseif ( location .eq. 2 ) then
c        min i max j node
         i = lower(0)
         j = upper(1)
         data(i,j) = -data(i+1,j-1) + ( data(i+1,j) + data(i,j-1) )
      elseif ( location .eq. 3 ) then
c        max i max j node
         i = upper(0)
         j = upper(1)
         data(i,j) = -data(i-1,j-1) + ( data(i,j-1) + data(i-1,j) )
      else
      endif
      return
      end
c***********************************************************************

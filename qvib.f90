!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: qvib.f

      subroutine qvib (a,b,m,w,x)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine constructs an m point particle-in-a-box
!     (midpoint) quadrature rule in the interval a < x < b.
!     -----------------------------------------------------------------
!
      dimension w(m),x(m)
!
      dx = (b-a)/m
      do k = 1,m
         w(k) = dx
         x(k) = a+(k-0.5d0)*dx
      enddo
      return
      end

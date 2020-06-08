!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: qrot.f

      subroutine qrot (m,w,x)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine constructs an m point Gauss-Legendre
!     quadrature rule in the interval -1 < x < 1.
!     -----------------------------------------------------------------
!
      dimension w(m),x(m)

!     dimension a2(m),a3(m)
      allocatable a2(:),a3(:)

      allocate (a2(m),a3(m))

!
      do j = 1,m
         a2(j) = (j+j-1.d0)/j
         a3(j) = (j-1.d0)/j
      enddo
      n = (m+1)/2
      pi = acos(-1.d0)
      do k = 1,n
         z = cos(pi*(k-0.25d0)/(m+0.5d0))
         do i = 1,6
            p2 = 0.d0
            p1 = 1.d0
            do j = 1,m
               p3 = p2
               p2 = p1
               p1 = a2(j)*z*p2-a3(j)*p3
            enddo
            p2 = m*(p2-z*p1)/(1.d0-z*z)
            z = z-p1/p2
         enddo
         x(k) = -z
         x(m+1-k) = +z
         w(k) = 2.d0/((1.d0-z*z)*p2*p2)
         w(m+1-k) = w(k)
      enddo
      return
      end

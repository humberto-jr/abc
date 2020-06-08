!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: prot.f

      subroutine prot (p,jmax,kmax,x)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine uses the algorithm in Section 6.8 of Numerical
!     Recipes to calculate an array of normalised associated Legendre
!     polynomials P(j,k,x) = sqrt(2*pi)*Y(j,k;acos(x),0).
!     -----------------------------------------------------------------
!
      parameter (jtop = 100)
      dimension p(0:jmax,0:kmax)
      dimension a(0:jtop,0:jtop),b(0:jtop,0:jtop),c(0:jtop,0:jtop)
      save init,a,b,c
      data init /0/
!
      if (init .eq. 0) then
         do j = 0,jtop
            do k = 0,j-2
               a(j,k) = (j+j-1.d0)/(j-k)
               b(j,k) = (j+k-1.d0)/(j-k)
            enddo
            fac = j+0.5d0
            c(j,0) = sqrt(fac)
            do k = 1,j
               fac = fac/((j+k)*(j-k+1))
               c(j,k) = sqrt(fac)
            enddo
         enddo
         init = 1
      endif
!
      if (jmax .gt. jtop) stop 'prot 0'
      if (abs(x) .gt. 1.d0) stop 'prot 1'
!
      do j = 0,jmax
         do k = j+1,kmax
            p(j,k) = 0.d0
         enddo
      enddo
      p(0,0) = 1.d0
      if (kmax .gt. 0) then
         sx = -sqrt((1.d0-x)*(1.d0+x))
         do k = 1,min(jmax,kmax)
            p(k,k) = (2*k-1)*sx*p(k-1,k-1)
         enddo
      endif
      if (jmax .gt. 0) then
         p(1,0) = x
         do k = 1,min(jmax-1,kmax)
            p(k+1,k) = (2*k+1)*x*p(k,k)
         enddo
      endif
      do k = 0,kmax
         do j = k+2,jmax
            p(j,k) = a(j,k)*x*p(j-1,k)-b(j,k)*p(j-2,k)
         enddo
      enddo
      do k = 0,kmax
         do j = k,jmax
            p(j,k) = c(j,k)*p(j,k)
         enddo
      enddo
      return
      end

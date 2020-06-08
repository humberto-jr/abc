!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: pvib.f

      subroutine pvib (a,x,b,c,n,m,p,mode)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine calculates m normalised vibrational
!     eigenfunctions in a < x < b (if mode = 0), or their
!     first derivatives with respect to x (if mode = 1).
!     -----------------------------------------------------------------
!
      dimension c(n,m),p(m)

!     dimension q(0:n)
      allocatable q(:)

      allocate (q(0:n))

!
!     primitive sine basis functions
!
      pi = acos(-1.d0)
      scale = sqrt(2.d0/(b-a))
      shift = pi/(b-a)
      y = shift*(x-a)
      cosy = cos(y)
      recur = 2*cosy
      if (mode .eq. 0) then
         q(0) = 0.d0
         q(1) = scale*sin(y)
         do j = 2,n
            q(j) = recur*q(j-1)-q(j-2)
         enddo
      else
         scale = shift*scale
         q(0) = scale
         q(1) = scale*cosy
         do j = 2,n
            q(j) = recur*q(j-1)-q(j-2)
         enddo
         do j = 2,n
            q(j) = j*q(j)
         enddo
      endif
!
!     contracted vibrational eigenfunctions
!
      call dgemv ('t',n,m,1.d0,c,n,q(1),1,0.d0,p,1)
      return
      end

!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: pvibv.f

      subroutine pvibv (a,x,b,c,n,m,p,l,mode)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     Vector version of subroutine pvib for l different values of x.
!     -----------------------------------------------------------------
!
      dimension x(l),c(n,m),p(l,m)

!     dimension q(l,0:n)
      allocatable q(:,:)

      allocate (q(l,0:n))

!
!     primitive sine basis functions
!
      pi = acos(-1.d0)
      shift = pi/(b-a)
      scale = sqrt(2.d0/(b-a))
      if (mode .eq. 0) then
         do k = 1,l
            y = shift*(x(k)-a)
            recur = 2*cos(y)
            q(k,0) = 0.d0
            q(k,1) = scale*sin(y)
            do j = 2,n
               q(k,j) = recur*q(k,j-1)-q(k,j-2)
            enddo
         enddo
      else
         scale = shift*scale
         do k = 1,l
            y = shift*(x(k)-a)
            cosy = cos(y)
            recur = 2*cosy
            q(k,0) = scale
            q(k,1) = scale*cosy
            do j = 2,n
               q(k,j) = recur*q(k,j-1)-q(k,j-2)
            enddo
            do j = 2,n
               q(k,j) = j*q(k,j)
            enddo
         enddo
      endif
!
!     contracted vibrational eigenfunctions
!
      call dgemm ('n','n',l,m,n,1.d0,q(1,1),l,c,n,0.d0,p,l)
      return
      end

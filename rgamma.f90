!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: rgamma.f

      function rgamma(x,odd,even)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     Direct fortran translation of Temme's algol routine for computing
!     rgamma = 1/Gamma(1-x), along with its odd and even parts, for
!     abs(x) .le. 0.5. [ N.M.Temme, J Comput Phys 19 (1975) 324-337 ]
!     -----------------------------------------------------------------
!
      dimension b(12)
      data b / -0.283876542276024d0, -0.076852840844786d0, &
               +0.001706305071096d0, +0.001271927136655d0, &
               +0.000076309597586d0, -0.000004971736704d0, &
               -0.000000865920800d0, -0.000000033126120d0, &
               +0.000000001745136d0, +0.000000000242310d0, &
               +0.000000000009161d0, -0.000000000000170d0 /
      save b
!
      x2 = x*x*8.d0
      alfa = -0.000000000000001d0
      beta = 0.d0
      do i = 12,2,-2
         beta = -(2*alfa+beta)
         alfa = -beta*x2-alfa+b(i)
      enddo
      even = (beta/2.d0+alfa)*x2-alfa+0.921870293650453d0
      alfa = -0.000000000000034d0
      beta = 0.d0
      do i = 11,1,-2
         beta = -(2*alfa+beta)
         alfa = -beta*x2-alfa+b(i)
      enddo
      odd = 2*(alfa+beta)
      rgamma = odd*x+even
      return
      end

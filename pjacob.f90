!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: pjacob.f

      function pjacob (n,a,b,x)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     Jacobi polynomial p(n,a,b;x)
!     Abramowitz and Stegun eq. (22.7.1)
!     -----------------------------------------------------------------
!
      parameter (zero = 0.0d0)
      parameter (half = 0.5d0)
      parameter (one  = 1.0d0)
      parameter (two  = 2.0d0)
!
      if (n .eq. 0) then
        fp = one
      else
        f = one
        apa = a+a
        apb = a+b
        amb = a-b
        apbamb = apb*amb
        apbp1 = apb+one
        apbp2 = apb+two
        onek = zero
        twok = zero
        fp = half*(amb+apbp2*x)
        do k = 1,n-1
          onek = onek+one
          twok = twok+two
          a1 = (twok+two)*(onek+apbp1)*(twok+apb)
          a2 = (twok+apbp1)*apbamb
          a3 = (twok+apb)*(twok+apbp1)*(twok+apbp2)
          a4 = (twok+apa)*(onek+b)*(twok+apbp2)
          fm = f
          f = fp
          fp = ((a2+a3*x)*f-a4*fm)/a1
        enddo
      endif
      pjacob = fp
      return
      end

!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: rotmel.f

      function rotmel (jtot,ipar,k,m,beta)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This function uses redrot to calculate a parity adapted reduced
!     rotation matrix element with angular momentum jtot, parity ipar,
!     angular momentum projections k and m, and argument beta.
!     -----------------------------------------------------------------
!
      rotmel = 0.d0
      if (k.gt.jtot .or. k.lt.0) return
      if (m.gt.jtot .or. m.lt.0) return
      sign = ipar*(-1)**(jtot+k)
      rj = jtot
      rk = k
      rm = m
      rotmel = redrot(rj,rk,rm,beta)+sign*redrot(rj,-rk,rm,beta)
      if (k .eq. 0) rotmel = sqrt(0.5d0)*rotmel
      if (m .eq. 0) rotmel = sqrt(0.5d0)*rotmel
      return
      end

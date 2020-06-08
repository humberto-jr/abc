!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: redrot.f

      function redrot (rj,rk,rm,beta)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This function uses eq. (4.1.23) of Edmonds
!     to calculate the reduced rotation matrix element
!     d(j,k,m;beta) = <jk|exp(+i*beta*Jy/hbar)|jm>.
!     -----------------------------------------------------------------
!
      parameter (zero = 0.0d0)
      parameter (half = 0.5d0)
      parameter (one  = 1.0d0)
      parameter (two  = 2.0d0)
!
!     half integer angular momenta
!
      sj = half*nint(two*rj)
      sk = half*nint(two*rk)
      sm = half*nint(two*rm)
!
!     projection ranges
!
      redrot = zero
      if (sk.gt.sj .or. sk.lt.-sj)  return
      if (sm.gt.sj .or. sm.lt.-sj)  return
      if (mod(sj-sk,one) .ne. zero) return
      if (mod(sj-sm,one) .ne. zero) return
!
!     reflection symmetries
!
      if (sk+sm .ge. zero) then
        if (sk-sm .ge. zero) then
          tk = sk
          tm = sm
          isign = 0
        else
          tk = sm
          tm = sk
          isign = sk-sm
        endif
      else
        if (sk-sm .ge. zero) then
          tk = -sm
          tm = -sk
          isign = 0
        else
          tk = -sk
          tm = -sm
          isign = sk-sm
        endif
      endif
!
!     evaluation
!
      n = sj-tk
      ia = tk-tm
      ib = tk+tm
      a = ia
      b = ib
      beta2 = half*beta
      cosb2 = cos(beta2)
      sinb2 = sin(beta2)
      cosb = (cosb2-sinb2)*(cosb2+sinb2)
      d1 = pjacob(n,a,b,cosb)
      d2 = cosb2**ib*sinb2**ia
      d3 = d1*d2
      d4 = d3*d3
      ti = tm
      do i = 1,ia
         ti = ti+one
         d4 = d4*(sj+ti)/(sj-ti+one)
      enddo
      d4 = sqrt(d4)
      redrot = sign(d4,d3)
      if (mod(isign,2) .ne. 0) redrot = -redrot
      return
      end

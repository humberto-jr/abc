!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: rbessk.f

      subroutine rbessk (ell,x,ck,dk,ek)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     Modified Riccati-Bessel function of the third kind
!     and its first derivative with respect to x:
!
!     k(ell,x) = ck BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LSTHPOT.
!     d/dx k(ell,x) = dk BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LST
!     -----------------------------------------------------------------
!
      if (x.le.0.0d0 .or. ell.lt.-0.5d0) stop 'rbessk 0'
      v = ell+0.5d0
      call mbessk (v,x,ck,dk,ek)
      pi = acos(-1.d0)
      ex = 0.5d0*dlog(pi*x/2.d0)
      dk = dk+ck/(2.d0*x)
      sk = sqrt(ck*ck+dk*dk)
      ck = ck/sk
      dk = dk/sk
      ek = ek+dlog(sk)+ex
      return
      end

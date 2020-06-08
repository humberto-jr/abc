!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: parset.f

      subroutine parset
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine uses the input parameters emax and jmax
!     to determine rmin, smax, mvi, mro, and nvi.
!     -----------------------------------------------------------------
!
!     common blocks
!
      double precision mass,mtot,mred
      common /arrays/ mro,mvi,nvi,n
      common /masses/ mass(3),mtot,mred
      common /inputs/ emax,mtr
      common /ranges/ rmin,rmax,smax
      common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
      common /scales/ scale(3),rmlmda
!
!     local arrays
!
      parameter (m = 1000)

!     dimension v(m)
      allocatable v(:)

!     dimension smin(3),vmid(3)
      allocatable smin(:),vmid(:)

      allocate (v(m))
      allocate (smin(3),vmid(3))

!
!     smax
!
      ds = 0.01d0
      smax = 0.d0
      scut = 5.d0
      ered = rmlmda*emax
      do ia = 1,3
         lmid = 1
         do l = 1,m
            sa = l*ds
            call potenl (100.d0,sa,0.d0,v(l),ia)
            if (v(l) .lt. v(lmid)) lmid = l
         enddo
         vmid(ia) = v(lmid)
         ll = 0
         do l = 1,lmid
            if (v(l) .lt. ered) go to 1
            ll = l
         enddo
   1     sum = 0.d0
         lmin = 0
         do l = ll,1,-1
            lmin = l
            sum = sum+ds*sqrt(v(l)-ered)
            if (sum .gt. scut) go to 2
         enddo
   2     smin(ia) = lmin*ds
         lr = m+1
         do l = m,lmid,-1
            if (v(l) .lt. ered) go to 3
            lr = l
         enddo
   3     sum = 0.0d0
         lmax = m+1
         do l = lr,m
            lmax = l
            sum = sum+ds*sqrt(v(l)-ered)
            if (sum .gt. scut) go to 4
         enddo
   4     smax = max(smax,lmax*ds)
      enddo
!
!     rmin
!
      pi = acos(-1.d0)
      rm2 = 0.d0
      rm2 = rm2+(1.d0-mass(1)/mtot)*smin(1)**2
      rm2 = rm2+(1.d0-mass(2)/mtot)*smin(2)**2
      rm2 = rm2+(1.d0-mass(3)/mtot)*smin(3)**2
      rmin = sqrt(rm2)
!
!     nvi
!
      nvi = 0
      range = 2.d0*smax/pi
      do ia = 1,3
         if (ered .gt. vmid(ia)) then
            wave = sqrt(ered-vmid(ia))
            nvia = range*wave
            nvi = max(nvi,nvia)
         endif
      enddo
!
!     mro and mvi
!
      mro = 3*(2*(jmax+1)+(nvi+1))/4
      mvi = 3*(2*(nvi+1)+(jmax+1))/4
!
!     parameter summary
!
      write (6,61) rmin,mro,mvi,nvi,smax
      if (rmin .ge. rmax) stop 'parset 1'
      return
  61  format(/1x,'PARSET:'/1x,70('-')/1x, &
      '  rmin = ',f32.2/1x, &
      '   mro = ',i32/1x, &
      '   mvi = ',i32/1x, &
      '   nvi = ',i32/1x, &
      '  smax = ',f32.2/1x,70('-'))
      end

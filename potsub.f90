!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: potsub.f

      subroutine potsub (r,vev)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine chooses which potential to use on the basis
!     of the atomic masses in common /masses/, and also ensures
!     that the potential is called with the bond lengths in the
!     correct order.
!     -----------------------------------------------------------------
!
      dimension r(3),s(3),m(3)
      double precision mass,mtot,mred
      common /masses/ mass(3),mtot,mred
!
      imax = 1
      imin = 1
      do i = 1,3
         m(i) = nint(mass(i))
         if (m(i) .gt. m(imax)) then
            imax = i
         else if (m(i) .lt. m(imin)) then
            imin = i
         endif
      enddo
      if (imax .eq. imin) then
         imin = 1
         imid = 2
         imax = 3
      else
         imid = 6-imax-imin
      endif
      s(1) = r(imax)
      s(2) = r(imid)
      s(3) = r(imin)
      mmax = m(imax)
      mmid = m(imid)
      mmin = m(imin)
!     if (mmin .lt. 1) stop 'potsub 1'
!     if (mmid .gt. 2) stop 'potsub 2'
!     if (mmax .le. 2) then
!        call hh2pot (s,vev)
!     else if (mmax .eq. 19) then
!        call fh2pot (s,vev)
!     else if (mmax .eq. 35) then
!        call clh2pt (s,vev)
!     else if (mmax .eq. 37) then
!        call clh2pt (s,vev)
!     else
!        stop 'potsub 3'
!     endif
      call pes(s(1),s(2),s(3),vev)
      return
      end

!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: putrec.f

      subroutine putrec (a,n,nrg)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine writes a matrix to the scratch file on unit 10.
!     -----------------------------------------------------------------
!
      dimension a(n,n)
      irec = nrg+1
      write (unit=10,rec=irec) a
      return
      end

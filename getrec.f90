!     File created at Fri Jun  5 21:58:56 PDT 2020
!     Original source code: getrec.f

      subroutine getrec (a,n,nrg)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine reads a matrix from the scratch file on unit 10.
!     -----------------------------------------------------------------
!
      dimension a(n,n)
      irec = nrg+1
      read (unit=10,rec=irec) a
      return
      end

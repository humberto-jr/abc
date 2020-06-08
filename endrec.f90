!     File created at Fri Jun  5 21:58:56 PDT 2020
!     Original source code: endrec.f

      subroutine endrec
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine closes and deletes the scratch file on unit 10.
!     -----------------------------------------------------------------
!
      close (unit=10,status='delete')
      return
      end

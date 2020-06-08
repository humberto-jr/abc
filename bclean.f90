!     File created at Fri Jun  5 21:58:55 PDT 2020
!     Original source code: bclean.f

      subroutine bclean (cvi,nvi,n,smin,smax,mvi)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine removes the phase ambiguity from the vibrational
!     expansion coefficient array cvi by ensuring that each vibrational
!     eigenfunction is positive at its inner turning point.
!     -----------------------------------------------------------------
!
      dimension cvi(nvi,n)

!     dimension pvi(n,mvi)
      allocatable pvi(:,:)

      allocate (pvi(n,mvi))

!
      ds = (smax-smin)/mvi
      do k = 1,mvi
         sa = smin+(k-0.5d0)*ds
         call pvib (smin,sa,smax,cvi,nvi,n,pvi(1,k),0)
      enddo
      do 1 j = 1,n
         pmax = 0.d0
         do k = 1,mvi
            pmax = max(pmax,abs(pvi(j,k)))
         enddo
         ptest = 0.1d0*pmax
         do k = 1,mvi
            if (abs(pvi(j,k)) .gt. ptest) then
               if (pvi(j,k) .lt. 0.d0) then
                  do i = 1,nvi
                     cvi(i,j) = -cvi(i,j)
                  enddo
               endif
               go to 1
            endif
         enddo
   1  continue
      return
      end

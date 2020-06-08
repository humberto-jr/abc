!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: overlp.f

      subroutine overlp (rhoa,cvia,rhob,cvib,ilev,jlev,klev,s)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine constructs the sector-to-sector overlap matrix
!     s between sectors centred at hyperradii rhoa and rhob.
!     Compare with subroutine couple.
!     -----------------------------------------------------------------
!
!     common blocks
!
      common /arrays/ mro,mvi,nvi,n
      common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
!
!     input arrays
!
      dimension cvia(nvi,n),cvib(nvi,n)
      dimension ilev(n),jlev(n),klev(n)
      dimension s(n,n)
!
      do j = 1,n
         do i = 1,n
            s(i,j) = 0.d0
         enddo
      enddo
      na = 3-iabs(jpar)
      do ia = 1,na
         call overdi (rhoa,cvia,rhob,cvib,ilev,jlev,klev,s,ia)
         do ib = 1,na
            if (ib .ne. ia) then
             call overex (rhoa,cvia,rhob,cvib,ilev,jlev,klev,s,ia,ib)
            endif
         enddo
      enddo
      if (na .eq. 2) then
         ia = 2
         ib = 3
         call overex (rhoa,cvia,rhob,cvib,ilev,jlev,klev,s,ia,ib)
      endif
      return
      end

!     File created at Fri Jun  5 21:58:55 PDT 2020
!     Original source code: arrang.f

      subroutine arrang (ilev,n,jpar,ia,nla,nha,na)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine calculates arrangement channel indices nla, nha
!     and na such that the na = nha-nla+1 channels in arrangement ia
!     occur between nla and nha (assuming that the channels have been
!     sorted in a,v,j order by subroutine basort).
!     -----------------------------------------------------------------
!
      dimension ilev(n)
      dimension nlo(3),nhi(3)
!
      do i = 1,3
         nlo(i) = n+1
         nhi(i) = 0
      enddo
      do j = 1,n
         nlo(ilev(j)) = min(nlo(ilev(j)),j)
         nhi(ilev(j)) = max(nhi(ilev(j)),j)
      enddo
      if (jpar .ne. 0) then
         nlo(3) = nlo(2)
         nhi(3) = nhi(2)
      endif
      nla = nlo(ia)
      nha = nhi(ia)
      na = nha-nla+1
      return
      end

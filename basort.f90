!     File created at Fri Jun  5 21:58:55 PDT 2020
!     Original source code: basort.f

      subroutine basort (cvi,nvi,cro,eint,ilev,jlev,klev,nlev,n)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine sorts the channels in the input
!     basis set in order of increasing a,v,j,k.
!     -----------------------------------------------------------------
!
      dimension cvi(nvi,n),cro(3,n)
      dimension eint(n),ilev(n),jlev(n),klev(n),nlev(n)
!
      do j = 1,n-1
         k = j
         do i = j+1,n
            if (ilev(i) .lt. ilev(k)) k = i
            if (ilev(i) .gt. ilev(k)) go to 1
            if (nlev(i) .lt. nlev(k)) k = i
            if (nlev(i) .gt. nlev(k)) go to 1
            if (jlev(i) .lt. jlev(k)) k = i
            if (jlev(i) .gt. jlev(k)) go to 1
            if (klev(i) .lt. klev(k)) k = i
   1        continue
         enddo
         if (k .ne. j) then
            do i = 1,nvi
               cswap = cvi(i,j)
               cvi(i,j) = cvi(i,k)
               cvi(i,k) = cswap
            enddo
            do i = 1,3
               cswap = cro(i,j)
               cro(i,j) = cro(i,k)
               cro(i,k) = cswap
            enddo
            eswap = eint(j)
            eint(j) = eint(k)
            eint(k) = eswap
            iswap = ilev(j)
            ilev(j) = ilev(k)
            ilev(k) = iswap
            jswap = jlev(j)
            jlev(j) = jlev(k)
            jlev(k) = jswap
            kswap = klev(j)
            klev(j) = klev(k)
            klev(k) = kswap
            nswap = nlev(j)
            nlev(j) = nlev(k)
            nlev(k) = nswap
         endif
      enddo
      return
      end

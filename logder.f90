!     File created at Fri Jun  5 21:58:56 PDT 2020
!     Original source code: logder.f

      subroutine logder (drho,ered,nnrg,y,z,t,e,n,np,npp)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     Constant reference potential log derivative propagator,
!     with nnrg energies propagated simultaneously.
!     -----------------------------------------------------------------
!
!     input arrays
!
      dimension ered(nnrg)
      dimension y(n,n),z(n,n),t(n,n),e(n)
!
!     local arrays
!
!     dimension y1(n),y2(n),y3(n),y4(n)
      allocatable y1(:),y2(:),y3(:),y4(:)

      allocate (y1(n),y2(n),y3(n),y4(n))

!
!     loop over energies
!
      do nrg = 1,nnrg
!
!        constant reference potential propagators
!
         if (npp .le. n) then
            h = drho
            eps = 1.d-16/(h*h)
            do i = 1,np
               psq = e(i)-ered(nrg)
               if (psq .gt. eps) then
                  p = sqrt(psq)
                  y1(i) = p/tanh(p*h)
                  y2(i) = p/sinh(p*h)
               else if (psq .lt. -eps) then
                  p = sqrt(-psq)
                  y1(i) = p/tan(p*h)
                  y2(i) = p/sin(p*h)
               else
                  y1(i) = 1.d0/h
                  y2(i) = 1.d0/h
               endif
               y3(i) = y2(i)
               y4(i) = y1(i)
            enddo
         endif
!
!        log derivative propagation
!
         if (npp .le. 0) then
            do j = 1,np
               do i = 1,np
                  y(i,j) = 0.d0
               enddo
               y(j,j) = y4(j)
            enddo
         else if (npp .le. n) then
            call getrec (y,n,nrg)
            call dgemm ('n','n',npp,np,npp,1.d0,y,n,t,n,0.d0,z,n)
            call dgemm ('t','n',np, np,npp,1.d0,t,n,z,n,0.d0,y,n)
            do j = 1,np
               do i = 1,j
                  y(i,j) = 0.5d0*(y(i,j)+y(j,i))
               enddo
               y(j,j) = y(j,j)+y1(j)
            enddo
            call syminv (y,n,np,ierr)
            if (ierr .ne. 0) stop 'logder 1'
            do j = 1,np
               do i = 1,j
                  y(i,j) = -y3(i)*y(i,j)*y2(j)
                  y(j,i) = y(i,j)
               enddo
               y(j,j) = y(j,j)+y4(j)
            enddo
         else if (npp .gt. n) then
            call getrec (y,n,nrg)
            call dgemm ('n','t',np,n,np,1.d0,y,n,t,n,0.d0,z,n)
            call dgemm ('n','n',n, n,np,1.d0,t,n,z,n,0.d0,y,n)
            do j = 1,np
               do i = 1,j
                  y(i,j) = 0.5d0*(y(i,j)+y(j,i))
                  y(j,i) = y(i,j)
               enddo
            enddo
         endif
         call putrec (y,n,nrg)
      enddo
      return
      end

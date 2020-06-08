!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: rbesjy.f

      subroutine rbesjy (ell,x,cj,dj,ej,cy,dy,ey)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     Riccati-Bessel functions of fractional order
!     and their first derivatives with respect to x:
!
!     j(ell,x) = cj BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LSTHPOT.
!     y(ell,x) = cy BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LSTHPOT.
!     d/dx j(ell,x) = dj BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LST
!     d/dx y(ell,x) = dy BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LST
!     -----------------------------------------------------------------
!
      if (x.le.0.0d0 .or. ell.lt.-0.5d0) stop 'rbesjy 0'
      v = ell+0.5d0
      call bessjy (v,x,cj,dj,ej,cy,dy,ey)
      pi = acos(-1.d0)
      ex = 0.5d0*dlog(pi*x/2.d0)
      dj = dj+cj/(2.d0*x)
      sj = sqrt(cj*cj+dj*dj)
      cj = cj/sj
      dj = dj/sj
      ej = ej+dlog(sj)+ex
      dy = dy+cy/(2.d0*x)
      sy = sqrt(cy*cy+dy*dy)
      cy = cy/sy
      dy = dy/sy
      ey = ey+dlog(sy)+ex
      return
      end

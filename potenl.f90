!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: potenl.f

      subroutine potenl (ra,sa,cosa,va,ia)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine returns the potential energy surface va as a
!     function of mass-scaled Jacobi coordinates in arrangement ia.
!     -----------------------------------------------------------------
!
      dimension r(3)
      double precision mass,mtot,mred
      common /masses/ mass(3),mtot,mred
      common /scales/ scale(3),rmlmda
!
      ib = ia+1
      if (ib .gt. 3) ib = 1
      ic = 6-ia-ib
      rap = ra/scale(ia)
      cap = 2.0d0*cosa*rap
      sap = scale(ia)*sa
      sbp = mass(ib)/(mass(ib)+mass(ic))*sap
      scp = mass(ic)/(mass(ib)+mass(ic))*sap
      r(ia) = sap
      r(ib) = sqrt(rap**2-cap*sbp+sbp**2)
      r(ic) = sqrt(rap**2+cap*scp+scp**2)
      call potsub (r,vev)
      va = rmlmda*vev
      return
      end

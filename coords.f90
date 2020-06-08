!     File created at Fri Jun  5 21:58:55 PDT 2020
!     Original source code: coords.f

      subroutine coords (ta,cosa,ia,tb,cosb,ib,xab)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine calculates ib-arrangement from ia-arrangement
!     Delves hyperangular coordinates, along with the angle xab
!     between the body-frame z-axes of the two arrangements.
!     -----------------------------------------------------------------
!
      double precision mass,mtot,mred
      common /masses/ mass(3),mtot,mred
      common /scales/ scale(3),rmlmda
!
      if (ib .eq. ia) then
         tb = ta
         cosb = cosa
         xab = 0.d0
      else
         sab = 1.d0/(scale(ia)*scale(ib))
         cab = -mred*sab/mass(6-ia-ib)
         if (ia-ib.eq.1 .or. ib-ia.eq.2) sab = -sab
         ra = cos(ta)
         sa = sin(ta)
         rsca = ra*sa*cosa
         rb = sqrt((cab*ra)**2-2.d0*cab*sab*rsca+(sab*sa)**2)
         sb = sqrt((sab*ra)**2+2.d0*cab*sab*rsca+(cab*sa)**2)
         rscb = cab*sab*(ra-sa)*(ra+sa)+rsca*(cab-sab)*(cab+sab)
         cosb = rscb/(rb*sb)
         cosb = max(cosb,-1.d0)
         cosb = min(cosb,+1.d0)
         tb = atan(sb/rb)
         cxab = (rb*cab+sb*cosb*sab)/ra
         cxab = max(cxab,-1.d0)
         cxab = min(cxab,+1.d0)
         xab = acos(cxab)
         if (ia-ib.eq.1 .or. ib-ia.eq.2) xab = -xab
      endif
      return
      end

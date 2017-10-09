c-------------------------------------------------------------
c    A routine to find primitive variables when the EOS is
c    given by a table.
c-------------------------------------------------------------
      subroutine table_inversion(retval,rholoc,ploc,uloc,tloc,
     &     vxloc,vyloc,vzloc,csloc,dloc,eloc,sxloc,syloc,szloc,
     &     Ye,rho0,gamxx,gamxy,gamxz,gamyy,gamyz,gamzz,
     &     gamup_xx,gamup_xy,gamup_xz,gamup_yy,
     &     gamup_yz,gamup_zz,lapse,detgam,shiftx,shifty,shiftz)
      implicit none
      integer retval
      real*8 rholoc, ploc, uloc, tloc, vxloc, vyloc, vzloc
      real*8 dloc, eloc, sxloc, syloc, szloc
      real*8 Ye, rho0, csloc
      real*8 gamxx,gamxy,gamxz,gamyy,gamyz,gamzz
      real*8 gamup_xx, gamup_xy, gamup_xz
      real*8 gamup_yy, gamup_yz, gamup_zz
      real*8 shiftx, shifty, shiftz
      real*8 lapse, detgam, sg, betas
      real*8 gupxs, gupys, gupzs
      real*8 supx, supy, supz
      real*8 utloc
      real*8 rhohloc   ! density times specific enthalpy.
      real*8 c1, c2
      real*8 w, eps, wbar    ! our independent variables.
      real*8 x(2)
      integer NTRIAL
      real*8 TOLF, TOLX
      parameter (NTRIAL = 20)
      parameter (TOLF = 1.d-10, TOLX = 1.d-10)

      sg = sqrt(detgam)

      ! Compute rhs quantities.
      c1 = gamup_xx * sxloc**2 + 2.d0*gamup_xy * sxloc*syloc +
     &     2.d0*gamup_xz * sxloc*szloc + gamup_yy * syloc**2 + 
     &     2.d0*gamup_yz * syloc*szloc + gamup_zz * szloc**2
      c1 = c1 / detgam

      betas = shiftx*sxloc + shifty*syloc + shiftz*szloc
      c2 = eloc - betas
      c2 = -c2 / (lapse * sqrt(detgam))

      ! Calculate initial guesses.
      w   = dloc/(sg*rholoc)
      eps = uloc/rholoc

      ! Make sure things are within range.
      if (w**2 < 1.d0) then
         w = 1.d0
         wbar = 0.d0
      end if
      wbar = sqrt(w**2-1.d0)

      if (wbar == 0.d0 .and. ((sxloc .ne. 0.d0) .or. 
     &     (syloc .ne. 0.d0) .or. (szloc .ne. 0.d0))) then

         ! We need to do something about the wbar estimate so 
         ! that it doesn't vanish.  Assume small velocities.
         rhohloc = rholoc + ploc + uloc
         vxloc = sxloc / (lapse*sg*rhohloc)
         vyloc = sxloc / (lapse*sg*rhohloc)
         vzloc = sxloc / (lapse*sg*rhohloc)

         utloc = lapse**2 - gamxx*(shiftx+vxloc)**2 -
     &        2.d0*gamxy*(shiftx+vxloc)*(shifty+vyloc) -
     &        2.d0*gamxz*(shiftx+vxloc)*(shiftz+vzloc) - 
     &        gamyy*(shifty+vyloc)**2 -
     &        2.d0*gamyz*(shifty+vyloc)*(shiftz+vzloc) -
     &        gamzz*(shiftz+vzloc)**2 

         if (utloc < 0.d0) then
            write(*,*) 'your inversion fix went superluminal!'
            stop
         end if

         utloc = 1.d0/sqrt(utloc)

         w = lapse*utloc
         wbar = sqrt(w**2-1.d0)

      end if

      x(1) = wbar
      x(2) = eps
      
      call mnewt(NTRIAL,x,2,TOLX,TOLF,retval,c1,c2,dloc,sg,Ye,tloc,rho0)

      if (retval > 0) then
         write(*,*) 'table inversion:  mnewt failed to converge.'
         stop
      end if

      wbar = x(1)
      eps  = x(2)

      ! Calculate everybody else.
      w = sqrt(1.d0+wbar**2)
      rholoc = dloc / (sg * w)
      uloc = eps * rholoc
      call eos_wrapper(ploc,rholoc,eps,Ye,tloc,rho0,csloc)
      rhohloc = rholoc + ploc + uloc

      gupxs = gamup_xx*sxloc + gamup_xy*syloc + gamup_xz*szloc
      gupys = gamup_xy*sxloc + gamup_yy*syloc + gamup_yz*szloc
      gupzs = gamup_xz*sxloc + gamup_yz*syloc + gamup_zz*szloc

      supx = shiftx*eloc/lapse**2 + gupxs - shiftx*betas/lapse**2
      supy = shifty*eloc/lapse**2 + gupys - shifty*betas/lapse**2
      supz = shiftz*eloc/lapse**2 + gupzs - shiftz*betas/lapse**2

      vxloc = (lapse*gupxs/sg - shiftx*ploc)/(rhohloc * w**2)
      vyloc = (lapse*gupys/sg - shifty*ploc)/(rhohloc * w**2)
      vzloc = (lapse*gupzs/sg - shiftz*ploc)/(rhohloc * w**2)
      end subroutine

c     A version which carries information about the local position
      subroutine table_inversion_loc(retval,rholoc,ploc,uloc,tloc,
     &     vxloc,vyloc,vzloc,csloc,dloc,eloc,sxloc,syloc,szloc,
     &     Ye,rho0,gamxx,gamxy,gamxz,gamyy,gamyz,gamzz,
     &     gamup_xx,gamup_xy,gamup_xz,gamup_yy,
     &     gamup_yz,gamup_zz,lapse,detgam,shiftx,shifty,shiftz,
     &     il,jl,kl)
      implicit none
      integer retval,il,jl,kl
      real*8 rholoc, ploc, uloc, tloc, vxloc, vyloc, vzloc
      real*8 dloc, eloc, sxloc, syloc, szloc
      real*8 Ye, rho0, csloc
      real*8 gamxx,gamxy,gamxz,gamyy,gamyz,gamzz
      real*8 gamup_xx, gamup_xy, gamup_xz
      real*8 gamup_yy, gamup_yz, gamup_zz
      real*8 shiftx, shifty, shiftz
      real*8 lapse, detgam, sg, betas
      real*8 gupxs, gupys, gupzs
      real*8 supx, supy, supz
      real*8 utloc
      real*8 rhohloc   ! density times specific enthalpy.
      real*8 c1, c2
      real*8 w, eps, wbar    ! our independent variables.
      real*8 x(2)
      integer NTRIAL
      real*8 TOLF, TOLX
      parameter (NTRIAL = 20)
      parameter (TOLF = 1.d-10, TOLX = 1.d-10)

      sg = sqrt(detgam)

      ! Compute rhs quantities.
      c1 = gamup_xx * sxloc**2 + 2.d0*gamup_xy * sxloc*syloc +
     &     2.d0*gamup_xz * sxloc*szloc + gamup_yy * syloc**2 + 
     &     2.d0*gamup_yz * syloc*szloc + gamup_zz * szloc**2
      c1 = c1 / detgam

      betas = shiftx*sxloc + shifty*syloc + shiftz*szloc
      c2 = eloc - betas
      c2 = -c2 / (lapse * sqrt(detgam))

      ! Calculate initial guesses.
      w   = dloc/(sg*rholoc)
      eps = uloc/rholoc

      ! Make sure things are within range.
      if (w**2 < 1.d0) then
         w = 1.d0
         wbar = 0.d0
      end if
      wbar = sqrt(w**2-1.d0)

      if (wbar == 0.d0 .and. ((sxloc .ne. 0.d0) .or. 
     &     (syloc .ne. 0.d0) .or. (szloc .ne. 0.d0))) then

         ! We need to do something about the wbar estimate so 
         ! that it doesn't vanish.  Assume small velocities.
         rhohloc = rholoc + ploc + uloc
         vxloc = sxloc / (lapse*sg*rhohloc)
         vyloc = sxloc / (lapse*sg*rhohloc)
         vzloc = sxloc / (lapse*sg*rhohloc)

         utloc = lapse**2 - gamxx*(shiftx+vxloc)**2 -
     &        2.d0*gamxy*(shiftx+vxloc)*(shifty+vyloc) -
     &        2.d0*gamxz*(shiftx+vxloc)*(shiftz+vzloc) - 
     &        gamyy*(shifty+vyloc)**2 -
     &        2.d0*gamyz*(shifty+vyloc)*(shiftz+vzloc) -
     &        gamzz*(shiftz+vzloc)**2 

         if (utloc < 0.d0) then
            write(*,*) 'your inversion fix went superluminal!'
            stop
         end if

         utloc = 1.d0/sqrt(utloc)

         w = lapse*utloc
         wbar = sqrt(w**2-1.d0)

      end if

      x(1) = wbar
      x(2) = eps
      
      if (wbar > 0) then
         call mnewt_loc(NTRIAL,x,2,TOLX,TOLF,retval,c1,c2,dloc,sg,Ye,
     &        tloc,rho0,il,jl,kl)
         
         if (retval > 0) then
            write(*,*) 'table inversion:  mnewt failed to converge.'
            stop
         end if
         
         wbar = x(1)
         eps  = x(2)
      else
         ! Invert analytically assuming v^i = 0.
         rholoc = dloc / sg
         eps = c2 / rholoc - 1.d0
      end if

      ! Calculate everybody else.
      w = sqrt(1.d0+wbar**2)
      rholoc = dloc / (sg * w)
      uloc = eps * rholoc
c      call eos_wrapper(ploc,rholoc,eps,Ye,tloc,rho0,csloc)
      call eos_wrapper_loc(ploc,rholoc,eps,Ye,tloc,rho0,csloc,il,jl,kl)
      rhohloc = rholoc + ploc + uloc

      gupxs = gamup_xx*sxloc + gamup_xy*syloc + gamup_xz*szloc
      gupys = gamup_xy*sxloc + gamup_yy*syloc + gamup_yz*szloc
      gupzs = gamup_xz*sxloc + gamup_yz*syloc + gamup_zz*szloc

      supx = shiftx*eloc/lapse**2 + gupxs - shiftx*betas/lapse**2
      supy = shifty*eloc/lapse**2 + gupys - shifty*betas/lapse**2
      supz = shiftz*eloc/lapse**2 + gupzs - shiftz*betas/lapse**2

      vxloc = (lapse*gupxs/sg - shiftx*ploc)/(rhohloc * w**2)
      vyloc = (lapse*gupys/sg - shifty*ploc)/(rhohloc * w**2)
      vzloc = (lapse*gupzs/sg - shiftz*ploc)/(rhohloc * w**2)
      end subroutine

c--------------------------------------------------------------------
c     User-supplied subroutine required by numerical recipes
c     Newton-Raphson method.
c
c     fvec(1) = (rho h)^2 wbar^2 (wbar^2+1) - c1 = 0
c     fvec(2) = (rho h) (wbar^2+1) - P - c2      = 0
c
c     Independent variables are wbar (where w = \alpha u^t and 
c     wbar^2 = w^2-1), and \epsilon, which 
c     is the specific internal energy
c--------------------------------------------------------------------
      subroutine usrfun(x,n,NP,fvec,fjac,c1,c2,dloc,sg,Ye,temp,rho0)
      integer n, NP
      real*8 x(n),fjac(NP,NP),fvec(NP)
      real*8 c1, c2
      real*8 dloc,sg,Ye,temp,rho0

      ! Internals 
      real*8 rho, rhoh, P, w, eps, wbar
      real*8 dpdw, dpdeps, dw, deps
      real*8 wp,wm,epsp,epsm
      real*8 rhop,rhom
      real*8 P_wp,P_wm,P_epsp,P_epsm
      real*8 t, t_wp, t_wm, t_epsp, t_epsm
      real*8 cs_temp
      
      wbar = x(1)
      eps  = x(2)

      w = sqrt(1.d0+wbar**2)
      dw = w * 0.01
      deps = eps * 0.01
      
      wp   = w   + dw/2.d0
      wm   = w   - dw/2.d0
      epsp = eps + deps/2.d0
      epsm = eps - deps/2.d0

      rho  = dloc / (sg * w)
      rhop = dloc / (sg * wp)
      rhom = dloc / (sg * wm)

      ! Make five calls to the eos for our difference star.
      ! Maybe you better check first to make sure that all
      ! of your values are on the table.
      ! Or put this in the wrapper?
      ! We will not be using the temperatures of sounds speeds
      ! that come from these calls.

      t      = temp
      t_wp   = temp
      t_wm   = temp
      t_epsp = temp
      t_epsm = temp

      call eos_wrapper(P,rho,eps,Ye,t,rho0,cs_temp)
      call eos_wrapper(P_wp,rhop,eps,Ye,t_wp,rho0,cs_temp)
      call eos_wrapper(P_wm,rhom,eps,Ye,t_wm,rho0,cs_temp)
      call eos_wrapper(P_epsp,rho,epsp,Ye,t_epsp,rho0,cs_temp)
      call eos_wrapper(P_epsm,rho,epsm,Ye,t_epsm,rho0,cs_temp)

      dpdw   = (P_wp   - P_wm)  /dw
      dpdeps = (P_epsp - P_epsm)/deps
      rhoh = rho + P + rho*eps

      fvec(1) = (rhoh*wbar)**2*(wbar**2+1.d0) - c1
      fvec(2) = rhoh*(wbar**2+1.d0) - P - c2
      
      fjac(1,1) = 2.d0*rhoh*wbar**3*w*(-dloc*(eps+1.d0)/(sg*w**2)+dpdw)+ 
     &     2.d0*rhoh**2*wbar*(1.d0+2*wbar**2)
      fjac(1,2) = 2.d0*rhoh*wbar**2*(dpdeps*w**2+dloc/sg*w)

      fjac(2,1) = -(wbar/w)*dloc/sg*(1.d0+eps) + dpdw*wbar**3/w + 
     &     2.d0*wbar*rhoh
      fjac(2,2) = dpdeps*wbar**2 + dloc/sg*w

      end subroutine

c     Branson messing around.  A version with location information.
      subroutine usrfun_loc(x,n,NP,fvec,fjac,c1,c2,dloc,sg,Ye,temp,rho0,
     &    il,jl,kl)
      integer n, NP
      real*8 x(n),fjac(NP,NP),fvec(NP)
      real*8 c1, c2
      real*8 dloc,sg,Ye,temp,rho0
      integer il,jl,kl

      ! Internals 
      real*8 rho, rhoh, P, w, eps, wbar
      real*8 dpdw, dpdeps, dw, deps
      real*8 wp,wm,epsp,epsm
      real*8 rhop,rhom
      real*8 P_wp,P_wm,P_epsp,P_epsm
      real*8 t, t_wp, t_wm, t_epsp, t_epsm
      real*8 cs_temp
      
      wbar = x(1)
      eps  = x(2)

      w = sqrt(1.d0+wbar**2)
      dw = w * 0.01
      deps = eps * 0.01
      
      wp   = w   + dw/2.d0
      wm   = w   - dw/2.d0
      epsp = eps + deps/2.d0
      epsm = eps - deps/2.d0

      rho  = dloc / (sg * w)
      rhop = dloc / (sg * wp)
      rhom = dloc / (sg * wm)

      ! Make five calls to the eos for our difference star.
      ! Maybe you better check first to make sure that all
      ! of your values are on the table.
      ! Or put this in the wrapper?
      ! We will not be using the temperatures of sounds speeds
      ! that come from these calls.

      t      = temp
      t_wp   = temp
      t_wm   = temp
      t_epsp = temp
      t_epsm = temp

      if (il == 103 .and. jl==1 .and. kl==1) then
         write(*,*) 'before p'
         write(*,*) 'rho  = ', rho
         write(*,*) 'eps  = ', eps
         write(*,*) 'temp = ', t
      end if
      call eos_wrapper_loc(P,rho,eps,Ye,t,rho0,cs_temp,il,jl,kl)
      call eos_wrapper_loc(P_wp,rhop,eps,Ye,t_wp,rho0,cs_temp,il,jl,kl)
      call eos_wrapper_loc(P_wm,rhom,eps,Ye,t_wm,rho0,cs_temp,il,jl,kl)
      call eos_wrapper_loc(P_epsp,rho,epsp,Ye,t_epsp,rho0,cs_temp,il,jl,
     &     kl)
      call eos_wrapper_loc(P_epsm,rho,epsm,Ye,t_epsm,rho0,cs_temp,il,jl,
     &     kl)

      dpdw   = (P_wp   - P_wm)  /dw
      dpdeps = (P_epsp - P_epsm)/deps
      rhoh = rho + P + rho*eps

      fvec(1) = (rhoh*wbar)**2*(wbar**2+1.d0) - c1
      fvec(2) = rhoh*(wbar**2+1.d0) - P - c2
      
      fjac(1,1) = 2.d0*rhoh*wbar**3*w*(-dloc*(eps+1.d0)/(sg*w**2)+dpdw)+ 
     &     2.d0*rhoh**2*wbar*(1.d0+2*wbar**2)
      fjac(1,2) = 2.d0*rhoh*wbar**2*(dpdeps*w**2+dloc/sg*w)

      fjac(2,1) = -(wbar/w)*dloc/sg*(1.d0+eps) + dpdw*wbar**3/w + 
     &     2.d0*wbar*rhoh
      fjac(2,2) = dpdeps*wbar**2 + dloc/sg*w

      end subroutine

      subroutine mnewt(ntrial,x,n,tolx,tolf,retval,c1,c2,dloc,sg,Ye,
     &     temp,rho0)
      integer n, ntrial, NP, retval
      real*8 tolf,tolx,x(n)
      real*8 c1, c2, dloc, sg, Ye, temp, rho0
      parameter (NP=15)
c     USES lubksb, ludcmp, usrfun

      integer i,k,indx(NP)
      real*8 d,errf,errx,fjac(NP,NP),fvec(NP),p(NP)

      retval = 0
      do k=1,ntrial
         
c     Branson messing around.
         do i=1,n
            fvec(i)=0.d0
         end do

         call usrfun(x,n,NP,fvec,fjac,c1,c2,dloc,sg,Ye,temp,rho0)
         errf=0.d0
         
         do i=1,n
            errf=errf+abs(fvec(i))
         enddo

         if (errf .le. tolf) return

         do i=1,n
            p(i)=-fvec(i)
         enddo

c        Solve linear equations using LU decomposition.
         call ludcmp(fjac,n,NP,indx,d)
         call lubksb(fjac,n,NP,indx,p)
         errx = 0.d0

         do i=1,n
            errx=errx+abs(p(i))
            x(i)=x(i)+p(i)
         enddo
         if (errx .le. tolx) return
      end do

      if (k==ntrial) retval = 1
      return

      end 

      subroutine mnewt_loc(ntrial,x,n,tolx,tolf,retval,c1,c2,dloc,sg,Ye,
     &     temp,rho0,il,jl,kl)
      integer n, ntrial, NP, retval
      real*8 tolf,tolx,x(n)
      real*8 c1, c2, dloc, sg, Ye, temp, rho0
      integer il, jl, kl
      parameter (NP=15)
c     USES lubksb, ludcmp, usrfun

      integer i,k,indx(NP)
      real*8 d,errf,errx,fjac(NP,NP),fvec(NP),p(NP)

c      write(*,*) 'inside mnewt_loc at ', il, ' ', jl, ' ', kl

      retval = 0
      do k=1,ntrial
         
c     Branson messing around.
         do i=1,n
            fvec(i)=0.d0
         end do

         call usrfun_loc(x,n,NP,fvec,fjac,c1,c2,dloc,sg,Ye,temp,rho0,il,
     &        jl,kl)
         errf=0.d0
         
         do i=1,n
            errf=errf+abs(fvec(i))
         enddo

         if (errf .le. tolf) return

         do i=1,n
            p(i)=-fvec(i)
         enddo

c        Solve linear equations using LU decomposition.
         call ludcmp_loc(fjac,n,NP,indx,d,il,jl,kl)
         call lubksb(fjac,n,NP,indx,p)
         errx = 0.d0

         do i=1,n
            errx=errx+abs(p(i))
            x(i)=x(i)+p(i)
         enddo

         ! Branson messing around.
         if (x(2) < 0.d0) x(2) = 0.d0

         if (errx .le. tolx) return
      end do

      if (k==ntrial) retval = 1
      return

      end 

      subroutine ludcmp(a,n,np,indx,d)
      integer n,np,indx(n),NMAX
      real*8 d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.d-20)
      integer i,imax,j,k
      real*8 aamax,dum,sum,vv(NMAX)

      d=1.d0
      do i=1,n
         aamax=0.d0
         do j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
         end do
         if (aamax .eq. 0.d0) pause 'singular matrix in ludcmp'
         vv(i)=1.d0/aamax
      end do

      do j=1,n
         do i=1,j-1
            sum=a(i,j)
            do k=1,i-1
               sum = sum - a(i,k)*a(k,j)
            end do
            a(i,j)=sum
         end do
         aamax=0.d0

         do i=j,n
            sum=a(i,j)
            do k=1,j-1
               sum=sum-a(i,k)*a(k,j)
            end do
            a(i,j)=sum
            dum=vv(i)*abs(sum)
            if (dum.ge.aamax) then
               imax=i
               aamax=dum
            end if
         end do

         if (j.ne.imax) then
            do k=1,n
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
            end do
            d=-d
            vv(imax)=vv(j)
         end if

         indx(j)=imax
         
         if (a(j,j) .eq. 0.d0) a(j,j)=TINY

         if (j.ne.n) then
            dum = 1.d0/a(j,j)
            do i=j+1,n
               a(i,j)=a(i,j)*dum
            end do
         end if

      end do

      return
      end 

      subroutine ludcmp_loc(a,n,np,indx,d,il,jl,kl)
      integer n,np,indx(n),NMAX
      real*8 d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.d-20)
      integer i,imax,j,k
      integer il,jl,kl
      real*8 aamax,dum,sum,vv(NMAX)

      d=1.d0
      do i=1,n
         aamax=0.d0
         do j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
         end do
         if (aamax .eq. 0.d0) then
            write(*,*)'singular at', il, ' ', jl, ' ', kl
            stop
         end if
         vv(i)=1.d0/aamax
      end do

      do j=1,n
         do i=1,j-1
            sum=a(i,j)
            do k=1,i-1
               sum = sum - a(i,k)*a(k,j)
            end do
            a(i,j)=sum
         end do
         aamax=0.d0

         do i=j,n
            sum=a(i,j)
            do k=1,j-1
               sum=sum-a(i,k)*a(k,j)
            end do
            a(i,j)=sum
            dum=vv(i)*abs(sum)
            if (dum.ge.aamax) then
               imax=i
               aamax=dum
            end if
         end do

         if (j.ne.imax) then
            do k=1,n
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
            end do
            d=-d
            vv(imax)=vv(j)
         end if

         indx(j)=imax
         
         if (a(j,j) .eq. 0.d0) a(j,j)=TINY

         if (j.ne.n) then
            dum = 1.d0/a(j,j)
            do i=j+1,n
               a(i,j)=a(i,j)*dum
            end do
         end if

      end do

      return
      end 

      subroutine lubksb(a,n,np,indx,b)
      integer n,np,indx(n)
      real*8 a(np,np), b(n)
      integer i,ii,j,ll
      real*8 sum

      ii=0
      do i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if (ii.ne.0) then
            do j=ii,i-1
               sum = sum-a(i,j)*b(j)
            end do
         else if (sum .ne. 0.d0) then
            ii=i
         end if
         b(i)=sum
         
      end do
      do i=n,1,-1
         sum=b(i)
         do j=i+1,n
            sum = sum - a(i,j)*b(j)
         end do
         b(i)=sum/a(i,i)
      end do
      return
      end
            

c-------------------------------------------------------------
c    A routine to find primitive variables when the EOS is
c    given by a table.
c-------------------------------------------------------------
      subroutine table_inversion(retval,rholoc,ploc,uloc,tloc,
     &     vxloc,vyloc,vzloc,csloc,dloc,eloc,sxloc,syloc,szloc,
     &     Ye,rho0,gamxx,gamxy,gamxz,gamyy,gamyz,gamzz,
     &     gamup_xx,gamup_xy,gamup_xz,gamup_yy,
     &     gamup_yz,gamup_zz,lapse,detgam,shiftx,shifty,shiftz,
     &     il,jl,kl,cap,floor,deltat1,deltat2)
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
      real*8 deltat1,deltat2
      logical cap, floor
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
         call mnewt(NTRIAL,x,2,TOLX,TOLF,retval,c1,c2,dloc,sg,Ye,
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
      call eos_wrapper(ploc,rholoc,eps,Ye,tloc,rho0,csloc,il,jl,kl,
     &     cap,floor,deltat1,deltat2)
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
c     for the 1d GR hydro code.  Note that we pass in tau instead
c     of e.
      subroutine table_inversion_1d(retval,rholoc,ploc,uloc,tloc,
     &     vloc,csloc,dloc,tauloc,srloc,
     &     Ye,rho0,psi,lapse,shift,xloc,il,cap,floor,deltat1,deltat2)
      implicit none
      integer retval,il,jl,kl
      real*8 rholoc, ploc, uloc, tloc, vloc
      real*8 dloc, tauloc, srloc
      real*8 Ye, rho0, csloc
      real*8 shift, psi, psi4
      real*8 lapse, sg 
      real*8 supx, supy, supz
      real*8 utloc, xloc
      real*8 rhohloc   ! density times specific enthalpy.
      real*8 c1, c2
      real*8 w, eps, wbar    ! our independent variables.
      real*8 x(2)
      real*8 deltat1,deltat2,factor,lorentz_max,vmp,vmm
      logical cap,floor
      integer NTRIAL
      real*8 TOLF, TOLX, TINY
      parameter (TINY = 1.d-15)
      parameter (NTRIAL = 20)
      parameter (TOLF = 1.d-10, TOLX = 1.d-10)

      lorentz_max = 10.d0

      jl = 0
      kl = 0

      sg = psi**6 * xloc**2
      psi4 = psi**4

      ! Compute rhs quantities.
      c1 = srloc**2 / (sg * psi4)
      c2 = tauloc / sg

      ! Calculate initial guesses.
      w   = dloc/(sg*rholoc)
      eps = uloc/rholoc

c$$$      if (il==1) then
c$$$         write(*,*) 'w = ', w
c$$$         write(*,*) 'sg = ', sg
c$$$         write(*,*) 'rholoc = ', rholoc
c$$$      end if

      ! Make sure things are within range.
      if (w**2 < 1.d0) then
         w = 1.d0
         wbar = 0.d0
      end if
      wbar = sqrt(w**2-1.d0)

c     Branson messing around.
c$$$      if (il==1) then
c$$$         write(*,*) 'eps = ', eps
c$$$         write(*,*) 'wbar = ', wbar
c$$$      end if

      if (wbar == 0.d0 .and. srloc.ne.0.d0) then
         ! We need to do something about the wbar estimate so 
         ! that it doesn't vanish.  Assume small velocities.
         rhohloc = rholoc + ploc + uloc
         vloc = srloc / (lapse*sg*rhohloc)
         utloc = lapse**2 - psi4*(shift+vloc)**2
         if (utloc < 0.d0) then
            write(*,*) 'your inversion fix went superluminal!'
            stop
         end if
         utloc = 1.d0/sqrt(utloc)
         w = lapse*utloc
         if ( (w**2-1.d0) > TINY) then
            wbar = sqrt(w**2-1.d0)
         else
            wbar = 0.d0
         end if
      end if

      x(1) = wbar
      x(2) = eps
      
      if (wbar > 0) then
         call mnewt(NTRIAL,x,2,TOLX,TOLF,retval,c1,c2,dloc,sg,Ye,
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
c     Branson messing around.
      call eos_wrapper(ploc,rholoc,eps,Ye,tloc,rho0,csloc,il,jl,kl,
     &     cap,floor,deltat1,deltat2)
c     Recalculate u in case you just floored the energy.
      uloc = eps * rholoc
      rhohloc = rholoc + ploc + uloc
      vloc = - shift + lapse*srloc/(psi4*sg*rhohloc*w**2)

      if ( (cap.eqv..true.) .or. (floor.eqv..true.)) then
         ! We need to recalculate tau

         ! Note to self.  Don't you think it would be better
         ! to write a little utility subroutine to calculate u^t
         ! instead of inlining it every time?
         utloc  = -lapse**2 + psi4*(shift**2 + 
     &        2.d0*shift*vloc + vloc**2)
         if (utloc .ge. 0.d0) then
            write(*,*) 'superluminality inside table_inversion_1d 
     &           at  = ', il
            factor = (1.d0-1.d0/lorentz_max**2)
            factor = sqrt(factor)
            ! Find the maximum allowed ingoing and outgoing
            ! velocities--for our given speed limit.
            vmp    = lapse/psi**2*factor - shift
            vmm    = -lapse/psi**2*factor - shift
            ! Choose the one closest to our guess.
            factor = abs(vmp - vloc)/abs(vmm - vloc)
            if (factor > 1) then
               vloc = vmm
            else
               vloc = vmp
            end if
            utloc = lorentz_max/lapse
         else
            utloc = 1.d0/sqrt(-utloc)
         end if

         tauloc = sg*rhohloc*(lapse*utloc)**2 - sg*ploc
      end if

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
      subroutine usrfun(x,n,NP,fvec,fjac,c1,c2,dloc,sg,Ye,temp,rho0,
     &    il,jl,kl,cap,floor)
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
      real*8 cs_temp, eps1, eps2
      real*8 deltat1,deltat2
      real*8 w_fac, eps_fac, w_fac_min, eps_fac_min
      logical cap1,floor1,finished_trying,cap2,floor2
      logical cap, floor
      
      deltat1 = 0.d0
      deltat2 = 0.d0

      wbar = x(1)
      eps  = x(2)

      w = sqrt(1.d0+wbar**2)
      rho  = dloc / (sg * w)
      t      = temp

      w_fac = 0.01
      eps_fac = 0.01
      w_fac_min = 1.d-5
      eps_fac_min = 1.d-5

      cap1 = .false.
      cap2 = .false.
      floor1 = .false.
      floor2 = .false.

      ! Make five calls to the eos for our difference star.
      ! Maybe you better check first to make sure that all
      ! of your values are on the table.
      ! Or put this in the wrapper?
      ! We will not be using the temperatures of sounds speeds
      ! that come from these calls.
      call eos_wrapper(P,rho,eps,Ye,t,rho0,cs_temp,il,jl,kl,
     &     cap,floor,deltat1,deltat2)

      if ( (cap.eqv..true.) ) then
         write(*,*) 'central value capped at ', il,
     &        ' ', jl, ' ', kl
         write(*,*) 'w = ', w
         write(*,*) 'eps = ', eps
         x(1) = wbar
         x(2) = eps
         return
      end if

      if ( (floor.eqv..true.)) then
         write(*,*) 'central value floored at ', il,
     &        ' ', jl, ' ', kl
         write(*,*) 'w = ', w
         write(*,*) 'eps = ', eps
         x(1) = wbar
         x(2) = eps
         return
      end if

      finished_trying = .false.
      do while (finished_trying .eqv. .false.)

         if (w_fac .lt. w_fac_min) then
            write(*,*) 'quitting:  w_fac too small at', 
     &           il, ' ', jl, ' ', kl
            write(*,*) 'eps = ', eps
            write(*,*) 'wbar = ', wbar
            stop
         end if

         dw = w * w_fac
         wp   = w   + dw/2.d0
         wm   = w   - dw/2.d0
         rhop = dloc / (sg * wp)
         rhom = dloc / (sg * wm)
         t_wp   = temp
         t_wm   = temp
         eps1 = eps
         eps2 = eps

         call eos_wrapper(P_wp,rhop,eps1,Ye,t_wp,rho0,cs_temp,
     &        il,jl,kl,cap1,floor1,deltat1,deltat2)
         call eos_wrapper(P_wm,rhom,eps2,Ye,t_wm,rho0,cs_temp,
     &        il,jl,kl,cap2,floor2,deltat1,deltat2)

c     Branson messing around.
c$$$         if (il==1) then
c$$$            write(*,*) 'cap1 = ', cap1
c$$$            write(*,*) 'cap2 = ', cap2
c$$$            write(*,*) 'floor1 = ', floor1
c$$$            write(*,*) 'floor2 = ', floor2
c$$$            write(*,*) 'rhop = ', rhop
c$$$            write(*,*) 'rhom = ', rhom
c$$$         end if

         if ( (cap1 .eqv. .false.) .and.
     &        (cap2 .eqv. .false.) .and.
     &        (floor1 .eqv. .false.) .and.
     &        (floor2 .eqv. .false.) ) then
            finished_trying = .true.
         else
            w_fac = w_fac / 2.d0
         end if
      end do

      finished_trying = .false.
      do while (finished_trying .eqv. .false.)

         if (eps_fac .lt. eps_fac_min) then
            write(*,*) 'quitting:  eps_fac too small at', 
     &           il, ' ', jl, ' ', kl
         end if

         deps = eps * eps_fac
         epsp = eps + deps/2.d0
         epsm = eps - deps/2.d0
         t_epsp = temp
         t_epsm = temp
         
         call eos_wrapper(P_epsp,rho,epsp,Ye,t_epsp,rho0,cs_temp,
     &        il,jl,kl,cap1,floor1,deltat1,deltat2)
         call eos_wrapper(P_epsm,rho,epsm,Ye,t_epsm,rho0,cs_temp,
     &        il,jl,kl,cap2,floor2,deltat1,deltat2)
         
         if ( (cap1 .eqv. .false.) .and.
     &        (cap2 .eqv. .false.) .and.
     &        (floor1 .eqv. .false.) .and.
     &        (floor2 .eqv. .false.) ) then
            finished_trying = .true.
         else
            eps_fac = eps_fac / 2.d0
         end if
      end do

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
     &     temp,rho0,il,jl,kl)
      integer n, ntrial, NP, retval
      real*8 tolf,tolx,x(n)
      real*8 c1, c2, dloc, sg, Ye, temp, rho0
      integer il, jl, kl
      parameter (NP=15)
      logical cap, floor
c     USES lubksb, ludcmp, usrfun

      integer i,k,indx(NP)
      real*8 d,errf,errx,fjac(NP,NP),fvec(NP),p(NP)

      retval = 0
      do k=1,ntrial
         
c     Branson messing around.
         do i=1,n
            fvec(i)=0.d0
         end do
         
         call usrfun(x,n,NP,fvec,fjac,c1,c2,dloc,sg,Ye,temp,rho0,il,
     &        jl,kl,cap,floor)
         
         if ((cap .eqv. .true.) .or. (floor .eqv. .true.)) then
            return
         else 
            errf=0.d0
            do i=1,n
               errf=errf+abs(fvec(i))
            enddo

            if (errf .le. tolf) return
            
            do i=1,n
               p(i)=-fvec(i)
            enddo

c     Solve linear equations using LU decomposition.
            call ludcmp(fjac,n,NP,indx,d,il,jl,kl)
            call lubksb(fjac,n,NP,indx,p)
            errx = 0.d0
            
            do i=1,n
               errx=errx+abs(p(i))
               x(i)=x(i)+p(i)
            enddo
            
c     Branson messing around.
            if (x(2) < 0.d0) x(2) = 0.d0
            
            if (errx .le. tolx) return

         end if
      end do  

      if (k==ntrial) retval = 1

      return
      end 

      subroutine ludcmp(a,n,np,indx,d,il,jl,kl)
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
            

c-------------------------------------------------------------
c     For this wrapper, the user provides density, temperature
c     and Ye.
c-------------------------------------------------------------
      subroutine eos_wrapper_temp(p,rho,energy,Ye,temp,rho0,cs,
     &     il,jl,kl)
      implicit none
      real*8 p, rho, energy, ye, cs
      integer inputz, il,jl,kl
      real*8 dhyx,zhtx,q,delty
      real*8 csx,temp
      real*8 eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,gmx
      real*8 var,mevcon,muex,dpdex,abar,zbar
      logical keyerr, cap, floor
      real*8 rho0,c
      real*8 deltat1,deltat2
      parameter (c = 2.99792458d10)

      keyerr = .false.
      inputz = 1

c     these variables are irrelevant.
      deltat1 = 0.d0
      deltat2 = 0.d0

c     convert units
      rho = rho * rho0
      var = temp
      
      call find_interpLSeos(inputz,var,temp,rho,ye,
     &     eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,
     &     gmx,dhyx,zhtx,csx,muex,dpdex,keyerr,il,jl,kl,
     &     cap,floor,deltat1,deltat2)

      p = prx

c     Convert back to code units.
      rho = rho / rho0
      energy = eix / c**2
      p = p / (rho0 * c**2)
      cs = csx

      return
      end subroutine eos_wrapper_temp

c-------------------------------------------------------------
c     For this wrapper, the user provides density, pressure
c     and Ye.  This is used for pressure depletion runs.
c-------------------------------------------------------------
      subroutine eos_wrapper_press(p,rho,energy,Ye,temp,rho0,cs,
     &     il,jl,kl,cap,floor,deltat1,deltat2)
      implicit none
      real*8 p, rho, energy, ye, cs
      integer inputz,il,jl,kl
      real*8 dhyx,zhtx,q,delty
      real*8 csx,temp
      real*8 eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,gmx
      real*8 var,mevcon,muex,dpdex,abar,zbar
      real*8 deltat1,deltat2
      logical keyerr,cap,floor
      real*8 rho0,c
      parameter (c = 2.99792458d10)

      keyerr = .false.
      inputz = 4

c     convert units
      rho = rho * rho0
      p = p * c**2*rho0
      var = p
      
      call find_interpLSeos(inputz,var,temp,rho,ye,
     &     eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,
     &     gmx,dhyx,zhtx,csx,muex,dpdex,keyerr,il,jl,kl,cap,floor,
     &     deltat1,deltat2)

      p = prx

c     Convert back to code units.
      rho = rho / rho0
      energy = eix / c**2
      p = p / (rho0 * c**2)
      cs = csx

      return
      end subroutine

c---------------------------------------------------------------------
c     For this wrapper, user provides internal energy, density, and Ye.
c---------------------------------------------------------------------
      subroutine eos_wrapper(p,rho,energy,Ye,temp,rho0,cs,il,jl,kl,
     &     cap,floor,deltat1,deltat2)
      implicit none
      real*8 p, rho, energy, ye, cs
      integer inputz
      real*8 dhyx,zhtx,q,delty
      real*8 csx,temp
      real*8 eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,gmx
      real*8 var,mevcon,muex,dpdex,abar,zbar
      real*8 deltat1,deltat2
      logical keyerr,cap,floor
      real*8 rho0,c
      integer il,jl,kl
      parameter (c = 2.99792458d10)

      keyerr = .false.
      inputz = 2

c     convert units
      rho = rho * rho0
      energy = energy * c**2
      var = energy
      
      call find_interpLSeos(inputz,var,temp,rho,ye,
     &     eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,
     &     gmx,dhyx,zhtx,csx,muex,dpdex,keyerr,il,jl,kl,
     &     cap,floor,deltat1,deltat2)

      p = prx

c     Convert back to code units.
      rho = rho / rho0
      energy = eix / c**2
      p = p / (rho0 * c**2)
      cs = csx 

      return
      end subroutine eos_wrapper

c-----------------------------------------------------------------------
c     Branson has modified this version so that, if the temperature
c     ends up being outside the desired range, it will be capped or
c     floored.  The appropriate logical flag will be triggered and
c     passed back to the calling function.
c-----------------------------------------------------------------------
      subroutine find_interpLSeos(inputz,var,temp,rho,ye,
     .     eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,
     .     gmx,dhyx,zhtx,csx,muex,dpdex,keyerr,il,jl,kl,cap,floor,
     .     deltat1,deltat2)
      implicit none
      double precision temp,rho,ye,dhyx,zhtx,q,p,delty,qinter,energy
      double precision rl,ql,alpha,beta,delta,linter,q0,q1,q2,qq
      double precision dhy0,zht0,dhy1,zht1,dhy2,zht2,yl0,yl1,yl2,csx
      double precision eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,gmx
      double precision entropy,f,f1,f2,var,mevcon,muex,dpdex,rhloc,c
      double precision deltat1,deltat2,tf,tc, exp_tf, exp_tc
      parameter (mevcon=0.95655684d18)
      parameter (c = 2.99792458d10)
      integer i,j,k,jr,jq,jy,nn,inputz,il,jl,kl
      logical keyerr,cap,floor
      common /quantity/ qq(50)
      include 'table_all_params.ecomp-2'
      temp   = temp/1.160445d10
      if(ye.lt.0.052d0)ye=0.05d0
      if(inputz.eq.1)then
         do nn=1,numel
            call findthis(nn,temp,rho,ye,f)
            qq(nn)=f
         enddo
         eix = qq(1)
         enx = qq(3)
      else if(inputz.eq.2)then
         energy=var
         energy=energy/mevcon - 9.3d0
c     Check whether energy is within range.
c     The \deltat T variables adjust the range.
         exp_tf = t1 + deltat1*(t2-t1)
         exp_tc = t2 - deltat2*(t2-t1)
         tf = 10.d0**exp_tf
         tc = 10.d0**exp_tc
         nn = 1
         call findthis(nn,tf,rho,ye,f1)
         call findthis(nn,tc,rho,ye,f2)
         cap = .false.
         floor = .false.
         if (energy < f1) then

c     Branson messing around.
c$$$            if (il==1) then
c$$$               write(*,*) 'rho = ', rho
c$$$               write(*,*) 'energy = ', energy
c$$$               write(*,*) 'f1 = ', f1
c$$$               write(*,*) 'f2 = ', f2
c$$$            end if

            floor = .true.
            temp = tf
            energy = f1

         else if (energy > f2) then
            cap = .true.
            temp = tc
            energy = f2
         else
            call findtemp(inputz,energy,temp,rho,ye,keyerr,il,jl,kl)
         end if
         do nn=1,numel
            call findthis(nn,temp,rho,ye,f)
            qq(nn)=f
         enddo
         eix = energy
         enx = qq(3)
      else if(inputz.eq.3)then
         entropy=var
         call findtemp(inputz,entropy,temp,rho,ye,keyerr,il,jl,kl)
         do nn=1,numel
            call findthis(nn,temp,rho,ye,f)
            qq(nn)=f
         enddo
         enx=entropy
         eix=qq(1)
      else if(inputz.eq.4)then
         p=var
         p=p/1.60217733d33
c     Check whether pressure is within range.
c     The \deltat T variables adjust the range.
         exp_tf = t1 + deltat1*(t2-t1)
         exp_tc = t2 - deltat2*(t2-t1)
         tf = 10.d0**exp_tf
         tc = 10.d0**exp_tc
         nn = 2
         call findthis(nn,tf,rho,ye,f1)
         call findthis(nn,tc,rho,ye,f2)
         cap = .false.
         floor = .false.
         if (p < f1) then
            floor = .true.
            temp = tf
            p = f1
         else if (p > f2) then
            cap = .true.
            temp = tc
            p = f2
         else
            call findtemp(inputz,p,temp,rho,ye,keyerr,il,jl,kl)
         end if
         do nn=1,numel
            call findthis(nn,temp,rho,ye,f)
            qq(nn)=f
         enddo
         prx = p
         enx = qq(3)
         eix = qq(1)
      endif   
      prx   = qq(2)
      cvx   = qq(4)
      xnx   = qq(5)
      xpx   = qq(6)
      xax   = qq(7)
      xhx   = qq(8)
      zax   = qq(9)
      awx   = qq(10)
      mhx   = qq(11)
      gmx   = qq(12)
      dhyx  = qq(13)
      zhtx  = qq(14)
      muex  = qq(15)
      dpdex = qq(16)
c
      temp = temp*1.160445d10
      prx  = prx*1.60217733d33
      eix  = (eix + 9.3d0)*mevcon
      cvx  = cvx/1.60217733d-6*1.381d-16
      rhloc = rho*(c**2+eix) + prx
      csx  = sqrt(prx/rhloc*gmx)
c
      if(keyerr)then
         print*,'Problem in Findtemp: keyerr',keyerr
         return
      endif

      return
      end

c---------------------------------------------------------------------
      subroutine shen_initialize
      implicit none
      integer jy,jr,jt,irec,irecl
c     include 'table_all_params.shen'
      include 'table_all_params.ecomp-2'
      irecl=8*ny*nr*nt

      open(10,file='newtable18018050.shen-grid.2.den1.dat',
     &     form='unformatted',access='direct',recl=irecl)
      do irec=1,numel
         read(10,rec=irec)(((table(irec,jy,jr,jt),jt = 1,nt),
     &        jr = 1,nr),jy = 1,ny)
      enddo

      call fix_table

      return
      end
c----------------------------------------------------------------------
c     Branson wrote this function to interpolate over some nans which
c     are inherent in the table.  Instead of writing over the data file
c     itself, we can just call this once and for all each time we 
c     call shen_initialize.
c----------------------------------------------------------------------
      subroutine fix_table
      implicit none
      integer jy,jr,jt,i,j,k,klow,kup,l
      real*8 gmxloc,gmxlow,gmxup
      include 'table_all_params.ecomp-2'

c     Go to each density and Ye point.
      do j=1,ny
         do i=1,nr
            do k=1,nt
               gmxloc = table(12,j,i,k)
               if (gmxloc/gmxloc .ne. 1.d0) then
c     Search for the nearest non-NaN points
c     in temperature.
                  gmxok(k) = .false.
                  do l=k-1,1,-1
                     gmxlow = table(12,j,i,l)
                     if (gmxok(l)) then
                        klow = l
                        exit
                     end if
                  end do

                  do l=k+1,nt
                     gmxup = table(12,j,i,l)
                     if (gmxup/gmxup .eq. 1.d0) then
                        kup = l
                        exit
                     end if
                  end do
c     Do linear interpolation in log T.
                  gmxloc = gmxlow + (gmxup-gmxlow)/(kup-klow)*(k-klow)
                  table(12,j,i,k) = gmxloc
               else
                  gmxok(k) = .true.
               end if
            end do
         end do
      end do

      end subroutine fix_table
c----------------------------------------------------------------------
      subroutine findtemp(inputz,var,temp,rho,ye,keyerr,il,jl,kl)
      implicit none
      double precision ep,temp,rho,ye,e,ep1,ep2,dt,depdt,dd,f,var
      double precision q,p,delty,qinter,energy,yl0,yl1,yl2,told
      double precision rl,ql,alpha,beta,delta,linter,q0,q1,q2,qq
      double precision temp0,xxx,deltat1,deltat2
      integer i,jy,jr,jq,nn,icount,inputz,index,nnnstep,il,jl,kl
      logical keyerr,gobisect
      common /indices/ jy,jr,jq
      common /xxpass/ index,nnnstep,icount
c     include 'table_all_params.shen'
      include 'table_all_params.ecomp-2'
      
      temp0=temp
      if(inputz.eq.2)nn=1
      if(inputz.eq.3)nn=3
      if(inputz.eq.4)nn=2
      icount=0
      gobisect=.false.
      do i = 1,15
         call findthis(nn,temp,rho,ye,f)
         if(jq.lt.2)then
            gobisect=.true.
            goto 10
         endif
         ep1    = f
         dt     = 0.01d0*temp
         temp   = temp+dt
         call findthis(nn,temp,rho,ye,f)
         ep2    = f
         temp   = temp-dt
         depdt  = (ep2-ep1)/dt
         dd     = (var-ep1)/depdt
         icount = icount + 1
         temp   = temp+dd
         if(temp.le.0.d0) temp=temp0+0.1d0*temp0
c     Branson:  Changing the tolerance here.
c         if(abs(dd/temp).lt.1.d-8)goto 10
         if(abs(dd/temp).lt.1.d-6)goto 10
      enddo
 10   continue
      if(abs(temp-temp0)/temp0.gt.0.5d0.or.icount.eq.15.or.gobisect)then
         call bisection(inputz,var,temp0,rho,ye,il,jl,kl)
         temp=temp0
      endif   
      return
      end
c-------------------------------------------------------------------------
c     Given the pressure, Ye, and temperature, we need to find the 
c     rest-mass density.  This is useful when integrating up the 
c     TOV equations.
c-------------------------------------------------------------------------
      subroutine finddens_wrapper(p,rho,energy,Ye,temp,rho0,cs,
     &     il,jl,kl)
      implicit none
      real*8 p, rho, energy, ye, cs
      integer il,jl,kl
      real*8 dhyx,zhtx,q,delty
      real*8 csx,temp
      integer inputz
      real*8 eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,gmx
      real*8 var,mevcon,muex,dpdex,abar,zbar
      logical keyerr,cap,floor
      real*8 deltat1, deltat2
      real*8 rho0,c
      parameter (c = 2.99792458d10)

      keyerr = .false.
c     These deltat variables are not relevant in the present case.
      deltat1 = 0.d0
      deltat2 = 0.d0

c     convert units
      p = p * (rho0 * c**2)
      p = p/1.60217733d33
      temp   = temp/1.160445d10
      rho = rho * rho0

      call finddens(p,temp,rho,ye,keyerr,il,jl,kl)

      temp   = temp * 1.160445d10

c     Now that you have density, find the rest of the quantities.
      inputz = 1
      var = temp
      
      call find_interpLSeos(inputz,var,temp,rho,ye,
     &     eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,
     &     gmx,dhyx,zhtx,csx,muex,dpdex,keyerr,il,jl,kl,
     &     cap,floor,deltat1,deltat2)

c     Convert back to code units.
      rho = rho / rho0
      energy = eix / c**2
      p = prx / (rho0 * c**2)
      cs = csx

      return
      end subroutine finddens_wrapper
c--------------------------------------------------------------------
      subroutine finddens(p,temp,rho,ye,keyerr,il,jl,kl)
      implicit none
      double precision p,temp,rho,ye,p1,p2,dr,dpdr,dd,f
      double precision rho0
      integer i,jy,jr,jq,nn,icount,inputz,index,nnnstep,il,jl,kl
      logical keyerr,gobisect
      common /indices/ jy,jr,jq
      common /xxpass/ index,nnnstep,icount
      include 'table_all_params.ecomp-2'
      
      rho0=rho
      nn = 2
      icount=0
      gobisect=.false.
      do i = 1,15
         call findthis(nn,temp,rho,ye,f)
         if(jq.lt.2)then
            gobisect=.true.
            goto 10
         endif
         p1    = f
         dr    = 0.01d0*rho
         rho   = rho+dr
         call findthis(nn,temp,rho,ye,f)
         p2    = f
         rho   = rho-dr
         dpdr  = (p2-p1)/dr
         dd    = (p-p1)/dpdr
         icount = icount + 1
         rho   = rho+dd

         if(rho.le.0.d0) rho=rho0+0.1d0*rho0
         if(abs(dd/rho).lt.1.d-8)goto 10
      enddo
 10   continue
      if(icount.eq.15.or.gobisect)then
         call bisect_dens(p,temp,rho0,ye,il,jl,kl)
         rho=rho0
      endif   
      return
      end
c----------------------------------------------------------------------
      subroutine bisection(inputz,var,temp,rho,ye,il,jl,kl)
      implicit none
      double precision ep,temp,rho,ye,e,ep1,ep2,dt,depdt,dd,f,var
      double precision q,p,delty,qinter,energy,yl0,yl1,yl2,told
      double precision rl,ql,alpha,beta,delta,linter,q0,q1,q2,qq
      double precision temp1,temp2,f1,f2,dtemp,tmid,fmid,temp00
      integer i,jy,jr,jq,nn,icount,inputz,index,nnnstep,mcount
      integer il,jl,kl
      common /indices/ jy,jr,jq
      common /xxpass/ index,nnnstep,icount
      logical cap, floor
      include 'table_all_params.ecomp-2'
      cap = .false.
      floor = .false.
      if(inputz.eq.2)nn=1
      if(inputz.eq.3)nn=3
      if(inputz.eq.4)nn=2
      temp00=temp
      mcount=0
      temp1=temp-0.1d0*temp
      temp2=temp+0.1d0*temp
      call findthis(nn,temp,rho,ye,f)
      f=f-var
      call findthis(nn,temp1,rho,ye,f1)
      call findthis(nn,temp2,rho,ye,f2)
      f1=f1-var
      f2=f2-var
 2    continue
      if(f1*f2.ge.0.d0)then
         mcount=mcount+1
         temp1=temp1-0.1d0*temp1
         temp2=temp2+0.1d0*temp2
         call findthis(nn,temp1,rho,ye,f1)
         call findthis(nn,temp2,rho,ye,f2)
         f1=f1-var
         f2=f2-var
         if(mcount.le.20)goto 2
         if(mcount.gt.20)then
c     Time to try the whole range.
            temp1 = 10.d0**t1
            temp2 = 10.d0**t2
            call findthis(nn,temp1,rho,ye,f1)
            call findthis(nn,temp2,rho,ye,f2)
            f1=f1-var
            f2=f2-var
            if (f1*f2 .ge. 0.d0) then
               write(*,*) 'findtemp giving up at', il, ' ', 
     &              jl, ' ', kl
               stop
            end if
         endif
      endif
 28   format(I4,I5,I5,9e13.5)
      if(f1.lt.0.d0)then
         temp=temp1
         dtemp=temp2-temp1
      else
         temp=temp2
         dtemp=temp1-temp2
      endif 
      do i=1,40
         dtemp=dtemp*0.5d0
         tmid=temp+dtemp
         call findthis(nn,tmid,rho,ye,fmid)
         fmid=fmid-var
         if(fmid.le.0.d0)temp=tmid
         if(abs(dtemp).lt.1.d-6)goto 10
      enddo
      temp=temp00
 10   continue
      temp=temp
      return
      end
c-----------------------------------------------------------------
c     This is similar to the bisection subroutine for finding the
c     temperature, but here it tries to find the rest-mass density
c     given P, Ye, and T.
c-----------------------------------------------------------------
      subroutine bisect_dens(p,temp,rho,ye,il,jl,kl)
      implicit none
      double precision temp,rho,ye,f, p
      double precision rho1,rho2,f1,f2,drho,rmid,fmid,rho00
      integer i,jy,jr,jq,nn,icount,inputz,index,nnnstep,mcount
      integer il,jl,kl
      common /indices/ jy,jr,jq
      common /xxpass/ index,nnnstep,icount
      include 'table_all_params.ecomp-2'
      nn = 2
      rho00=rho
      mcount=0
      rho1=rho-0.1d0*rho
      rho2=rho+0.1d0*rho
      call findthis(nn,temp,rho,ye,f)
      f=f-p
      call findthis(nn,temp,rho1,ye,f1)
      call findthis(nn,temp,rho2,ye,f2)
      f1=f1-p
      f2=f2-p
 2    continue
      if(f1*f2.ge.0.d0)then
         mcount=mcount+1
         rho1=rho1-0.1d0*rho1
         rho2=rho2+0.1d0*rho2
         call findthis(nn,temp,rho1,ye,f1)
         call findthis(nn,temp,rho2,ye,f2)
         f1=f1-p
         f2=f2-p
         if(mcount.le.40)goto 2
         if(mcount.gt.40)then
            rmid=rho00
            write(*,*) 'finddens giving up at ', il, ' ', jl, ' ', kl
            stop
         endif
      endif
 28   format(I4,I5,I5,9e13.5)
      if(f1.lt.0.d0)then
         rho=rho1
         drho=rho2-rho1
      else
         rho=rho2
         drho=rho1-rho2
      endif 
      do i=1,40
         drho=drho*0.5d0
         rmid=rho+drho
         call findthis(nn,temp,rmid,ye,fmid)
         fmid=fmid-p
         if(fmid.le.0.d0)rho=rmid
         if(abs(drho).lt.1.d-6)goto 10
      enddo
      rho=rho00
 10   continue
      rho=rho
      return
      end
c----------------------------------------------------------------------
      subroutine offtable
      implicit none
      integer jy,jr,jq
      common /indices/ jy,jr,jq
      include 'table_all_params.ecomp-2'
      if(jy.le.2)jy=2
      if(jr.le.2)jr=2
      if(jq.le.2)jq=2
      if(jy.ge.ny)jy=ny-1
      if(jr.ge.nr)jr=nr-1
      if(jq.ge.nt)jq=nt-1
      return
      end

c----------------------------------------------------------------------
      subroutine print_range(rho,ye,rho0)
      implicit none
      double precision rho, ye
      double precision temp, eps
      double precision mevcon, c
      double precision rho0
      parameter (mevcon=0.95655684d18)
      parameter (c = 2.99792458d10)
      integer i
      include 'table_all_params.ecomp-2'

      rho = rho * rho0
      do i = 1, nt
         temp = dfloat(i-1)/nt*(t2 - t1) + t1
         temp = 10.d0**temp 
         eps = 0.d0
         call findthis(1,temp,rho,ye,eps)
         temp = temp * 1.160445d10
         eps = (eps + 9.3d0)*mevcon / c**2
         write(*,*) i, temp, eps
      end do
      rho = rho / rho0
      end subroutine
c----------------------------------------------------------------------
      subroutine findthis (nn,temp,rho,ye,f)
      implicit none
      double precision ep,temp,rho,ye,e,ep1,ep2,dt,depdt,dd,f,tl
      double precision q,p,delty,qinter,energy,yl0,yl1,yl2,told
      double precision rl,ql,alpha,beta,delta,linter,q0,q1,q2,qq
      double precision rl0,rl1,rl2,tl0,tl1,tl2,t,u
      integer jy,jr,jq,nn,icount
      common /indices/ jy,jr,jq
      include 'table_all_params.ecomp-2'
      rl    = dlog10(rho)
      tl    = dlog10(temp)
      alpha = t1+(rl-r1)/(r2-r1)*(t12-t1)
      beta  = t2-t1+((t22-t12)-(t2-t1))*(rl-r1)/(r2-r1)
      ql    = (tl - alpha)/beta
      delta = (rl-r1)/(r2-r1)*dfloat(nr-1)
      delty = (ye-y1)/(y2-y1)*dfloat(ny-1)
      jy    = 1 + idint(delty)
      jr    = 1 + idint(delta)
      jq = 1 + idint(dfloat(nt-1)*ql)

      if(jy.le.2.or.jy.ge.ny-1.or.jr.le.2.or.jr.ge.nr-1.or.
     .   jq.le.2.or.jq.ge.nt-1) call offtable
      p     = delta - (jr-1)
      q     = dfloat(nt-1)*ql - (jq-1)
      yl0 = y1+(y2-y1)*dfloat(jy-1-1)/dfloat(ny-1)
      yl1 = y1+(y2-y1)*dfloat(jy-1)/dfloat(ny-1)
      yl2 = y1+(y2-y1)*dfloat(jy+1-1)/dfloat(ny-1)
c     if(nn.eq.5.or.nn.eq.6.or.nn.eq.7.or.nn.eq.8)then
      if(nn.ne.20) then
         q0 = (1.D0-P)*(1.d0-Q)*table(nn,jy-1,jr,jq)
     .        + P*(1.D0-Q)*table(nn,jy-1,jr+1,jq)
     .        + Q*(1.D0-P)*table(nn,jy-1,jr,jq+1)
     .        + P*Q*table(nn,jy-1,jr+1,jq+1)
         q1 = (1.D0-P)*(1.d0-Q)*table(nn,jy,jr,jq)
     .        + P*(1.D0-Q)*table(nn,jy,jr+1,jq)
     .        + Q*(1.D0-P)*table(nn,jy,jr,jq+1)
     .        + P*Q*table(nn,jy,jr+1,jq+1)
         q2 = (1.D0-P)*(1.d0-Q)*table(nn,jy+1,jr,jq)
     .        + P*(1.D0-Q)*table(nn,jy+1,jr+1,jq)
     .        + Q*(1.D0-P)*table(nn,jy+1,jr,jq+1)
     .        + P*Q*table(nn,jy+1,jr+1,jq+1)
         f  = linter(ye,yl1,q1,yl2,q2)
         if(nn.eq.5.or.nn.eq.6.or.nn.eq.7.or.nn.eq.8) 
     &      f  = min(max(f,0.d0),1.d0)
      else    
         q0 = 0.5D0*Q*(Q-1.D0)*table(nn,jy-1,jr,jq-1)
     .        + 0.5D0*P*(P-1.D0)*table(nn,jy-1,jr-1,jq)
     .        + (1.D0+P*Q-P*P-Q*Q)*table(nn,jy-1,jr,jq)
     .        + 0.5D0*P*(P-2.D0*Q+1.D0)*table(nn,jy-1,jr+1,jq)
     .        + 0.5D0*Q*(Q-2.D0*P+1.D0)*table(nn,jy-1,jr,jq+1)
     .        + P*Q*table(nn,jy-1,jr+1,jq+1)
         q1 = 0.5D0*Q*(Q-1.D0)*table(nn,jy,jr,jq-1)
     .        + 0.5D0*P*(P-1.D0)*table(nn,jy,jr-1,jq)
     .        + (1.D0+P*Q-P*P-Q*Q)*table(nn,jy,jr,jq)
     .        + 0.5D0*P*(P-2.D0*Q+1.D0)*table(nn,jy,jr+1,jq)
     .        + 0.5D0*Q*(Q-2.D0*P+1.D0)*table(nn,jy,jr,jq+1)
     .        + P*Q*table(nn,jy,jr+1,jq+1)
         q2 = 0.5D0*Q*(Q-1.D0)*table(nn,jy+1,jr,jq-1)
     .        + 0.5D0*P*(P-1.D0)*table(nn,jy+1,jr-1,jq)
     .        + (1.D0+P*Q-P*P-Q*Q)*table(nn,jy+1,jr,jq)
     .        + 0.5D0*P*(P-2.D0*Q+1.D0)*table(nn,jy+1,jr+1,jq)
     .        + 0.5D0*Q*(Q-2.D0*P+1.D0)*table(nn,jy+1,jr,jq+1)
     .        + P*Q*table(nn,jy+1,jr+1,jq+1)
         f  = qinter(ye,yl0,q0,yl1,q1,yl2,q2)
      endif
      return
      end
c----------------------------------------------------------------------
      double precision function qinter(x,x1,y1,x2,y2,x3,y3)
      implicit none
      double precision x,x1,y1,x2,y2,x3,y3,y
      qinter = y1*(x-x2)*(x-x3)/(x1-x2)/(x1-x3)+
     .         y2*(x-x1)*(x-x3)/(x2-x1)/(x2-x3)+
     .         y3*(x-x1)*(x-x2)/(x3-x1)/(x3-x2)
      return
      end
c----------------------------------------------------------------------
      double precision function linter(x,x0,y0,x1,y1)
      implicit none
      double precision x0,y0,x1,y1,x,y
      linter = y0+(y1-y0)/(x1-x0)*(x-x0)
      return
      end

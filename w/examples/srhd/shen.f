      subroutine eos_wrapper(p,rho,energy,Ye,temp,rho0,cs)
      implicit none
      real*8 p, rho, energy, ye, cs
      integer inputz
      real*8 dhyx,zhtx,q,delty
      real*8 csx,temp
      real*8 eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,gmx
      real*8 var,mevcon,muex,dpdex,abar,zbar
      logical keyerr
      real*8 rho0,c
      parameter (c = 2.99792458d10)

      keyerr = .false.
      inputz = 2

c     convert units
      rho = rho * rho0
      energy = energy * c**2
      var = energy
      
      call find_interpLSeos(inputz,var,temp,rho,ye,
     &     eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,
     &     gmx,dhyx,zhtx,csx,muex,dpdex,keyerr)

      p = prx

c     Convert back to code units.
      rho = rho / rho0
      energy = energy / c**2
      p = p / (rho0 * c**2)
c      cs = csx / c
      cs = csx

      return
      end subroutine eos_wrapper

c     Branson messing around.  A version with info about location
      subroutine eos_wrapper_loc(p,rho,energy,Ye,temp,rho0,cs,il,jl,kl)
      implicit none
      real*8 p, rho, energy, ye, cs
      integer inputz
      real*8 dhyx,zhtx,q,delty
      real*8 csx,temp
      real*8 eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,gmx
      real*8 var,mevcon,muex,dpdex,abar,zbar
      logical keyerr
      real*8 rho0,c
      integer il,jl,kl
      parameter (c = 2.99792458d10)

      keyerr = .false.
      inputz = 2

c     convert units
      rho = rho * rho0
      energy = energy * c**2
      var = energy
      
      call find_interpLSeos_loc(inputz,var,temp,rho,ye,
     &     eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,
     &     gmx,dhyx,zhtx,csx,muex,dpdex,keyerr,il,jl,kl)

      p = prx

c     Convert back to code units.
      rho = rho / rho0
      energy = energy / c**2
      p = p / (rho0 * c**2)
c      cs = csx / c
      cs = csx 

      return
      end subroutine eos_wrapper_loc

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     The stuff after this is from Burrows.
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine driver (nnstep,inputz)
      implicit none
      double precision e,t,r,y,eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx
      double precision mhx,gmx,yein,rin,tin,sin,ein,pin,gamin,dhyx,zhtx
      double precision rhonew,rhoold,dv,pnew,pold,thep,told,unew,uold
      double precision sinit,rhoinit,var,emu,ehat,aw,ye,za,xh,xxn,xa
      double precision xpro,massden,temp0,etae,cv,yl,etanu,emnu,znu
      double precision dhy,zht,dxndt,dxndy,dxprodt,dxprody,detotdye
      double precision dpdr,rhotrans,drhotrans,prh,eih,enh,prl,eil,enl
      double precision linter,tl,gl,th,gh,muex,csx,dpdex,abar,zbar
      double precision rhom, rhop
      integer i,i2,inputz,nnstep,nnx,nrmax,index,nnnstep,in,icount
      integer nrmaxp, nrmax3, maxcore
      integer iminus,iplus,inputzk,iflag,ndex
      logical keyerr
      include 'specttnr'
      common /eosin/ yein(nrmax),rin(nrmax),tin(nrmax),sin(nrmax),
     .  ein(nrmax),pin(nrmax),gamin(nrmax)
      common /themp/ emu(nrmax), ehat(nrmax),aw(nrmax),ye(nrmax),
     *  za(nrmax), xh(nrmax), xxn(nrmax), xa(nrmax), xpro(nrmax),
     *  massden(nrmax), temp0(nrmax), etae(nrmax),
     *  cv(nrmax), yl(nrmax), etanu(nrmax), emnu(nrmax), znu(nrmax),
     *  dhy(nrmax), zht(nrmax), dxndt(nrmax), dxndy(nrmax),
     *  dxprodt(nrmax), dxprody(nrmax), detotdye(nrmax) 
      common /size/ nnx
      common /xxpass/ index,nnnstep,icount
      common /drivela/ idrive
      integer idrive
      nnnstep   = nnstep
      keyerr=.false.
c
c     rhom = 1.39d14
c     rhop = 1.92d14
      rhom = 1.39d16
      rhop = 1.92d16
c     Branson:  I commented these out because I don't have ndex....
c      iminus = min(nnx,ndex(rhom,rin,nnx) + 2)  !  Could be + 1
c      iplus  = max(1,ndex(rhop,rin,nnx) - 1)  !  -1 not necessary
c
      do i=1,nnx
         index=i
         t = tin(i)
         r = rin(i)
         y = yein(i)
         if(inputz.eq.1)var=t*1.160445d10
         if(inputz.eq.2)var=ein(i)
         if(inputz.eq.3)var=sin(i)
c
         iflag = 0
         if((r .ge. rhom) .and. (r .lt. rhop) .and. (iminus .ne. iplus) 
     &        .and. (iplus .gt. 5)) then
           inputzk = inputz
           t = log(tin(iminus)) + (log(tin(iplus))-log(tin(iminus)))*
     &         (log(r)-log(rin(iminus)))/(log(rin(iplus))
     &          -log(rin(iminus)))
           t = exp(t)
c          write(6,*) t,r,iminus,iplus
           var = t*1.160445d10
           inputz = 1
           iflag = 1
         endif
c
         call getxs(abar,zbar)
         t=t*1.160445d10
         if(idrive .eq. 1) then
           call find_interp(inputz,var,t,r,y,abar,zbar,
     .        eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,
     .        gmx,dhyx,zhtx,csx,muex,dpdex,keyerr)
         else if(idrive .eq. 0) then
           call find_interpLSeos(inputz,var,t,r,y,
     .        eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,
     .        gmx,dhyx,zhtx,csx,muex,dpdex,keyerr)
         else
           call find_interpHeos(inputz,var,t,r,abar,zbar,
     .        eix,prx,enx,cvx,gmx,csx,muex,dpdex,keyerr)
         endif
c
         if(iflag .eq. 1) inputz = inputzk
c
         t=t/1.160445d10
         tin   (i) = t
         ein   (i) = eix
         pin   (i) = prx
         sin   (i) = enx
         gamin (i) = gmx
         cv    (i) = cvx*1.60217733d-6/1.381d-16
         xxn   (i) = xnx   
         xpro  (i) = xpx
         xa    (i) = xax
         xh    (i) = xhx
         aw    (i) = min (max (50.d0, awx), 500.d0)
         za    (i) = zax
         ehat  (i) = mhx
         dhy   (i) = dhyx
         zht   (i) = zhtx
         emu   (i) = muex
      enddo
      return
      end
c---------------------------------------------------------------------
      subroutine getxs(abar,zbar)
      include 'hydrocom'
      double precision xxm,xxa,xxz,ymass,abar,zbar
      integer index,nnnstep,icount,ionmax
      parameter (ionmax=4)
      dimension xxm(ionmax),xxa(ionmax),xxz(ionmax),ymass(ionmax)
      common /xxpass/ index,nnnstep,icount
      xxm(1) = xnuc(index,2)
      xxm(2) = xnuc(index,3)
      xxm(3) = xnuc(index,4)
      xxm(4) = xnuc(index,5)      
      xxa(1) = 1.d0
      xxz(1) = 1.d0
      xxa(2) = 4.d0
      xxz(2) = 2.d0
      xxa(3) = 16.d0
      xxz(3) = 8.d0
      xxa(4) = 56.d0
      xxz(4) = 26.d0
      call azbar(xxm,xxa,xxz,ionmax,ymass,abar,zbar)
      return
      end
c------------------------------------------------------------------------
c      subroutine driver (ind,nstep,inputz,var,t,r,y,
c     .  eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,
c     .  gmx,dhyx,zhtx,muex,dpdex)
c
c      implicit none
c
c      double precision var,t,r,y,abar,zbar,
c     .  eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,
c     .  gmx,dhyx,zhtx,muex,dpdex,csx
c      integer ind,nstep,inputz,index,nnstep,icount
c      common /xxpass/ index,nnstep,icount
c
c      nnstep = nstep
c      index  = ind
c      call azbar(xxm,xxa,xxz,ionmax,ymass,abar,zbar)
c      abar=56.d0
c      zbar=26.d0
c      t=t*1.160445d10
c      call find_interpLSeos(inputz,var,t,r,y,
c     .    eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,
c     .    gmx,dhyx,zhtx,csx,muex,dpdex)
c      call find_interp(inputz,var,t,r,y,abar,zbar,
c     .     eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,
c     .     gmx,dhyx,zhtx,csx,muex,dpdex)
c
c      t=t/1.160445d10
c      return
c      end
c-------------------------------------------------------------------------------
c                              FIND_INTERP
c-------------------------------------------------------------------------------
      subroutine find_interp(inputz,var,temp,rho,ye,abar,zbar,
     .    eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,
     .    gmx,dhyx,zhtx,csx,muex,dpdex,keyerr)
      implicit none
      double precision var,temp,rho,ye,linter,eix,prx,enx,cvx,xnx,xpx
      double precision xax,xhx,zax,awx,mhx,gmx,dhyx,zhtx,prh,prl,enh,enl
      double precision rtrans,drtrans,eih,eil,th,tl,gmh,gml,ch,cl,sound
      double precision abar,zbar,xxm,xxa,xxz,ymass,ptot,stot,etot,cv,dd
      double precision temp0,etot1,dedt,cv1,csx,csl,csh,ematch,gl,gh
      double precision prm,eim,cvm,csm,enm,gmm,ttm,qinter,r1,r2,r3
      double precision muex,dpdex,mum,mul,muh,dxh,dxl,dxm,abarl,zbarl
      double precision lseix,heix
      integer inputz,ionmax,index,i,nnnstep,icount
      logical keyerr
c
      ematch = 0.0887d0*4.d18
c
c .. transition density and dtransition density for interpolation between
c     lattimer-swesty and helmholtz eos
c
c     rtrans  = 1.0d5
      rtrans  = 1.0d4
      drtrans = 3.d3
c
c     call find_interpLSeos(inputz,var,temp,rtrans,ye,
c    .     eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,
c    .     gmx,dhyx,zhtx,csx,muex,dpdex,keyerr)
c     lseix=eix
c     call find_interpHeos(inputz,var,temp,rtrans,
c    .     abar,zbar,eix,prx,enx,cvx,gmx,csx,muex,dpdex,keyerr)
c     heix=eix
c     ematch=lseix-heix
c
      if(rho.ge.rtrans+drtrans)then
c
c .. if rho is high enough, call the lattimer-swesty eos only
c
         call find_interpLSeos(inputz,var,temp,rho,ye,
     .    eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,
     .    gmx,dhyx,zhtx,csx,muex,dpdex,keyerr)
c
      elseif(rho.gt.rtrans-drtrans.and.rho.lt.rtrans+drtrans)then
c
c .. if rho is in transition region, call lattimer-swesty and then
c     the helmholtz eos.  then, interpolate the important quantities.
c
         call find_interpLSeos(inputz,var,temp,rtrans,ye,
     .    eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,
     .    gmx,dhyx,zhtx,csx,muex,dpdex,keyerr)
c
         abarl=awx
         zbarl=zax*awx


         prl  = prx
         eil  = eix
         enl  = enx
         tl   = temp
         gl   = gmx
         cl   = cvx
         csl  = csx
         mul  = muex
         dxl  = dpdex
c
         if(inputz.eq.2)var=var-ematch
         call find_interpHeos(inputz,var,temp,rtrans,
     .        abar,zbar,eix,prx,enx,cvx,gmx,csx,muex,dpdex,keyerr)
         if(inputz.eq.2)var=var+ematch
c
         prh  = prx
         eih  = eix + ematch
         enh  = enx
         th   = temp
         gh   = gmx
         ch   = cvx
         csh  = csx
         muh  = muex
         dxh  = dpdex
c
         prm  = 0.5d0*(prl + prh)
         eim  = 0.5d0*(eil + eih)
         enm  = 0.5d0*(enl + enh)
         ttm  = 0.5d0*(tl  + th )
         gmm  = 0.5d0*(gl  + gh )
         cvm  = 0.5d0*(cl  + ch )
         csm  = 0.5d0*(csl + csh)
         mum  = 0.5d0*(mul + muh)
         dxm  = 0.5d0*(dxl + dxh)
c
         call find_interpLSeos(inputz,var,temp,rtrans+drtrans,ye,
     .    eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,
     .    gmx,dhyx,zhtx,csx,muex,dpdex,keyerr)
c
         prl  = prx
         eil  = eix
         enl  = enx
         tl   = temp
         gl   = gmx
         cl   = cvx
         csl  = csx
         mul  = muex
         dxl  = dpdex
c
         if(inputz.eq.2)var=var-ematch
         call find_interpHeos(inputz,var,temp,rtrans-drtrans,
     .        abar,zbar,eix,prx,enx,cvx,gmx,csx,muex,dpdex,keyerr)
         if(inputz.eq.2)var=var+ematch
c
         prh  = prx
         eih  = eix + ematch
         enh  = enx
         th   = temp
         gh   = gmx
         ch   = cvx
         csh  = csx
         muh  = muex
         dxh  = dpdex
c
         r1 = rtrans-drtrans
         r2 = rtrans
         r3 = rtrans+drtrans
c
         prx   = qinter(rho,r1,prh,r2,prm,r3,prl)
         if(inputz .eq. 2) then
           eix = var
         else
           eix   = qinter(rho,r1,eih,r2,eim,r3,eil)
         endif
         enx   = qinter(rho,r1,enh,r2,enm,r3,enl)
         temp  = qinter(rho,r1,th,r2,ttm,r3,tl)
         gmx   = qinter(rho,r1,gh,r2,gmm,r3,gl)
         cvx   = qinter(rho,r1,ch,r2,cvm,r3,cl)
         csx   = qinter(rho,r1,csh,r2,csm,r3,csl)
         muex  = qinter(rho,r1,muh,r2,mum,r3,mul)
         dpdex = qinter(rho,r1,dxh,r2,dxm,r3,dxl)
c
      elseif(rho.le.rtrans-drtrans.and.
     .         rho*ye.gt.1.d-9.and.temp.gt.6.d4)then
c
c .. if rho less than rhotrans-drhotrans, just use the helmholtz eos.
c     but, if too low a temp or rho, call the ideal eos.
c
         if(inputz.eq.2)var=var-ematch
         call find_interpHeos(inputz,var,temp,rho,abar,zbar,
     .        eix,prx,enx,cvx,gmx,csx,muex,dpdex,keyerr)
         if(inputz.eq.2)var=var+ematch
         eix = eix  + ematch
c
      else
c
c .. if rho*ye less than 4.d-10 and temp less than 2d4, call ideal eos.
c
         if(inputz.eq.2)var=var-ematch
         call IDEALeos(inputz,var,temp,rho,abar,zbar,
     .        prx,eix,enx,cvx,gmx,csx,dpdex,keyerr)
         if(inputz.eq.2)var=var+ematch
         eix = eix + ematch
c
      endif   
c
      return
      end
c----------------------------------------------------------------------------
      subroutine IDEALeos(inputz,var,temp,rho,abar,zbar,
     .     ptot,etot,stot,cv,gmx,sound,dpde,keyerr)
      implicit double precision (a-z)
      logical keyerr
      integer inputz
      if(temp.le.0.d0.or.rho.le.0.d0)then
         print*,'Bad T or Rho in IDEALeos: T,Rho',temp,rho
         keyerr=.true.
         return
      endif
      temp0=temp
      if(inputz.eq.1)then
         call find_LOWeos(temp,rho,abar,zbar,
     .       ptot,etot,stot,cv,sound,dpde)
      elseif(inputz.eq.2)then
         icount=0
         do i = 1,15  
            call find_LOWeos(temp,rho,abar,zbar,
     .           ptot,etot,stot,cv,sound,dpde)
            cv     = cv
            dd     = (var-etot)/(cv)
            icount = icount + 1
            temp   = temp+dd
            if(temp.le.0.d0) temp=temp0+0.1d0*temp0
            if(abs(dd/temp).lt.1.d-8)goto 10
         enddo
 10      continue
         call find_LOWeos(temp,rho,abar,zbar,
     .        ptot,etot,stot,cv,sound,dpde)
      endif
       stot = stot/(1.381d-16*6.022d23)
       gmx  = sound*sound*rho/ptot
      return
      end
c---------------------------------------------------------------------
c                         FIND_LOWEOS
c---------------------------------------------------------------------
c
c This is the low temperature, low density EOS.  It has radiation,
c  and ions and electrons treated as ideal gases.  
c
c This code has been adapted from the Helmholtz EOS found on Frank
c  Timmes' website.
c
c IN: 
c   index     : integer, zone number, used only for debugging
c   temp      : temperature in K
c   rho       : mass density in g cm^{-3}
c   abar,zbar : 
c
c OUT: 
c   ptot  : total pressure 
c   etot  : total internal energy
c   stot  : total entropy
c   cv    : the specific heat
c   sound : the adiabatic sound speed 
c
c---------------------------------------------------------------------
      subroutine find_LOWeos(temp,rho,abar,zbar,
     .   ptot,etot,stot,cv,sound,dpde)
      implicit double precision (a-z)
      parameter(pi=3.1415926535897932384d0,amu=1.6605402d-24)
      parameter(kerg=1.380658d-16,clight=2.99792458d10,avo=6.0221367d23)
      parameter(qe=4.8032068d-10,h=6.6260755d-27,ssol=5.67051d-5)
      asol    = 4.0d0*ssol/clight
      sioncon = (2.0d0*pi*amu*kerg)/(h*h)
      forth   = 4.0d0/3.0d0
      forpi   = 4.0d0*pi
      kergavo = kerg*avo 
      ikavo   = 1.0d0/kergavo
      asoli3  = asol/3.0d0
      light2  = clight*clight
c
      den     = rho
c
      ytot1   = 1.0d0/abar
      ye      = ytot1 * zbar
      deni    = 1.0d0/den
      tempi   = 1.0d0/temp 
      kt      = kerg * temp
      ktinv   = 1.0d0/kt
c
      prad    = asoli3*temp*temp*temp*temp
      dpraddd = 0.0d0
      dpraddt = 4.0d0*prad*tempi
      dpradda = 0.0d0
      dpraddz = 0.0d0
      erad    = 3.0d0*prad*deni
      deraddd = -erad*deni
      deraddt = 3.0d0*dpraddt*deni
      deradda = 0.0d0
      deraddz = 0.0d0
      srad    = (prad*deni+erad)*tempi
      dsraddd = (dpraddd*deni-prad*deni*deni+deraddd)*tempi
      dsraddt = (dpraddt*deni+deraddt-srad)*tempi
      dsradda = 0.0d0
      dsraddz = 0.0d0
c
      xni     = avo*ytot1*den
      dxnidd  = avo*ytot1
      dxnida  = -xni*ytot1
      pion    = xni*kt
      dpiondd = dxnidd*kt
      dpiondt = xni*kerg
      dpionda = dxnida*kt 
      dpiondz = 0.0d0
      eion    = 1.5d0*pion*deni
      deiondd = (1.5d0*dpiondd-eion)*deni
      deiondt = 1.5d0*dpiondt*deni
      deionda = 1.5d0*dpionda*deni
      deiondz = 0.0d0
      x       = abar*abar*sqrt(abar)*deni/avo
      s       = sioncon*temp
      z       = x*s*sqrt(s)
      y       = log(z)
c
      sion1    = (pion*deni+eion)
      sion2    = sion1*tempi
      sion3    = sion2+kergavo*ytot1*y
      sion     = sion3
c
      dsiondd = (dpiondd*deni-pion*deni*deni+deiondd)*tempi
     .            - kergavo*deni*ytot1
      dsiondt = (dpiondt*deni+deiondt)*tempi- 
     .           (pion*deni+eion)*tempi*tempi 
     .           +1.5d0*kergavo*tempi*ytot1
      x       = avo*kerg/abar
      dsionda = (dpionda*deni+deionda)*tempi 
     .            +kergavo*ytot1*ytot1*(2.5d0-y)
      dsiondz = 0.0d0
c
      call electroneos(pi,avo,temp,tempi,den,deni,ye,kt,kerg,kergavo,
     .             h,pe,ee,se,dpedt,deedt,dpedd)
c
      stot    = srad + sion + se
      ptot    = prad + pion + pe
      etot    = erad + eion + ee
      dptotdd = dpraddd + dpiondd + dpedd 
      dptotdt = dpraddt + dpiondt + dpedt
      cv      = deraddt + deiondt + deedt
      zz      = ptot*deni
      zzi     = den/ptot
      chit    = temp/ptot * dptotdt
      chid    = dptotdd*zzi
      x       = zz * chit/(temp * cv)
      gam3    = x + 1.0d0
      gam1    = chit*x + chid
      nabad   = x/gam1
      gam2    = 1.0d0/(1.0d0 - nabad)
      cp      = cv * gam1/chid
      z       = 1.0d0 + (etot + light2)*zzi
      sound   = clight * sqrt(gam1/z)
c
      temp    = temp
      dpde=dptotdt/cv
      return
      end
c-----------------------------------------------------------------------------
      subroutine electroneos(pi,avo,temp,tempi,den,deni,ye,kt,kerg,
     .             kergavo,h,pe,ee,se,dpedt,deedt,dpedd)
      implicit none
      double precision temp,tempi,den,deni,ye,kt,kerg,kergavo,h,pe,ee,z
      double precision se,dpedt,deedt,xni,avo,dxnidde,dpedd,aaa,pi,x,s,y
      xni     = avo * ye * den
      dxnidde = avo * ye
      dpedd   = dxnidde * kt
      pe      = xni * kt
      ee      = 1.5d0 * pe * deni
      dpedt   = xni * kerg
      deedt   = 1.5d0 * dpedt*deni
      aaa     = (2.0d0 * pi * 9.11d-28 * kerg)/(h*h)
      x       = deni/avo/ye
      s       = aaa*temp
      z       = x * 2.d0 * s * sqrt(s)
      y       = log(z)
      se      = (pe*deni + ee)*tempi + kergavo * ye * y
      return
      end
c--------------------------------------------------------------------------
      subroutine find_interpLSeos(inputz,var,temp,rho,ye,
     .    eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,
     .    gmx,dhyx,zhtx,csx,muex,dpdex,keyerr)
      implicit none
      double precision temp,rho,ye,dhyx,zhtx,q,p,delty,qinter,energy
      double precision rl,ql,alpha,beta,delta,linter,q0,q1,q2,qq
      double precision dhy0,zht0,dhy1,zht1,dhy2,zht2,yl0,yl1,yl2,csx
      double precision eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,gmx
      double precision entropy,f,var,mevcon,muex,dpdex, rhloc
      parameter (mevcon=0.95655684d18)
      integer i,j,k,jr,jq,jy,nn,inputz
      logical keyerr
      common /quantity/ qq(50)
      include 'table_all_params.ecomp-2'
      temp   = temp/1.160445d10
c      if(rho.lt.2.5d6.or.rho.gt.1.26d15.or.
c     .   temp.lt.0.16d0.or.temp.gt.39.d0.or.
c     .   ye.lt.0.05d0.or.ye.gt.0.51d0)then
c         print*,'Outside LSeos limits in find_interpLSeos'
c         print*,'T, Rho, Ye: ',temp,rho,ye
c         temp=temp*1.160445d10
c         keyerr=.true.
c         return
c      endif
      if(ye.lt.0.052d0)ye=0.05d0
c
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
         call findtemp(inputz,energy,temp,rho,ye,keyerr)
         do nn=1,numel
            call findthis(nn,temp,rho,ye,f)
            qq(nn)=f
         enddo
         eix = energy
         enx = qq(3)
      else if(inputz.eq.3)then
         entropy=var
         call findtemp(inputz,entropy,temp,rho,ye,keyerr)
         do nn=1,numel
            call findthis(nn,temp,rho,ye,f)
            qq(nn)=f
         enddo
         enx=entropy
         eix=qq(1)
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
      rhloc = rho*(1.d0+eix) + prx
      csx  = sqrt(prx/rhloc*gmx)
c
      if(keyerr)then
         print*,'Problem in Findtemp: keyerr',keyerr
         return
      endif

      return
      end


c     Branson messing around.  A version with location information.
      subroutine find_interpLSeos_loc(inputz,var,temp,rho,ye,
     .    eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,
     .    gmx,dhyx,zhtx,csx,muex,dpdex,keyerr,il,jl,kl)
      implicit none
      double precision temp,rho,ye,dhyx,zhtx,q,p,delty,qinter,energy
      double precision rl,ql,alpha,beta,delta,linter,q0,q1,q2,qq
      double precision dhy0,zht0,dhy1,zht1,dhy2,zht2,yl0,yl1,yl2,csx
      double precision eix,prx,enx,cvx,xnx,xpx,xax,xhx,zax,awx,mhx,gmx
      double precision entropy,f,var,mevcon,muex,dpdex,rhloc,c
      parameter (mevcon=0.95655684d18)
      parameter (c = 2.99792458d10)
      integer i,j,k,jr,jq,jy,nn,inputz,il,jl,kl
      logical keyerr
      common /quantity/ qq(50)
      include 'table_all_params.ecomp-2'
      temp   = temp/1.160445d10
c      if(rho.lt.2.5d6.or.rho.gt.1.26d15.or.
c     .   temp.lt.0.16d0.or.temp.gt.39.d0.or.
c     .   ye.lt.0.05d0.or.ye.gt.0.51d0)then
c         print*,'Outside LSeos limits in find_interpLSeos'
c         print*,'T, Rho, Ye: ',temp,rho,ye
c         temp=temp*1.160445d10
c         keyerr=.true.
c         return
c      endif
      if(ye.lt.0.052d0)ye=0.05d0
c
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
         call findtemp_loc(inputz,energy,temp,rho,ye,keyerr,il,jl,kl)
         do nn=1,numel
            call findthis(nn,temp,rho,ye,f)
            qq(nn)=f
         enddo
         eix = energy
         enx = qq(3)
      else if(inputz.eq.3)then
         entropy=var
         call findtemp(inputz,entropy,temp,rho,ye,keyerr)
         do nn=1,numel
            call findthis(nn,temp,rho,ye,f)
            qq(nn)=f
         enddo
         enx=entropy
         eix=qq(1)
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
      return
      end
c----------------------------------------------------------------------
      subroutine findtemp(inputz,var,temp,rho,ye,keyerr)
      implicit none
      double precision ep,temp,rho,ye,e,ep1,ep2,dt,depdt,dd,f,var
      double precision q,p,delty,qinter,energy,yl0,yl1,yl2,told
      double precision rl,ql,alpha,beta,delta,linter,q0,q1,q2,qq
      double precision temp0,xxx
      integer i,jy,jr,jq,nn,icount,inputz,index,nnnstep
      logical keyerr,gobisect
      common /indices/ jy,jr,jq
      common /xxpass/ index,nnnstep,icount
c     include 'table_all_params.shen'
      include 'table_all_params.ecomp-2'
      
      temp0=temp
      if(inputz.eq.2)nn=1
      if(inputz.eq.3)nn=3
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
         if(abs(dd/temp).lt.1.d-8)goto 10
      enddo
 10   continue
      if(abs(temp-temp0)/temp0.gt.0.5d0.or.icount.eq.15.or.gobisect)then
         call bisection(inputz,var,temp0,rho,ye)
         temp=temp0
      endif   
      return
      end

c     Branson messing around.  A version with location information.
      subroutine findtemp_loc(inputz,var,temp,rho,ye,keyerr,il,jl,kl)
      implicit none
      double precision ep,temp,rho,ye,e,ep1,ep2,dt,depdt,dd,f,var
      double precision q,p,delty,qinter,energy,yl0,yl1,yl2,told
      double precision rl,ql,alpha,beta,delta,linter,q0,q1,q2,qq
      double precision temp0,xxx
      integer i,jy,jr,jq,nn,icount,inputz,index,nnnstep,il,jl,kl
      logical keyerr,gobisect
      common /indices/ jy,jr,jq
      common /xxpass/ index,nnnstep,icount
c     include 'table_all_params.shen'
      include 'table_all_params.ecomp-2'
      
      temp0=temp
      if(inputz.eq.2)nn=1
      if(inputz.eq.3)nn=3
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
         if(abs(dd/temp).lt.1.d-8)goto 10
      enddo
 10   continue
      if(abs(temp-temp0)/temp0.gt.0.5d0.or.icount.eq.15.or.gobisect)then
         call bisection_loc(inputz,var,temp0,rho,ye,il,jl,kl)
         temp=temp0
      endif   
      return
      end
c----------------------------------------------------------------------
      subroutine bisection(inputz,var,temp,rho,ye)
      implicit none
      double precision ep,temp,rho,ye,e,ep1,ep2,dt,depdt,dd,f,var
      double precision q,p,delty,qinter,energy,yl0,yl1,yl2,told
      double precision rl,ql,alpha,beta,delta,linter,q0,q1,q2,qq
      double precision temp1,temp2,f1,f2,dtemp,tmid,fmid,temp00
      integer i,jy,jr,jq,nn,icount,inputz,index,nnnstep,mcount
      common /indices/ jy,jr,jq
      common /xxpass/ index,nnnstep,icount
c     include 'table_all_params.shen'
      include 'table_all_params.ecomp-2'
      if(inputz.eq.2)nn=1
      if(inputz.eq.3)nn=3
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
            tmid=temp00
            write(*,*) 'findtemp giving up'
            stop
c            goto 10
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
c      tmid=temp00
      temp=temp00
 10   continue
      temp=temp
      return
      end

c     Branson messing around.  A version with location information.
      subroutine bisection_loc(inputz,var,temp,rho,ye,il,jl,kl)
      implicit none
      double precision ep,temp,rho,ye,e,ep1,ep2,dt,depdt,dd,f,var
      double precision q,p,delty,qinter,energy,yl0,yl1,yl2,told
      double precision rl,ql,alpha,beta,delta,linter,q0,q1,q2,qq
      double precision temp1,temp2,f1,f2,dtemp,tmid,fmid,temp00
      integer i,jy,jr,jq,nn,icount,inputz,index,nnnstep,mcount
      integer il,jl,kl
      common /indices/ jy,jr,jq
      common /xxpass/ index,nnnstep,icount
c     include 'table_all_params.shen'
      include 'table_all_params.ecomp-2'
      if(inputz.eq.2)nn=1
      if(inputz.eq.3)nn=3
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
            tmid=temp00
            write(*,*) 'findtemp giving up at ', il, ' ', jl, ' ', kl
            stop
c            goto 10
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
c      tmid=temp00
      temp=temp00
 10   continue
      temp=temp
      return
      end
c----------------------------------------------------------------------
      subroutine offtable
      implicit none
      integer jy,jr,jq
      common /indices/ jy,jr,jq
c     include 'table_all_params.shen'
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
c     include 'table_all_params.shen'
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
c
c-----------------------------------------------------------------------
c          THE HELMHOLTZ EOS FROM TIMMES (EVERYTHING BELOW)
c-----------------------------------------------------------------------
c
      subroutine find_interpHeos(inputz,var,t,r,abar,zbar,
     .        eix,prx,enx,cvx,gmx,csx,muex,dpdex,keyerr)
      implicit none
      save
      include 'vector_eos.dek'
      double precision var,t,r,y,xxm,xxa,xxz,eix,prx,enx,cvx,xnx,xpx,xax
      double precision xhx,zax,awx,mhx,gmx,dhyx,zhtx,temp0,muex,dpdex
      double precision dt,e1,e2,dedt,dd,csx
      integer index,inputz,i,ionmax,icount,nnnstep
      double precision temp,den,abar,zbar
      logical keyerr
      temp_row(1) = t
      den_row(1)  = r
      abar_row(1) = abar
      zbar_row(1) = zbar
      jlo_eos = 1
      jhi_eos = 1
 11   format(5e14.5)
      if(inputz.eq.1)then
         call helmeos
         if(eosfail)then
            print*,'Problem in Heos, returning Keyerr'
            write(6,11)var,t,r,abar,zbar
            keyerr=.true.
            return
         endif
         t     = temp_row(1)
         eix   = etot_row(1)
         prx   = ptot_row(1)
         enx   = stot_row(1)/(1.381d-16*6.022d23)
         gmx   = cs_row(1)*cs_row(1)*r/prx
         cvx   = cv_row(1)
         csx   = cs_row(1)
         muex  = etaele_row(1)*t/1.160445d10
         dpdex = dpt_row(1)/det_row(1)
      elseif(inputz.eq.2)then
         temp0  = temp_row(1)
         icount = 0
         do i=1,20
            call helmeos
            if(eosfail)then
               print*,'Problem in Heos, returning Keyerr'
               write(6,11)var,t,r,abar,zbar
               keyerr=.true.
               return
            endif
            dd          = (var-etot_row(1))/cv_row(1)
            icount      = icount+1
            temp_row(1) = temp_row(1)+dd
            if(temp_row(1).le.0.d0) temp_row(1) = temp0+0.1d0*temp0
            if(abs(dd/temp_row(1)).lt.1.d-8)goto 10
         enddo
 10      continue
         if(abs(temp_row(1)-temp0)/temp0.gt.1.d0.or.icount.eq.20)then
            call helm_bisection(inputz,var,temp0,den_row(1),keyerr)
            temp_row(1)=temp0
         endif   
         call helmeos
         t     = temp_row(1)
         eix   = etot_row(1)
         prx   = ptot_row(1)
         enx   = stot_row(1)/(1.381d-16*6.022d23)
         gmx   = cs_row(1)*cs_row(1)*r/ptot_row(1)
         cvx   = cv_row(1)
         csx   = cs_row(1)
         muex  = etaele_row(1)*t/1.160445d10
         dpdex = dpt_row(1)/det_row(1)
c
         if(eosfail)then
            keyerr=.true.
         endif
c
      endif
      return
      end 
c----------------------------------------------------------------------
      subroutine helm_bisection(inputz,var,temp,rho,keyerr)
      implicit none
      double precision ep,temp,rho,ye,e,ep1,ep2,dt,depdt,dd,f,var
      double precision q,p,delty,qinter,energy,yl0,yl1,yl2,told
      double precision rl,ql,alpha,beta,delta,linter,q0,q1,q2,qq
      double precision temp1,temp2,f1,f2,dtemp,tmid,fmid,temp00
      integer i,jy,jr,jq,nn,icount,inputz,index,nnnstep,mcount
      logical keyerr
      include 'vector_eos.dek'
      den_row(1)  = rho
      temp00      = temp
      mcount      = 0
      temp1       = temp-0.1d0*temp
      temp2       = temp+0.1d0*temp
      temp_row(1) = temp
      call helmeos
      f = etot_row(1)
      f = f-var
      temp_row(1) = temp1
      call helmeos
      f1 = etot_row(1)
      temp_row(1) = temp2
      call helmeos
      f2 = etot_row(1)
      f1 = f1-var
      f2 = f2-var
 2    continue
      if(f1*f2.ge.0.d0)then
         mcount = mcount+1
         temp1  = temp1/1.2d0
         temp2  = temp2*1.2d0
         temp_row(1) = temp1
         call helmeos
         f1 = etot_row(1)
         temp_row(1) = temp2
         call helmeos
         f2 = etot_row(1)
         f1 = f1-var
         f2 = f2-var
         if(mcount.le.9)goto 2
         if(mcount.gt.9)then
            tmid = temp00
            goto 10
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
         temp_row(1)=tmid
         call helmeos
         if(eosfail)then
            print*,'problem in helmbisection, returning temp'
            print*,temp,tmid,var
            temp=temp00
            keyerr=.true.
            return
         endif
         fmid=etot_row(1)-var
         if(fmid.le.0.d0)temp=tmid
         if(abs(dtemp).lt.1.d-6)goto 10
      enddo
c      tmid=temp00
      temp=temp00
 10   continue
      temp=temp
      return
      end
c---------------------------------------------------------------------
      subroutine azbar(xmass,aion,zion,ionmax,
     1                 ymass,abar,zbar)
      implicit none
      save
c..this routine calculates composition variables for an eos routine
c..input:
c..mass fractions     = xmass(1:ionmax)
c..number of nucleons = aion(1:ionmax)
c..charge of nucleus  = zion(1:ionmax)
c..number of isotopes = ionmax
c..output:
c..molar abundances        = ymass(1:ionmax), 
c..mean number of nucleons = abar
c..mean nucleon charge     = zbar
c..declare
      integer          i,ionmax
      double precision xmass(ionmax),aion(ionmax),zion(ionmax),
     1                 ymass(ionmax),abar,zbar,zbarxx,ytot1
      zbarxx  = 0.0d0
      ytot1   = 0.0d0
      do i=1,ionmax
       ymass(i) = xmass(i)/aion(i)
       ytot1    = ytot1 + ymass(i)
       zbarxx   = zbarxx + zion(i) * ymass(i)
      enddo
      abar   = 1.0d0/ytot1
      zbar   = zbarxx * abar
      return
      end
c------------------------------------------------------------------
      subroutine pretty_eos_out(whose)
      implicit none
      save
      include 'vector_eos.dek'
c..
c..writes a pretty output for the eos tester
c..
c..declare
      integer     j
      character*7 whose

c..popular formats
01    format(1x,t2,a,t11,'total',t24,'ion',t34,'e- & e+',
     1       t46,'radiation',t58,'coulomb')
02    format(1x,t2,a,1p6e12.4)
03    format(1x,t2,a6,1pe12.4,t22,a6,1pe12.4,
     1         t42,a6,1pe12.4,t62,a6,1pe12.4)

      do j=jlo_eos,jhi_eos

c..the input 
      write(6,03) 'temp =',temp_row(1),'den  =',den_row(1),
     1            'abar =',abar_row(1),'zbar =',zbar_row(1)
      write(6,*) ' ' 
c..and the output
c..first the totals from each of the components
      write(6,01)  whose
      write(6,02) 'pres =',
     1            ptot_row(j),pion_row(j),pele_row(j),
     2            prad_row(j),pcou_row(j)
      write(6,02) 'ener =',
     1            etot_row(j),eion_row(j),eele_row(j),
     2            erad_row(j),ecou_row(j)
      write(6,02) 'entr =',
     1            stot_row(j),sion_row(j),sele_row(j),
     2            srad_row(j),scou_row(j)
c..derivatives of the totals with respect to the input variables
      write(6,*)  ' '
      write(6,03) 'dp/dd=',dpd_row(j),'dp/dt=',dpt_row(j),
     1            'dp/da=',dpa_row(j),'dp/dz=',dpz_row(j)
      write(6,03) 'de/dd=',ded_row(j),'de/dt=',det_row(j),
     1            'de/da=',dea_row(j),'de/dz=',dez_row(j)
      write(6,03) 'ds/dd=',dsd_row(j),'ds/dt=',dst_row(j),
     1            'ds/da=',dsa_row(j),'ds/dz=',dsz_row(j)
c..derivatives of the electron-positron compoenets with
c..respect to the input variables
      write(6,*) ' ' 
      write(6,03) 'dpepd=',dpepd_row(j),'dpept=',dpept_row(j),
     1            'dpepa=',dpepa_row(j),'dpepz=',dpepz_row(j)
      write(6,03) 'deepd=',deepd_row(j),'deept=',deept_row(j),
     1            'deepa=',deepa_row(j),'deepz=',deepz_row(j)
      write(6,03) 'dsepd=',dsepd_row(j),'dsept=',dsept_row(j),
     1            'dsepa=',dsepa_row(j),'dsepz=',dsepz_row(j)
c..the thermodynamic consistency relations, these should all be
c..at the floating poiint limit of zero
      write(6,*) ' ' 
      write(6,03) 'maxw1=',dse_row(j),'maxw2=',dpe_row(j),
     1            'maxw3=',dsp_row(j)
c..number density of electrons, poistrons, matter electrons, and ions
      write(6,03) 'xne  =',xne_row(j),'xnp  =',xnp_row(j),
     1            'xnem =',xnem_row(j),'xni  =',xni_row(j)
c..derivatibves of the electron number density with 
c..respect to the input variables
      write(6,03) 'dxned=',dxned_row(j),'dxnet=',dxnet_row(j),
     1            'dxnea=',dxnea_row(j),'dxnez=',dxnez_row(j)
c..electron chemical potential, positron chemical potential
c..and derivatives of electron chemical potential with respect
c..to the input variables
      write(6,03) 'eta  =',etaele_row(j),'etap =',etapos_row(j)
      write(6,03) 'detad=',detad_row(j),'detat=',detat_row(j),
     1            'detaa=',detaa_row(j),'detaz=',detaz_row(j)
c..specific heats, and ratio of electostatic to thermal energy
      write(6,03) 'cp   =',cp_row(j),'cv   =',cv_row(j),
     1            'plasg=',plasg_row(j)
c..the 3 gammas and the sound speed
      write(6,03) 'gam1 =',gam1_row(j),'gam2 =',gam2_row(j),
     1            'gam3 =',gam3_row(j),'csond=',cs_row(j)
      write(6,*) ' '
      enddo
      return
      end
c------------------------------------------------------------------------
      subroutine helmeos
      implicit none
      save
      include 'vector_eos.dek'
c..given a temperature temp [K], density den [g/cm**3], and a composition 
c..characterized by abar and zbar, this routine returns most of the other 
c..thermodynamic quantities. of prime interest is the pressure [erg/cm**3], 
c..specific thermal energy [erg/gr], the entropy [erg/g/K], along with 
c..their derivatives with respect to temperature, density, abar, and zbar.
c..other quantites such the normalized chemical potential eta (plus its
c..derivatives), number density of electrons and positron pair (along 
c..with their derivatives), adiabatic indices, specific heats, and 
c..relativistically correct sound speed are also returned.
c..
c..this routine assumes planckian photons, an ideal gas of ions,
c..and an electron-positron gas with an arbitrary degree of relativity
c..and degeneracy. interpolation in a table of the helmholtz free energy
c..is used to return the electron-positron thermodynamic quantities.
c..all other derivatives are analytic.
c..
c..references: cox & giuli chapter 24 ; timmes & swesty apj 1999
c..declare
      double precision pi,amu,kerg,clight,avo,qe,h,ssol,asol
      parameter       (pi      = 3.1415926535897932384d0,
     1                  amu    = 1.6605402d-24,
     2                  kerg   = 1.380658d-16,
     3                  clight = 2.99792458d10, 
     4                  avo    = 6.0221367d23,
     5                  qe     = 4.8032068d-10,  
     6                  h      = 6.6260755d-27,
     7                  ssol   = 5.67051d-5,
     8                  asol   = 4.0d0 * ssol / clight)
      integer          i,j
      double precision x,y,zz,zzi,deni,tempi,xni,dxnidd,dxnida,
     1                 dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt,
     2                 dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt,
     3                 deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt,
     4                 dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion,
     5                 sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd,
     6                 dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp,
     7                 gam1,gam2,gam3,chit,chid,nabad,sound,etaele,
     8                 detadt,detadd,xnefer,dxnedt,dxnedd,s,
     9                 temp,den,abar,zbar,ytot1,ye,
     &                 sioncon,forth,forpi,kergavo,ikavo,asoli3,light2
      parameter        (sioncon = (2.0d0 * pi * amu * kerg)/(h*h),
     1                  forth   = 4.0d0/3.0d0,
     2                  forpi   = 4.0d0 * pi,
     3                  kergavo = kerg * avo, 
     4                  ikavo   = 1.0d0/kergavo,
     5                  asoli3  = asol/3.0d0,
     6                  light2  = clight * clight)
c..for the abar derivatives
      double precision dpradda,deradda,dsradda,
     1                 dpionda,deionda,dsionda,
     2                 dpepda,deepda,dsepda,
     3                 dpresda,denerda,dentrda,
     4                 detada,dxneda
c..for the zbar derivatives
      double precision dpraddz,deraddz,dsraddz,
     1                 dpiondz,deiondz,dsiondz,
     2                 dpepdz,deepdz,dsepdz,
     3                 dpresdz,denerdz,dentrdz,
     4                 detadz,dxnedz
c..for the tables, in general
      integer          imax,jmax
      parameter        (imax = 211, jmax = 71)
      double precision d(imax),t(jmax)
c..for the helmholtz free energy tables
      double precision f(imax,jmax),fd(imax,jmax),
     1                 ft(imax,jmax),fdd(imax,jmax),ftt(imax,jmax),
     2                 fdt(imax,jmax),fddt(imax,jmax),fdtt(imax,jmax),
     3                 fddtt(imax,jmax)
c..for the pressure derivative with density ables
      double precision dpdf(imax,jmax),dpdfd(imax,jmax),
     1                 dpdft(imax,jmax),dpdfdd(imax,jmax),
     2                 dpdftt(imax,jmax),dpdfdt(imax,jmax)
c..for chemical potential tables
      double precision ef(imax,jmax),efd(imax,jmax),
     1                 eft(imax,jmax),efdd(imax,jmax),eftt(imax,jmax),
     2                 efdt(imax,jmax)
c..for the number density tables
      double precision xf(imax,jmax),xfd(imax,jmax),
     1                 xft(imax,jmax),xfdd(imax,jmax),xftt(imax,jmax),
     2                 xfdt(imax,jmax)
c..for the interpolations
      integer          iat,jat
      double precision tlo,thi,tstp,tstpi,dlo,dhi,dstp,dstpi,
     1                 tsav,dsav,free,df_d,df_t,df_dd,df_tt,df_dt
      double precision dth,dt2,dti,dt2i,dd,dd2,ddi,dd2i,xt,xd,mxt,mxd,
     1                 si0t,si1t,si2t,si0mt,si1mt,si2mt,
     2                 si0d,si1d,si2d,si0md,si1md,si2md,
     3                 dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt,
     4                 dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md,
     5                 ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt,
     6                 ddsi0d,ddsi1d,ddsi2d,ddsi0md,ddsi1md,ddsi2md,
     7                 z,psi0,dpsi0,ddpsi0,psi1,dpsi1,ddpsi1,psi2,
     8                 dpsi2,ddpsi2,din,h5,fi(36),
     9                 xpsi0,xdpsi0,xpsi1,xdpsi1,h3,
     1                 w0t,w1t,w2t,w0mt,w1mt,w2mt,
     2                 w0d,w1d,w2d,w0md,w1md,w2md
c..for storing the differences
      double precision dt_sav(jmax),dt2_sav(jmax),
     1                 dti_sav(jmax),dt2i_sav(jmax),
     2                 dd_sav(imax),dd2_sav(imax),
     3                 ddi_sav(imax),dd2i_sav(imax)
c..for the coulomb corrections
      double precision dsdd,dsda,lami,inv_lami,lamida,lamidd,
     1                 plasg,plasgdd,plasgdt,plasgda,plasgdz,
     1                 a1,b1,c1,d1,e1,a2,b2,c2,
     3                 ecoul,decouldd,decouldt,decoulda,decouldz,
     4                 pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz,
     5                 scoul,dscouldd,dscouldt,dscoulda,dscouldz,
     6                 tmelt,tfermi,rhocond,z2,x1,x2,third,esqu
      parameter        (a1    = -0.898004d0, 
     1                  b1    =  0.96786d0, 
     2                  c1    =  0.220703d0, 
     3                  d1    = -0.86097d0,
     4                  e1    =  2.5269d0, 
     5                  a2    =  0.29561d0, 
     6                  b2    =  1.9885d0,    
     7                  c2    =  0.288675d0,
     8                  third = 1.0d0/3.0d0,
     9                  esqu  = qe * qe)
c..for initialization
      integer          ifirst
      data             ifirst/0/ 
c..quintic hermite polynomial statement functions
c..psi0 and its derivatives
      psi0(z)   = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
      dpsi0(z)  = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
      ddpsi0(z) = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)
c..psi1 and its derivatives
      psi1(z)   = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
      dpsi1(z)  = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
      ddpsi1(z) = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)
c..psi2  and its derivatives
      psi2(z)   = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
      dpsi2(z)  = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
      ddpsi2(z) = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)
c..biquintic hermite polynomial statement function
      h5(i,j,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)=
     1       fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t
     2     + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt
     3     + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t
     4     + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt
     5     + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t
     6     + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt
     7     + fi(13) *w1d*w0t   + fi(14) *w1md*w0t
     8     + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt
     9     + fi(17) *w2d*w0t   + fi(18) *w2md*w0t
     &     + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt
     1     + fi(21) *w1d*w1t   + fi(22) *w1md*w1t
     2     + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt
     3     + fi(25) *w2d*w1t   + fi(26) *w2md*w1t
     4     + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt
     5     + fi(29) *w1d*w2t   + fi(30) *w1md*w2t
     6     + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt
     7     + fi(33) *w2d*w2t   + fi(34) *w2md*w2t
     8     + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt
c..cubic hermite polynomial statement functions
c..psi0 & derivatives
      xpsi0(z)  = z * z * (2.0d0*z - 3.0d0) + 1.0
      xdpsi0(z) = z * (6.0d0*z - 6.0d0)
c..psi1 & derivatives
      xpsi1(z)  = z * ( z * (z - 2.0d0) + 1.0d0)
      xdpsi1(z) = z * (3.0d0*z - 4.0d0) + 1.0d0
c..bicubic hermite polynomial statement function
      h3(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) = 
     1       fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t 
     2     + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt
     3     + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t 
     4     + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt
     5     + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t 
     6     + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt
     7     + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t 
     8     + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt
c..popular format statements
01    format(1x,5(a,1pe11.3))
02    format(1x,a,1p4e16.8)
03    format(1x,4(a,1pe11.3))
04    format(1x,4(a,i4))
c..do this stuff once
      if (ifirst .eq. 0) then
       ifirst = 1
c..open the table
       write(6,*)
       print*,'........... Initializing HELMHOLTZ EPtable ...........'
      open(unit=2,file='../htables/helm_table.dat',
     .             status='old')
c..read the helmholtz free energy table
       tlo   = 4.0d0
       thi   = 11.0d0
       tstp  = (thi - tlo)/float(jmax-1)
       tstpi = 1.0d0/tstp
       dlo   = -10.0d0
       dhi   = 11.0d0
       dstp  = (dhi - dlo)/float(imax-1)
       dstpi = 1.0d0/dstp
       do j=1,jmax
        tsav = tlo + (j-1)*tstp
        t(j) = 10.0d0**(tsav)
        do i=1,imax
         dsav = dlo + (i-1)*dstp
         d(i) = 10.0d0**(dsav)
         read(2,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j),
     1            fddt(i,j),fdtt(i,j),fddtt(i,j)
        enddo
       enddo
c..read the pressure derivative with density table
       do j=1,jmax
        do i=1,imax
         read(2,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
        enddo
       enddo
c..read the electron chemical potential table
       do j=1,jmax
        do i=1,imax
         read(2,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
        enddo
       enddo
c..read the number density table
       do j=1,jmax
        do i=1,imax
         read(2,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
        enddo
       enddo
c..construct the temperature and density deltas and their inverses 
       do j=1,jmax-1
        dth          = t(j+1) - t(j)
        dt2         = dth * dth
        dti         = 1.0d0/dth
        dt2i        = 1.0d0/dt2
        dt_sav(j)   = dth
        dt2_sav(j)  = dt2
        dti_sav(j)  = dti
        dt2i_sav(j) = dt2i
       end do
       do i=1,imax-1
        dd          = d(i+1) - d(i)
        dd2         = dd * dd
        ddi         = 1.0d0/dd
        dd2i        = 1.0d0/dd2
        dd_sav(i)   = dd
        dd2_sav(i)  = dd2
        ddi_sav(i)  = ddi
        dd2i_sav(i) = dd2i
       enddo
       close(unit=2)
c       write(6,*)
c       write(6,*) 'finished reading eos table'
       print*,'...... Finished Helmholtz EPtable Initialization .....'
c       write(6,04) 'imax=',imax,' jmax=',jmax
c       write(6,03) 'temp(1)   =',t(1),' temp(jmax)   =',t(jmax)
c       write(6,03) 'ye*den(1) =',d(1),' ye*den(imax) =',d(imax)
       write(6,*)
      end if
c..start of vectorization loop, normal executaion starts here
      eosfail = .false.
      do j=jlo_eos,jhi_eos
       if (temp_row(j) .le. 0.0) stop 'temp less than 0 in helmeos'
       if (den_row(j)  .le. 0.0) stop 'den less than 0 in helmeos'
       temp  = temp_row(j)
       den   = den_row(j)
       abar  = abar_row(j)
       zbar  = zbar_row(j)
       ytot1 = 1.0d0/abar
       ye    = ytot1 * zbar
c..initialize
       deni    = 1.0d0/den
       tempi   = 1.0d0/temp 
       kt      = kerg * temp
       ktinv   = 1.0d0/kt
c..radiation section:
       prad    = asoli3 * temp * temp * temp * temp
       dpraddd = 0.0d0
       dpraddt = 4.0d0 * prad*tempi
       dpradda = 0.0d0
       dpraddz = 0.0d0
       erad    = 3.0d0 * prad*deni
       deraddd = -erad*deni
       deraddt = 3.0d0 * dpraddt*deni
       deradda = 0.0d0
       deraddz = 0.0d0
       srad    = (prad*deni + erad)*tempi
       dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
       dsraddt = (dpraddt*deni + deraddt - srad)*tempi
       dsradda = 0.0d0
       dsraddz = 0.0d0
c..ion section:
       xni     = avo * ytot1 * den
       dxnidd  = avo * ytot1
       dxnida  = -xni * ytot1
       pion    = xni * kt
       dpiondd = dxnidd * kt
       dpiondt = xni * kerg
       dpionda = dxnida * kt 
       dpiondz = 0.0d0
       eion    = 1.5d0 * pion*deni
       deiondd = (1.5d0 * dpiondd - eion)*deni
       deiondt = 1.5d0 * dpiondt*deni
       deionda = 1.5d0 * dpionda*deni
       deiondz = 0.0d0
       x       = abar*abar*sqrt(abar) * deni/avo
       s       = sioncon * temp
       z       = x * s * sqrt(s)
       y       = log(z)
       sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
       dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi
     1            - kergavo * deni * ytot1
       dsiondt = (dpiondt*deni + deiondt)*tempi - 
     1           (pion*deni + eion) * tempi*tempi 
     2           + 1.5d0 * kergavo * tempi*ytot1
       x       = avo*kerg/abar
       dsionda = (dpionda*deni + deionda)*tempi 
     1           + kergavo*ytot1*ytot1* (2.5d0 - y)
       dsiondz = 0.0d0
c..electron-positron section:
c..assume complete ionization 
       xnem    = xni * zbar
c..enter the table with ye*den
       din = ye*den
c..bomb proof the input
       if (temp .gt. t(jmax)) then
        write(6,01) 'temp=',temp,' t(jmax)=',t(jmax)
        write(6,*) 'temp too hot, off grid'       
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if
       if (temp .lt. t(1)) then
        write(6,01) 'temp=',temp,' t(1)=',t(1)
        write(6,*) 'temp too cold, off grid'
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if
       if (din  .gt. d(imax)) then
        write(6,01) 'den*ye=',din,' d(imax)=',d(imax)
        write(6,*) 'ye*den too big, off grid'
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if
       if (din  .lt. d(1)) then
        write(6,01) 'ye*den=',din,' d(1)=',d(1)
        write(6,*) 'ye*den too small, off grid'
        write(6,*) 'setting eosfail to true and returning'
        call flush(6)
        eosfail = .true.
        return
       end if
c..hash locate this temperature and density
       jat = int((log10(temp) - tlo)*tstpi) + 1
       jat = max(1,min(jat,jmax-1))
       iat = int((log10(din) - dlo)*dstpi) + 1
       iat = max(1,min(iat,imax-1))
c..access the table locations only once
       fi(1)  = f(iat,jat)
       fi(2)  = f(iat+1,jat)
       fi(3)  = f(iat,jat+1)
       fi(4)  = f(iat+1,jat+1)
       fi(5)  = ft(iat,jat)
       fi(6)  = ft(iat+1,jat)
       fi(7)  = ft(iat,jat+1)
       fi(8)  = ft(iat+1,jat+1)
       fi(9)  = ftt(iat,jat)
       fi(10) = ftt(iat+1,jat)
       fi(11) = ftt(iat,jat+1)
       fi(12) = ftt(iat+1,jat+1)
       fi(13) = fd(iat,jat)
       fi(14) = fd(iat+1,jat)
       fi(15) = fd(iat,jat+1)
       fi(16) = fd(iat+1,jat+1)
       fi(17) = fdd(iat,jat)
       fi(18) = fdd(iat+1,jat)
       fi(19) = fdd(iat,jat+1)
       fi(20) = fdd(iat+1,jat+1)
       fi(21) = fdt(iat,jat)
       fi(22) = fdt(iat+1,jat)
       fi(23) = fdt(iat,jat+1)
       fi(24) = fdt(iat+1,jat+1)
       fi(25) = fddt(iat,jat)
       fi(26) = fddt(iat+1,jat)
       fi(27) = fddt(iat,jat+1)
       fi(28) = fddt(iat+1,jat+1)
       fi(29) = fdtt(iat,jat)
       fi(30) = fdtt(iat+1,jat)
       fi(31) = fdtt(iat,jat+1)
       fi(32) = fdtt(iat+1,jat+1)
       fi(33) = fddtt(iat,jat)
       fi(34) = fddtt(iat+1,jat)
       fi(35) = fddtt(iat,jat+1)
       fi(36) = fddtt(iat+1,jat+1)
c..various differences
       xt  = max( (temp - t(jat))*dti_sav(jat), 0.0d0)
       xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
       mxt = 1.0d0 - xt
       mxd = 1.0d0 - xd
c..the six density and six temperature basis functions
       si0t =   psi0(xt)
       si1t =   psi1(xt)*dt_sav(jat)
       si2t =   psi2(xt)*dt2_sav(jat)
       si0mt =  psi0(mxt)
       si1mt = -psi1(mxt)*dt_sav(jat)
       si2mt =  psi2(mxt)*dt2_sav(jat)
       si0d =   psi0(xd)
       si1d =   psi1(xd)*dd_sav(iat)
       si2d =   psi2(xd)*dd2_sav(iat)
       si0md =  psi0(mxd)
       si1md = -psi1(mxd)*dd_sav(iat)
       si2md =  psi2(mxd)*dd2_sav(iat)
c..derivatives of the weight functions
       dsi0t =   dpsi0(xt)*dti_sav(jat)
       dsi1t =   dpsi1(xt)
       dsi2t =   dpsi2(xt)*dt_sav(jat)
       dsi0mt = -dpsi0(mxt)*dti_sav(jat)
       dsi1mt =  dpsi1(mxt)
       dsi2mt = -dpsi2(mxt)*dt_sav(jat)
       dsi0d =   dpsi0(xd)*ddi_sav(iat)
       dsi1d =   dpsi1(xd)
       dsi2d =   dpsi2(xd)*dd_sav(iat)
       dsi0md = -dpsi0(mxd)*ddi_sav(iat)
       dsi1md =  dpsi1(mxd)
       dsi2md = -dpsi2(mxd)*dd_sav(iat)
c..second derivatives of the weight functions
       ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
       ddsi1t =   ddpsi1(xt)*dti_sav(jat)
       ddsi2t =   ddpsi2(xt)
       ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
       ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
       ddsi2mt =  ddpsi2(mxt)
c       ddsi0d =   ddpsi0(xd)*dd2i_sav(iat)
c       ddsi1d =   ddpsi1(xd)*ddi_sav(iat)
c       ddsi2d =   ddpsi2(xd)
c       ddsi0md =  ddpsi0(mxd)*dd2i_sav(iat)
c       ddsi1md = -ddpsi1(mxd)*ddi_sav(iat)
c       ddsi2md =  ddpsi2(mxd)
c..the free energy
       free  = h5(iat,jat,
     1         si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
     2         si0d,   si1d,   si2d,   si0md,   si1md,   si2md)
c..derivative with respect to density
       df_d  = h5(iat,jat,
     1         si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
     2         dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)
c..derivative with respect to temperature
       df_t = h5(iat,jat,
     1         dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt,
     2         si0d,   si1d,   si2d,   si0md,   si1md,   si2md)
c..derivative with respect to density**2
c       df_dd = h5(iat,jat,
c     1         si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
c     2         ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)
c..derivative with respect to temperature**2
       df_tt = h5(iat,jat,
     1       ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt,
     2         si0d,   si1d,   si2d,   si0md,   si1md,   si2md)
c..derivative with respect to temperature and density
       df_dt = h5(iat,jat,
     1         dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt,
     2         dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)
c..now get the pressure derivative with density, chemical potential, and 
c..electron positron number densities
c..get the interpolation weight functions
       si0t   =  xpsi0(xt)
       si1t   =  xpsi1(xt)*dt_sav(jat)
       si0mt  =  xpsi0(mxt)
       si1mt  =  -xpsi1(mxt)*dt_sav(jat)
       si0d   =  xpsi0(xd)
       si1d   =  xpsi1(xd)*dd_sav(iat)
       si0md  =  xpsi0(mxd)
       si1md  =  -xpsi1(mxd)*dd_sav(iat)
c..derivatives of weight functions
       dsi0t  = xdpsi0(xt)*dti_sav(jat)
       dsi1t  = xdpsi1(xt)
       dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
       dsi1mt = xdpsi1(mxt)
       dsi0d  = xdpsi0(xd)*ddi_sav(iat)
       dsi1d  = xdpsi1(xd)
       dsi0md = -xdpsi0(mxd)*ddi_sav(iat)
       dsi1md = xdpsi1(mxd)
c..look in the pressure derivative only once
       fi(1)  = dpdf(iat,jat)
       fi(2)  = dpdf(iat+1,jat)
       fi(3)  = dpdf(iat,jat+1)
       fi(4)  = dpdf(iat+1,jat+1)
       fi(5)  = dpdft(iat,jat)
       fi(6)  = dpdft(iat+1,jat)
       fi(7)  = dpdft(iat,jat+1)
       fi(8)  = dpdft(iat+1,jat+1)
       fi(9)  = dpdfd(iat,jat)
       fi(10) = dpdfd(iat+1,jat)
       fi(11) = dpdfd(iat,jat+1)
       fi(12) = dpdfd(iat+1,jat+1)
       fi(13) = dpdfdt(iat,jat)
       fi(14) = dpdfdt(iat+1,jat)
       fi(15) = dpdfdt(iat,jat+1)
       fi(16) = dpdfdt(iat+1,jat+1)
c..pressure derivative with density
       dpepdd  = h3(iat,jat,
     1                 si0t,   si1t,   si0mt,   si1mt,
     2                 si0d,   si1d,   si0md,   si1md)
       dpepdd  = max(ye * dpepdd,0.0d0)
c..look in the electron chemical potential table only once
       fi(1)  = ef(iat,jat)
       fi(2)  = ef(iat+1,jat)
       fi(3)  = ef(iat,jat+1)
       fi(4)  = ef(iat+1,jat+1)
       fi(5)  = eft(iat,jat)
       fi(6)  = eft(iat+1,jat)
       fi(7)  = eft(iat,jat+1)
       fi(8)  = eft(iat+1,jat+1)
       fi(9)  = efd(iat,jat)
       fi(10) = efd(iat+1,jat)
       fi(11) = efd(iat,jat+1)
       fi(12) = efd(iat+1,jat+1)
       fi(13) = efdt(iat,jat)
       fi(14) = efdt(iat+1,jat)
       fi(15) = efdt(iat,jat+1)
       fi(16) = efdt(iat+1,jat+1)
c..electron chemical potential etaele
       etaele  = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2               si0d,   si1d,   si0md,   si1md)
c..derivative with respect to density
       x       = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2              dsi0d,  dsi1d,  dsi0md,  dsi1md)
       detadd  = ye * x
c..derivative with respect to temperature
       detadt  = h3(iat,jat,
     1              dsi0t,  dsi1t,  dsi0mt,  dsi1mt,
     2               si0d,   si1d,   si0md,   si1md)
c..derivative with respect to abar and zbar
      detada = -x * din * ytot1
      detadz =  x * den * ytot1
c..look in the number density table only once
       fi(1)  = xf(iat,jat)
       fi(2)  = xf(iat+1,jat)
       fi(3)  = xf(iat,jat+1)
       fi(4)  = xf(iat+1,jat+1)
       fi(5)  = xft(iat,jat)
       fi(6)  = xft(iat+1,jat)
       fi(7)  = xft(iat,jat+1)
       fi(8)  = xft(iat+1,jat+1)
       fi(9)  = xfd(iat,jat)
       fi(10) = xfd(iat+1,jat)
       fi(11) = xfd(iat,jat+1)
       fi(12) = xfd(iat+1,jat+1)
       fi(13) = xfdt(iat,jat)
       fi(14) = xfdt(iat+1,jat)
       fi(15) = xfdt(iat,jat+1)
       fi(16) = xfdt(iat+1,jat+1)
c..electron + positron number densities
      xnefer   = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2               si0d,   si1d,   si0md,   si1md)
c..derivative with respect to density
      x        = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2              dsi0d,  dsi1d,  dsi0md,  dsi1md)
      x = max(x,0.0d0)
      dxnedd   = ye * x
c..derivative with respect to temperature
      dxnedt   = h3(iat,jat,
     1              dsi0t,  dsi1t,  dsi0mt,  dsi1mt,
     2               si0d,   si1d,   si0md,   si1md)
c..derivative with respect to abar and zbar
      dxneda = -x * din * ytot1
      dxnedz =  x  * den * ytot1
c..the desired electron-positron thermodynamic quantities

c..dpepdd at high temperatures and low densities is below the
c..floating point limit of the subtraction of two large terms.
c..since dpresdd doesn't enter the maxwell relations at all, use the
c..bicubic interpolation done above instead of this one
       x       = din * din
       pele    = x * df_d
       dpepdt  = x * df_dt
c       dpepdd  = ye * (x * df_dd + 2.0d0 * din * df_d)
       s       = dpepdd/ye - 2.0d0 * din * df_d
       dpepda  = -ytot1 * (2.0d0 * pele + s * din)
       dpepdz  = den*ytot1*(2.0d0 * din * df_d  +  s)
       x       = ye * ye
       sele    = -df_t * ye
       dsepdt  = -df_tt * ye
       dsepdd  = -df_dt * x
       dsepda  = ytot1 * (ye * df_dt * din - sele)
       dsepdz  = -ytot1 * (ye * df_dt * den  + df_t)
       eele    = ye*free + temp * sele
       deepdt  = temp * dsepdt
       deepdd  = x * df_d + temp * dsepdd
       deepda  = -ye * ytot1 * (free +  df_d * din) + temp * dsepda
       deepdz  = ytot1* (free + ye * df_d * den) + temp * dsepdz
c..coulomb section:
c..initialize
       pcoul    = 0.0d0
       dpcouldd = 0.0d0
       dpcouldt = 0.0d0
       dpcoulda = 0.0d0
       dpcouldz = 0.0d0
       ecoul    = 0.0d0
       decouldd = 0.0d0
       decouldt = 0.0d0
       decoulda = 0.0d0
       decouldz = 0.0d0
       scoul    = 0.0d0
       dscouldd = 0.0d0
       dscouldt = 0.0d0
       dscoulda = 0.0d0
       dscouldz = 0.0d0
c..uniform background corrections only 
c..see yakovlev & shalybkov 1989 
c..lami is the average ion seperation
c..plasg is the plasma coupling parameter
       z        = forth * pi
       s        = z * xni
       dsdd     = z * dxnidd
       dsda     = z * dxnida
       lami     = 1.0d0/s**third
       inv_lami = 1.0d0/lami
       z        = -third * lami/s 
       lamidd   = z * dsdd
       lamida   = z * dsda
       plasg    = zbar*zbar*esqu*ktinv*inv_lami
       z        = -plasg * inv_lami 
       plasgdd  = z * lamidd
       plasgda  = z * lamida
       plasgdt  = -plasg*ktinv * kerg
       plasgdz  = 2.0d0 * plasg/zbar
c..fermi temperature, melting temperature and density condition
       y        = den * 1.0d-6
       x        = 1.009d0 * (y * ye)**third
       tfermi   = 5.930d9 * (sqrt(1.0d0 + x*x) - 1.0d0)
       tmelt    = 1.278d5 * zbar*zbar*(y*ytot1)**third
       rhocond  = abar*zbar
c..for any of these uniform background corrections to apply
c..the following conditions must hold
c..temp << t_fermi   electron gas is degenerate
c..temp >  t_melt    temperature exceeds the ion melting temperature
c..rho >>  abar*zbar electron gas is almost perfect, ionization full
       if ( (temp .lt. tfermi) .and.
     1      (temp .gt. 0.1d0 * tmelt)  .and.
     2      (den  .gt. 2.0d0 * rhocond)     
     3       ) then
c..yakovlev & shalybkov 1989 equations 82, 85, 86, 87
        if (plasg .ge. 1.0) then
         x        = plasg**(0.25d0) 
         y        = avo/abar * kerg 
         ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
         pcoul    = third * den * ecoul
         scoul    = -y * (3.0d0*b1*x - 5.0d0*c1/x
     1              + d1 * (log(plasg) - 1.0d0) - e1)
         y        = avo/abar*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
         decouldd = y * plasgdd 
         decouldt = y * plasgdt + ecoul/temp
         decoulda = y * plasgda - ecoul/abar
         decouldz = y * plasgdz
         y        = third * den
         dpcouldd = third * ecoul + y*decouldd
         dpcouldt = y * decouldt
         dpcoulda = y * decoulda
         dpcouldz = y * decouldz
         y        = -avo*kerg/(abar*plasg)*(0.75d0*b1*x+1.25d0*c1/x +d1)
         dscouldd = y * plasgdd
         dscouldt = y * plasgdt
         dscoulda = y * plasgda - scoul/abar
         dscouldz = y * plasgdz
c..yakovlev & shalybkov 1989 equations 102, 103, 104
        else if (plasg .lt. 1.0) then
         x        = plasg*sqrt(plasg)
         y        = plasg**b2
         z        = c2 * x - third * a2 * y
         pcoul    = -pion * z
         ecoul    = 3.0d0 * pcoul/den
         scoul    = -avo/abar*kerg*(c2*x -a2*(b2-1.0d0)/b2*y)
         s        = 1.5d0*c2*x/plasg - third*a2*b2*y/plasg
         dpcouldd = -dpiondd*z - pion*s*plasgdd
         dpcouldt = -dpiondt*z - pion*s*plasgdt
         dpcoulda = -dpionda*z - pion*s*plasgda
         dpcouldz = -dpiondz*z - pion*s*plasgdz
         s        = 3.0d0/den
         decouldd = s * dpcouldd - ecoul/den
         decouldt = s * dpcouldt
         decoulda = s * dpcoulda
         decouldz = s * dpcouldz
         s        = -avo*kerg/(abar*plasg)*(1.5d0*c2*x -a2*(b2-1.0d0)*y)
         dscouldd = s * plasgdd
         dscouldt = s * plasgdt
         dscoulda = s * plasgda - scoul/abar
         dscouldz = s * plasgdz
        end if
c..soften the edges of the coulomb correction turn-on
c..by doing linear interpolation to zero coulomb corrections
c..within reasonable limits of when the above expressions are valid
         z2 = 1.0d0
         if (den .gt. 2.0d0*rhocond .and. den .lt. 20.0d0*rhocond) then
          x  = den
          x1 = 2.0d0 * rhocond
          x2 = 20.0d0 * rhocond
          z2 = (x-x1)/(x2-x1)
         else if (temp .gt. 0.1d0 * tmelt  .and. temp .lt. tmelt) then
          x  = temp
          x1 = 0.1d0 * tmelt
          x2 = tmelt
          z2 = (x-x1)/(x2-x1)
         else if (temp .gt. 0.5d0 * tfermi  .and. temp .lt. tfermi) then
          x  = temp
          x1 = tfermi
          x2 = 0.5d0 * tfermi
          z2 = (x-x1)/(x2-x1)
         end if
         pcoul    = z2 * pcoul
         dpcouldd = z2 * dpcouldd
         dpcouldt = z2 * dpcouldt
         dpcoulda = z2 * dpcoulda
         dpcouldz = z2 * dpcouldz
         ecoul    = z2 * ecoul
         decouldd = z2 * decouldd
         decouldt = z2 * decouldt
         decoulda = z2 * decoulda
         decouldz = z2 * decouldz
         scoul  =   z2 * scoul
         dscouldd = z2 * dscouldd
         dscouldt = z2 * dscouldt
         dscoulda = z2 * dscoulda
         dscouldz = z2 * dscouldz
        end if
c..sum all the components
       pres    = prad + pion + pele + pcoul
       ener    = erad + eion + eele + ecoul
       entr    = srad + sion + sele + scoul
       dpresdd = dpraddd + dpiondd + dpepdd + dpcouldd 
       dpresdt = dpraddt + dpiondt + dpepdt + dpcouldt
       dpresda = dpradda + dpionda + dpepda + dpcoulda
       dpresdz = dpraddz + dpiondz + dpepdz + dpcouldz
       denerdd = deraddd + deiondd + deepdd + decouldd
       denerdt = deraddt + deiondt + deepdt + decouldt
       denerda = deradda + deionda + deepda + decoulda
       denerdz = deraddz + deiondz + deepdz + decouldz
       dentrdd = dsraddd + dsiondd + dsepdd + dscouldd
       dentrdt = dsraddt + dsiondt + dsepdt + dscouldt
       dentrda = dsradda + dsionda + dsepda + dscoulda
       dentrdz = dsraddz + dsiondz + dsepdz + dscouldz
c..the temperature and density exponents (c&g 9.81 9.82) 
c..the specific heat at constant volume (c&g 9.92)
c..the third adiabatic exponent (c&g 9.93)
c..the first adiabatic exponent (c&g 9.97) 
c..the second adiabatic exponent (c&g 9.105)
c..the specific heat at constant pressure (c&g 9.98) 
c..and relativistic formula for the sound speed (c&g 14.29)
       zz    = pres*deni
       zzi   = den/pres
       chit  = temp/pres * dpresdt
       chid  = dpresdd*zzi
       cv    = denerdt
       x     = zz * chit/(temp * cv)
       gam3  = x + 1.0d0
       gam1  = chit*x + chid
       nabad = x/gam1
       gam2  = 1.0d0/(1.0d0 - nabad)
       cp    = cv * gam1/chid
       z     = 1.0d0 + (ener + light2)*zzi
       sound = clight * sqrt(gam1/z)
c..maxwell relations; each is zero if the consistency is perfect
       x   = den * den
       dse = temp*dentrdt/denerdt - 1.0d0
       dpe = (denerdd*x + temp*dpresdt)/pres - 1.0d0
       dsp = -dentrdd*x/dpresdt - 1.0d0
c..store this row
        ptot_row(j)   = pres
        dpt_row(j)    = dpresdt
        dpd_row(j)    = dpresdd
        dpa_row(j)    = dpresda   
        dpz_row(j)    = dpresdz
        etot_row(j)   = ener
        det_row(j)    = denerdt
        ded_row(j)    = denerdd
        dea_row(j)    = denerda   
        dez_row(j)    = denerdz
        stot_row(j)   = entr 
        dst_row(j)    = dentrdt
        dsd_row(j)    = dentrdd
        dsa_row(j)    = dentrda        
        dsz_row(j)    = dentrdz
        prad_row(j)   = prad
        erad_row(j)   = erad
        srad_row(j)   = srad 
        pion_row(j)   = pion
        eion_row(j)   = eion
        sion_row(j)   = sion 
        xni_row(j)    = xni
        pele_row(j)   = pele
        ppos_row(j)   = 0.0d0
        dpept_row(j)  = dpepdt
        dpepd_row(j)  = dpepdd
        dpepa_row(j)  = dpepda  
        dpepz_row(j)  = dpepdz
        eele_row(j)   = eele
        epos_row(j)   = 0.0d0
        deept_row(j)  = deepdt
        deepd_row(j)  = deepdd
        deepa_row(j)  = deepda   
        deepz_row(j)  = deepdz
        sele_row(j)   = sele 
        spos_row(j)   = 0.0d0
        dsept_row(j)  = dsepdt 
        dsepd_row(j)  = dsepdd 
        dsepa_row(j)  = dsepda        
        dsepz_row(j)  = dsepdz
        xnem_row(j)   = xnem
        xne_row(j)    = xnefer
        dxnet_row(j)  = dxnedt
        dxned_row(j)  = dxnedd
        dxnea_row(j)  = dxneda
        dxnez_row(j)  = dxnedz
        xnp_row(j)    = 0.0d0
        etaele_row(j) = etaele
        detat_row(j)  = detadt
        detad_row(j)  = detadd
        detaa_row(j)  = detada
        detaz_row(j)  = detadz
        etapos_row(j) = 0.0d0
        pcou_row(j)   = pcoul
        ecou_row(j)   = ecoul
        scou_row(j)   = scoul 
        plasg_row(j)  = plasg
        dse_row(j)    = dse
        dpe_row(j)    = dpe
        dsp_row(j)    = dsp
        cv_row(j)     = cv
        cp_row(j)     = cp
        gam1_row(j)   = gam1
        gam2_row(j)   = gam2
        gam3_row(j)   = gam3
        cs_row(j)     = sound
c..end of vectorization loop
      enddo
      return
      end


c-------------------------------------------------------------
c    A program to test the tabular eos inversion procedure.
c-------------------------------------------------------------
      program eos_test
      implicit none
      real*8 rho, u, p, t, vx, vy, vz, eps, Ye, rhoh, cs
      real*8 ut, ux, uy, uz
      real*8 u_t, u_x, u_y, u_z
      real*8 d, e, sx, sy, sz
      real*8 rho0
      real*8 xi, delta
      real*8 gupxx, gupxy, gupxz, gupyy, gupyz, gupzz
      real*8 lapse, shiftx, shifty, shiftz, detgam
      real*8 shift_x, shift_y, shift_z
      real*8 ran1
      integer retval, idum
      real*8 gamma, rho_in, u_in, p_in, t_in
      real*8 vx_in, vy_in, vz_in
      integer eos_flag

      gamma = 2.d0
      eos_flag = 1

      ! Define metric variables.
      gupxx  = 1.d0
      gupxy  = 0.d0
      gupxz  = 0.d0
      gupyy  = 1.d0
      gupyz  = 0.d0
      gupzz  = 1.d0
      lapse  = 1.d0
      detgam = 1.d0
      shiftx = 0.d0
      shifty = 0.d0
      shiftz = 0.d0

      ! Choose density scale.
!      rho0 = 1.d14
      rho0 = 1.d5

      ! Choose some values for the primitives
c$$$      rho  = 1.d0
c$$$      eps  = 0.1
c$$$      Ye   = 0.5
c$$$      t    = 1.d11
c$$$      u    = rho * eps
c$$$      vx   = 0.1
c$$$      vy   = 0.02
c$$$      vz   = 0.2
      rho  = 10.d0
      eps  = 10.d0
      Ye   = 0.5
      t    = 2.6d10
      u    = rho * eps
      vx   = 0.0
      vy   = 0.0
      vz   = 0.0

      rho_in = rho
      u_in   = u
      vx_in  = vx
      vy_in  = vy
      vz_in  = vz
      
      call shen_initialize
      
      call print_range(rho,Ye,rho0)

      ! Call the eos wrapper to calculate remaining quantities.
      call eos_wrapper(p,rho,eps,Ye,t,rho0)

      t_in = t
      p_in = p

      write(*,*) 'Test of tabular EOS inversion'
      write(*,*) ' '
      write(*,*) '  Primitive variables: '
      write(*,*) '    rho = ', rho
      write(*,*) '    eps = ', eps
      write(*,*) '      u = ', u
      write(*,*) '     Ye = ', Ye
      write(*,*) '      t = ', t
      write(*,*) '      p = ', p
      write(*,*) '     vx = ', vx
      write(*,*) '     vy = ', vy
      write(*,*) '     vz = ', vz

      call compute_conservatives(d,e,sx,sy,sz,rho,p,u,vx,vy,vz,gamma, 
     &     eos_flag)

      write(*,*) ' '
      write(*,*) '  Conservative variables: '
      write(*,*) '    d  = ', d
      write(*,*) '    St = ', e
      write(*,*) '    Sx = ', sx
      write(*,*) '    Sy = ', sy
      write(*,*) '    Sz = ', sz

      ! Perturb primitives (and temperature) with some random
      ! deviations.
      delta = 1.d-2  ! the relative size of the perturbations
c      delta = 0.d0  ! the relative size of the perturbations
      idum = -1
      rho = rho*(1.d0+delta*(ran1(idum)-1.d0))
      eps = eps*(1.d0+delta*(ran1(idum)-1.d0))
      u   = rho * eps
      t   = t*(1.d0+delta*(ran1(idum)-1.d0))

      ! Write perturbed primitives to output.
      write(*,*) ' '
      write(*,*) '  Perturbed primitive variables: '
      write(*,*) '    rho = ', rho
      write(*,*) '    eps = ', eps
      write(*,*) '      u = ', u
      write(*,*) '      t = ', t
     
      ! Call the inversion scheme.
      call table_inversion(retval,rho,p,u,t,vx,vy,vz,cs,d,e,sx,sy,sz,
     &     Ye,rho0,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,lapse,
     &     detgam,shiftx,shifty,shiftz)

      ! Write results to output.
      write(*,*) ' '
      write(*,*) '  Results: '
      write(*,*) ' retval = ', retval
      write(*,*) '    rho = ', rho
      write(*,*) '    eps = ', u/rho
      write(*,*) '      u = ', u
      write(*,*) '      p = ', p
      write(*,*) '      t = ', t
      write(*,*) '     vx = ', vx
      write(*,*) '     vy = ', vy
      write(*,*) '     vz = ', vz
      write(*,*) '     cs = ', cs

      ! Error calculation
      write(*,*) ' '
      write(*,*) '  Errors:'
      write(*,*) '  drho = ', abs(rho-rho_in)/rho_in
      write(*,*) '  du   = ', abs(u-u_in)/u_in
      write(*,*) '  dp   = ', abs(p-p_in)/p_in
      write(*,*) '  dt   = ', abs(t-t_in)/t_in
      if (vx_in > 0.d0)   write(*,*) '  dvx  = ', abs(vx-vx_in)/vx_in
      if (vy_in > 0.d0)   write(*,*) '  dvy  = ', abs(vy-vy_in)/vy_in
      if (vz_in > 0.d0)   write(*,*) '  dvz  = ', abs(vz-vz_in)/vz_in
      
      end program

c     Branson typed this from the book.
      function ran1(idum)
      integer idum,ia,im,iq,ir,ntab,ndiv
      real*8 ran1, am, eps, rnmx
      parameter (ia=16807, im=2147483647,am=1.d0/im,iq=127773,ir=2836,
     &     ntab=32,ndiv=1+(im-1)/ntab,eps=1.2d-7,rnmx=1.d0-eps)
      integer j,k,iv(ntab),iy
      save iv,iy
      data iv /ntab*0/, iy /0/
      if (idum.le.0 .or. iy.eq.0) then
         idum = max(-idum,1)
         do j = ntab+8,1,-1
            k=idum/iq
            idum=ia*(idum-k*iq)-ir*k
            if (idum.lt.0) idum=idum+im
            if (j.le.ntab) iv(j)=idum
         end do
         iy = iv(1)
      endif
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      j=1+iy/ndiv
      iy=iv(j)
      iv(j)=idum
      ran1=min(am*iy,rnmx)
      return
      end

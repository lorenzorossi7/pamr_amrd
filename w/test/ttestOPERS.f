      program      ttest

      implicit     none

!-----------------------------------------------------------
!     Functions called
!-----------------------------------------------------------
      integer      iargc,      i4arg
      real*8       r8arg,      dmnrm2

      logical       ltrace 
      parameter   ( ltrace  = .true. ) 

      character*9   cdnm
      parameter   ( cdnm = 'ttest' )

      integer      level

      integer     maxN
      parameter ( maxN = 256 )
 
      real*8     fv    ( maxN, maxN, maxN )
      real*8     fv2   ( maxN, maxN, maxN )
      real*8     fc    ( maxN, maxN, maxN )
      real*8     fc2   ( maxN, maxN, maxN )

      real*8     chr_v ( maxN, maxN, maxN )
      real*8     chr_c ( maxN, maxN, maxN )

      real*8     xv ( maxN ),  xc ( maxN )
      real*8     yv ( maxN ),  yc ( maxN )
      real*8     zv ( maxN ),  zc ( maxN )

      real*8     xvf( maxN ),  xcf( maxN )
      real*8     yvf( maxN ),  ycf( maxN )
      real*8     zvf( maxN ),  zcf( maxN )

      real*8     xmin,   xmax, xminf, xmaxf
      real*8     ymin,   ymax, yminf, ymaxf
      real*8     zmin,   zmax, zminf, zmaxf 

      integer    type
      integer    order, do_ex
      real*8     ex

      integer    Nx,   Ny,   Nz 
      integer    Nxf,  Nyf,  Nzf 
      integer    Nx0,  Ny0,  Nz0 

      integer        getu, indlnb,do_out, argc, ret
      integer        get_int_param, get_real_param
      character*64   param_file,  optarg
      character*64   argv(64)
       
      integer        dim

!-----------------------------------------------------------
!     Parameter to describe the different data. 
!-----------------------------------------------------------

      real*8       x0,   y0,   z0
      real*8       sigx, sigy, sigz
      real*8       amp
      real*8       dx,   dy,   dz
      real*8       dxf,  dyf,  dzf

      integer      i,    j,    k 
      integer      i1,   j1,   k1
      integer      i2,   j2,   k2
      integer      ni,   nj,   nk

!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================

      do_ex = 0
      ex = 0.0d0

!-----------------------------------------------------------------------
!  Reading the name of the parameter file from the command line.
!-----------------------------------------------------------------------

      argc=iargc()
      do i=1,argc
        call getarg(i,argv(i))
      end do
      do i=1,argc
        param_file=argv(i)
      enddo 

!-----------------------------------------------------------------------
!  Reading the parameters from the parameter file.
!-----------------------------------------------------------------------

      ret = get_int_param (param_file, 'level',    level,     1)  
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter level'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        write(0,*) cdnm,': Ensure do you HAVE the parameter file ',
     &                       param_file 
        stop
      endif
      ret = get_int_param (param_file, 'type',   type,     1)  
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter type'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_int_param (param_file, 'dim',   dim,     1)
      if ( ret .ne. 1 ) then   
        write(0,*) cdnm,': Problems reading parameter dim'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_int_param (param_file, 'nx0',   nx0,     1)  
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter nx0'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_int_param (param_file, 'ny0',   ny0,    1)  
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter ny0'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_int_param (param_file, 'nz0',   nz0,    1)  
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter nz0'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_real_param (param_file, 'xmin',  xmin,   1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter xmin'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_real_param (param_file, 'xminf',  xminf,   1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter xminf'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_real_param (param_file, 'xmax',  xmax,   1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter xmax'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_real_param (param_file, 'xmaxf',  xmaxf,   1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter xmaxf'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_real_param (param_file, 'ymin',  ymin,   1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter ymin'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_real_param (param_file, 'yminf',  yminf,   1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter yminf'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_real_param (param_file, 'ymax',  ymax,   1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter ymax'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_real_param (param_file, 'ymaxf',  ymaxf,   1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter ymaxf'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_real_param (param_file, 'zmin',  zmin,   1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter zmin'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_real_param (param_file, 'zminf',  zminf,  1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter zminf'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_real_param (param_file, 'zmax',  zmax,   1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter zmax'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_real_param (param_file, 'zmaxf',  zmaxf,   1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter zmaxf'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_real_param(param_file, 'x0',      x0,       1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter x0'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_real_param(param_file, 'y0',      y0,       1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter y_1'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_real_param(param_file, 'z0',      z0,       1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter z_1'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_real_param(param_file, 'sigx',    sigx,       1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter sigx'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_real_param(param_file, 'sigy',      sigy,       1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter sigy'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_real_param(param_file, 'sigz',      sigz,       1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter sigz'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_real_param(param_file, 'amp',     amp,      1)  
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter amp'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_int_param (param_file, 'order',    order,     1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter order'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_int_param (param_file, 'i1',    i1,     1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter i1'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_int_param (param_file, 'j1',    j1,     1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter j1'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_int_param (param_file, 'k1',    k1,     1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter k1'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_int_param (param_file, 'i2',    i2,     1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter i2'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_int_param (param_file, 'j2',    j2,     1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter j2'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_int_param (param_file, 'k2',    k2,     1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter k2'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_int_param (param_file, 'ni',    ni,     1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter ni'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_int_param (param_file, 'nj',    nj,     1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter nj'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_int_param (param_file, 'nk',    nk,     1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter nk'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_int_param (param_file, 'nxf',   nxf,    1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter nxf'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_int_param (param_file, 'nyf',    nyf,     1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter nyf'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif
      ret = get_int_param (param_file, 'nzf',    nzf,     1)
      if ( ret .ne. 1 ) then
        write(0,*) cdnm,': Problems reading parameter nzf'
        write(0,*) cdnm,': Ensure it is on the parameter file'
        stop
      endif


!-----------------------------------------------------------------------
!  Trace the parameters read from the parameter file.
!-----------------------------------------------------------------------
      Nx = (Nx0-1)*2**(level)+1
      Ny = (Ny0-1)*2**(level)+1
      Nz = (Nz0-1)*2**(level)+1
      dx = (xmax-xmin)/(Nx-1) 
      dy = (ymax-ymin)/(Ny-1) 
      dz = (zmax-zmin)/(Nz-1) 
      if ( ltrace ) then
          write(0,*) cdnm,': level      = ',level
          write(0,*) cdnm,': xmin       = ',xmin
          write(0,*) cdnm,': xmax       = ',xmax
          write(0,*) cdnm,': ymin       = ',ymin
          write(0,*) cdnm,': ymax       = ',ymax
          write(0,*) cdnm,': zmin       = ',zmin
          write(0,*) cdnm,': zmax       = ',zmax
          write(0,*) cdnm,': xminf      = ',xminf
          write(0,*) cdnm,': xmaxf      = ',xmaxf
          write(0,*) cdnm,': yminf      = ',yminf
          write(0,*) cdnm,': ymaxf      = ',ymaxf
          write(0,*) cdnm,': zminf      = ',zminf
          write(0,*) cdnm,': zmaxf      = ',zmaxf
          write(0,*) cdnm,': sigx       = ',sigx
          write(0,*) cdnm,': sigy       = ',sigy
          write(0,*) cdnm,': sigz       = ',sigz
          write(0,*) cdnm,': Nx         = ',Nx
          write(0,*) cdnm,': Ny         = ',Ny
          write(0,*) cdnm,': Nz         = ',Nz
          write(0,*) cdnm,': dx         = ',dx
          write(0,*) cdnm,': dy         = ',dy
          write(0,*) cdnm,': dz         = ',dz
          write(0,*) cdnm,': amp        = ',amp
          write(0,*) cdnm,': order      = ',order
      endif

      if( level .lt. 0 ) go to 900

c-----------------------------------------------------------------------
c Definition of the coordinates for the vertices.
c-----------------------------------------------------------------------
      xv(1)  = xmin
      do i = 2, Nx-1
         xv(i) = xv(i-1) + dx 
      enddo
      xv(Nx) = xmax
      yv(1)  = ymin
      do j = 2, Ny-1
         yv(j) = yv(j-1) + dy
      enddo
      yv(Ny) = ymax
      zv(1)  = zmin
      do k = 2, Nz-1
         zv(k) = zv(k-1) + dz
      enddo
      zv(Nz) = zmax

c-----------------------------------------------------------------------
c Definition of the coordinates for the vertices.
c for the finer subgrid.
c-----------------------------------------------------------------------
      dxf = 0.5d0*dx
      dyf = 0.5d0*dy
      dzf = 0.5d0*dz
c-----------------------------------------------------------------------

      xvf(1)  = xmin + (i1-1)*dx - (i2-1)*dxf
      do i = 2, Nxf-1
         xvf(i) = xvf(i-1) + dxf
      enddo
      xvf(Nxf) = xvf(1) + ni*dx

      yvf(1)  = ymin + (j1-1)*dy - (j2-1)*dyf
      do j = 2, Nyf-1
         yvf(j) = yvf(j-1) + dyf
      enddo
      yvf(Nyf) = yvf(1) + nj*dy

      zvf(1)  = zmin + (k1-1)*dz - (k2-1)*dzf
      do k = 2, Nz-1
         zvf(k) = zvf(k-1) + dzf
      enddo 
      zvf(Nzf) = zv(1) + nk*dz

c-----------------------------------------------------------------------
c  Definition of the coordenates for the cell centres.
c-----------------------------------------------------------------------

      do i = 1, Nx-1
         xc(i) = xv(i) + 0.5d0*dx 
      enddo
      do j = 1, Ny-1
         yc(j) = yv(j) + 0.5d0*dy 
      enddo
      do k = 1, Nz-1
         zc(k) = zv(k) + 0.5d0*dz 
      enddo
      do i = 1, Nxf-1
         xcf(i) = xvf(i) + 0.5d0*dxf
      enddo
      do j = 1, Nyf-1
         ycf(j) = yvf(j) + 0.5d0*dyf
      enddo
      do k = 1, Nzf-1
         zcf(k) = zvf(k) + 0.5d0*dzf
      enddo
      if ( ltrace ) then
         write(0,*) cdnm,': dx  = ',dx 
         write(0,*) cdnm,': dxf = ',dxf
         write(0,*) cdnm,': dy  = ',dy 
         write(0,*) cdnm,': dyf = ',dyf
         write(0,*) cdnm,': dz  = ',dz 
         write(0,*) cdnm,': dzf = ',dzf
         write(0,*) cdnm,': xminf = ',xminf
         write(0,*) cdnm,': xmaxf = ',xmaxf
         write(0,*) cdnm,': yminf = ',yminf
         write(0,*) cdnm,': ymaxf = ',ymaxf
      endif
c-----------------------------------------------------------------------
c Trace out the coordinates.
c-----------------------------------------------------------------------
      if ( ltrace .and. .false. ) then  
        do i = 1, Nx 
           write(0,*) cdnm,': x(',i,') =',xv(i)
        enddo
        do j = 1, Ny 
           write(0,*) cdnm,': y(',j,') =',yv(j)
        enddo
        do k = 1, Nz
           write(0,*) cdnm,': z(',k,') =',zv(k)
        enddo
      endif
      if ( dim .eq. 1 ) then
          call initialize_fv(fc, amp, 
     &                       x0,     y0,     z0,
     &                       sigx,   sigy,   sigz,  
     &                       xc,     yc,     zc,
     &                       nx-1,   1,      1,
     &                       type)
      elseif ( dim .eq. 2 ) then
          call initialize_fv(fc, amp, 
     &                       x0,     y0,     z0,
     &                       sigx,   sigy,   sigz,  
     &                       xc,     yc,     zc,
     &                       nx-1,   ny-1,   1,
     &                       type)
      elseif ( dim .eq. 3 ) then
          call initialize_fv(fc, amp, 
     &                       x0,     y0,     z0,
     &                       sigx,   sigy,   sigz,  
     &                       xc,     yc,     zc,
     &                       nx-1,   ny-1,   nz-1,
     &                       type)
      endif  
      if ( ltrace ) then
         write(0,*) cdnm,': xminf = ',xminf
         write(0,*) cdnm,': i1 = ',i1
         write(0,*) cdnm,': xv(i1) = ',xv(i1)
         write(0,*) cdnm,': xvf(1) = ',xvf(1)
         write(0,*) cdnm,': ----'
         write(0,*) cdnm,': xc(i1) = ',xc(i1)
         write(0,*) cdnm,': xcf(1) = ',xcf(1)
         write(0,*) cdnm,': dx = ', dx
         write(0,*) cdnm,': dxf= ', dxf
      endif
      if ( dim .eq. 1 ) then
         call dminterp3d_c(fc,    fc2,
     &                     Nx-1,  Ny-1,  Nz-1 ,
     &                     Nxf-1, Nyf-1, Nzf-1,
     &                     i1,    j1,    k1,  
     &                     i2,    j2,    k2,  
     &                     ni-1,  1,     1, 
     &                     1,     1,     1,
     &                     2,     2,     2,
     &                     5,
     &                     chr_c, ex,    do_ex, dim)
      elseif ( dim .eq. 2 ) then 
         call dminterp3d_c(fc,    fc2,
     &                     Nx-1,  Ny-1,  Nz-1 ,
     &                     Nxf-1, Nyf-1, Nzf-1,
     &                     i1,    j1,    k1,  
     &                     i2,    j2,    k2,  
     &                     ni-1,  nj-1,  1, 
     &                     1,     1,     1,
     &                     2,     2,     2,
     &                     5,
     &                     chr_c, ex,    do_ex, dim)
      elseif ( dim .eq. 3 ) then 
         call dminterp3d_c(fc,    fc2,
     &                     Nx-1,  Ny-1,  Nz-1 ,
     &                     Nxf-1, Nyf-1, Nzf-1,
     &                     i1,    j1,    k1,  
     &                     i2,    j2,    k2,  
     &                     ni-1,  nj-1,  nk-1, 
     &                     1,     1,     1,
     &                     2,     2,     2,
     &                     5,
     &                     chr_c, ex,    do_ex, dim)
      endif   


c-----------------------------------------------------------------------
c Interpolating from the vertex to the cell centres.
c-----------------------------------------------------------------------
 
        if ( dim .eq. 1 ) then  
           call output_f(fc, 'fcccc', level,
     &                   xc,   yc,    zc,
     &                   nx-1, 1,  1, dim)

           call output_f(fc2,  'fcccc', level+1,
     &                   xcf,   ycf,    zcf,
     &                   nxf-1, 1,  1, dim)
        elseif ( dim .eq. 2 ) then 
           call output_f(fc, 'fcccc', level,
     &                   xc,   yc,    zc,
     &                   nx-1, ny-1,  1, dim)

           call output_f(fc2,  'fcccc', level+1,
     &                   xcf,   ycf,    zcf,
     &                   nxf-1, nyf-1, 1, dim)
        elseif ( dim .eq. 3 ) then 
           call output_f(fc,   'fcccc', level,
     &                   xc,    yc,     zc,
     &                   nx-1,  ny-1,   nz-1, dim)

           call output_f(fc2,  'fcccc', level+1,
     &                   xcf,   ycf,    zcf,
     &                   nxf-1, nyf-1,  nzf-1, dim)
       endif 
!-----------------------------------------------------------
!     Normal exit.
!-----------------------------------------------------------
      stop

!-----------------------------------------------------------
!     Usage exit.
!-----------------------------------------------------------
 900  continue
         write(0,*) 'usage: ttest <parameter file> '
      stop

      end


c----------------------------------------------------------------------
c numerical routines for the Special Relativistic Hydro Code.
c----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     Do one hydro step
c-----------------------------------------------------------------------
      subroutine hydro_1step(d,sx,sy,sz,e,
     &     d_p,sx_p,sy_p,sz_p,e_p,
     &     rho,vx,vy,vz,u,P,T,cs,x,y,z,gamma, 
     &     limiter,dt,Ye,rho0,Nx,Ny,Nz,
     &     fc_mask,d_fcs,sx_fcs,sy_fcs,sz_fcs,e_fcs,iter,eos_flag,
     &     phys_bdy,bc_type,N_bound,lorentz_max)
      implicit none
      integer Nx,Ny,Nz
c     ! Conserved variables
      real*8 d(Nx,Ny,Nz),d_p(Nx,Ny,Nz)
      real*8 sx(Nx,Ny,Nz),sx_p(Nx,Ny,Nz)
      real*8 sy(Nx,Ny,Nz),sy_p(Nx,Ny,Nz)
      real*8 sz(Nx,Ny,Nz),sz_p(Nx,Ny,Nz)
      real*8 e(Nx,Ny,Nz),e_p(Nx,Ny,Nz)
c     ! Primitive variables
      real*8 vx(Nx,Ny,Nz),vy(Nx,Ny,Nz),vz(Nx,Ny,Nz)
      real*8 u(Nx,Ny,Nz),rho(Nx,Ny,Nz)
      real*8 x(Nx),y(Ny),z(Nz),dt,gamma,Ye, rho0
      real*8 P(Nx,Ny,Nz), T(Nx,Ny,Nz), cs(Nx,Ny,Nz)
      integer eos_flag
      integer phys_bdy(6)
      real*8 lorentz_max
      
c     ! Flux correction stuff.
      integer fc_mask(Nx,Ny,Nz)
      integer lfc_mask, iter
      real*8 d_fcs(Nx,Ny,Nz), e_fcs(Nx,Ny,Nz)
      real*8 sx_fcs(Nx,Ny,Nz), sy_fcs(Nx,Ny,Nz), sz_fcs(Nx,Ny,Nz)
      real*8 lsign              ! the local sign with which fluxes enter the 
                                ! flux correction.  
      include 'fc_mask.inc'

c     ! Local arrays

      real*8 rhol(Nx,Ny,Nz),rhor(Nx,Ny,Nz)
      real*8 vxl(Nx,Ny,Nz),vxr(Nx,Ny,Nz)
      real*8 vyl(Nx,Ny,Nz),vyr(Nx,Ny,Nz)
      real*8 vzl(Nx,Ny,Nz),vzr(Nx,Ny,Nz)
      real*8 Pl(Nx,Ny,Nz),Pr(Nx,Ny,Nz)
      real*8 fd(Nx,Ny,Nz),fsx(Nx,Ny,Nz)
      real*8 fsy(Nx,Ny,Nz),fsz(Nx,Ny,Nz)
      real*8 fe(Nx,Ny,Nz)
      
      integer i,j,k, limiter, dir_x, dir_y, dir_z
      real*8 dx,dy,dz
      integer il,jl,kl
      integer N_bound           ! the number of points to treat as physical boundary.
                                ! This had better be the same as AMR_bdy_width_c.
      integer bc_type           ! 0 for copy, 1 for slanty.

      dx=(x(2)-x(1))
      dy=(y(2)-y(1))
      dz=(z(2)-z(1))

      dir_x = 1
      dir_y = 2
      dir_z = 3

      call cons_bc(bc_type,d,e,sx,sy,sz,Nx,Ny,Nz,N_bound,phys_bdy)

c     ! First call primitive variable inversion
      if (eos_flag == 0) then
         call find_primitives(d,sx,sy,
     &        sz,e,rho,vx,vy,vz,u,gamma, 
     &        Nx,Ny,Nz)
      else if (eos_flag == 1) then
         call find_primitives_table(d,sx,sy,
     &        sz,e,rho,vx,vy,vz,u,P,T,cs,Ye,rho0, 
     &        Nx,Ny,Nz)
      end if

c     ! Use EOS to calculate pressure
      if (eos_flag == 0) then
         do k=1,Nz
            do j=1,Ny
               do i=1,Nx
                  P(i,j,k) = (gamma-1.d0)*u(i,j,k)
               end do
            end do 
         end do
      end if

      do k=1,Nz
         do j=1,Ny
            do i=1,Nx
               fd(i,j,k)  = 0.d0
               fsx(i,j,k) = 0.d0
               fsy(i,j,k) = 0.d0
               fsz(i,j,k) = 0.d0
               fe(i,j,k)  = 0.d0
            end do
         end do 
      end do

      call hydro_advection(dir_x,limiter,rho,vx,vy,vz,P,u,cs, 
     &     fd,fsx,fsy,fsz,fe,gamma,dx,dy,dz,Nx,Ny,Nz,eos_flag,
     &     lorentz_max)

c     ! Assume zero flux through outer boundaries
      do k=1,Nz
         do j=1,Ny
            do i=1,Nx
               if (i==1) then
                  d_p(i,j,k)  =  - 1.d0/dx * fd(i+1,j,k)
                  sx_p(i,j,k) =  - 1.d0/dx * fsx(i+1,j,k)
                  sy_p(i,j,k) =  - 1.d0/dx * fsy(i+1,j,k)
                  sz_p(i,j,k) =  - 1.d0/dx * fsz(i+1,j,k)
                  e_p(i,j,k)  =  - 1.d0/dx * fe(i+1,j,k)
               else if (i==Nx) then
                  d_p(i,j,k)  = 1.d0/dx * fd(i,j,k)
                  sx_p(i,j,k) = 1.d0/dx * fsx(i,j,k)
                  sy_p(i,j,k) = 1.d0/dx * fsy(i,j,k)
                  sz_p(i,j,k) = 1.d0/dx * fsz(i,j,k)
                  e_p(i,j,k)  = 1.d0/dx * fe(i,j,k)
               else
                  d_p(i,j,k)  = - 1.d0/dx * (fd(i+1,j,k)-fd(i,j,k)) 
                  sx_p(i,j,k) = - 1.d0/dx * (fsx(i+1,j,k)-fsx(i,j,k)) 
                  sy_p(i,j,k) = - 1.d0/dx * (fsy(i+1,j,k)-fsy(i,j,k)) 
                  sz_p(i,j,k) = - 1.d0/dx * (fsz(i+1,j,k)-fsz(i,j,k)) 
                  e_p(i,j,k)  = - 1.d0/dx * (fe(i+1,j,k)-fe(i,j,k)) 
               end if
c     !        Save flux corrections if necessary.
c     !        We first cast the real value of the mask to an integer.
               lfc_mask = fc_mask(i,j,k)
               lsign = -1.0
               if (iand(lfc_mask,SIGN_FLAG) .ne. 0) lsign = 1.0
c     !        Note that lsign = 1 corresponds to the
c     !        natural sign, i.e., the sign with which
c     !        the fluxes contribute for evolution.
c     !        This natural sign should be used for type 
c     !        B cells.  

               if ((iand(lfc_mask,PXFLAG).ne.0) .and. 
     &                 iter==2) then

c     !           saving corrections at + x
                  d_fcs(i,j,k)  = d_fcs(i,j,k)  - 
     &                 lsign * dt/dx * fd(i+1,j,k)
                  sx_fcs(i,j,k) = sx_fcs(i,j,k) - 
     &                 lsign * dt/dx * fsx(i+1,j,k)
                  sy_fcs(i,j,k) = sy_fcs(i,j,k) - 
     &                 lsign * dt/dx * fsy(i+1,j,k)
                  sz_fcs(i,j,k) = sz_fcs(i,j,k) - 
     &                 lsign * dt/dx * fsz(i+1,j,k)
                  e_fcs(i,j,k)  = e_fcs(i,j,k)  - 
     &                 lsign * dt/dx * fe(i+1,j,k)
               end if
               if ((iand(lfc_mask,MXFLAG).ne.0) .and. 
     &              iter==2) then
c     !           saving corrections at - x
                  d_fcs(i,j,k)  = d_fcs(i,j,k)  +
     &                 lsign * dt/dx * fd(i,j,k)
                  sx_fcs(i,j,k) = sx_fcs(i,j,k) + 
     &                 lsign * dt/dx * fsx(i,j,k)
                  sy_fcs(i,j,k) = sy_fcs(i,j,k) + 
     &                 lsign * dt/dx * fsy(i,j,k)
                  sz_fcs(i,j,k) = sz_fcs(i,j,k) + 
     &                 lsign * dt/dx * fsz(i,j,k)
                  e_fcs(i,j,k)  = e_fcs(i,j,k)  + 
     &                 lsign * dt/dx * fe(i,j,k)
               end if
               
            end do
         end do
      end do
      
      call hydro_advection(dir_y,limiter,rho,vx,vy,vz,P,u,cs, 
     &     fd,fsx,fsy,fsz,fe,gamma,dx,dy,dz,Nx,Ny,Nz,eos_flag,
     &     lorentz_max)
      
c     ! add contribution
      do k=1,Nz
         do j=1,Ny
            do i=1,Nx
               if (j==1) then
                  d_p(i,j,k)  = d_p(i,j,k)  - 1.d0/dy * fd(i,j+1,k)
                  sx_p(i,j,k) = sx_p(i,j,k) - 1.d0/dy * fsx(i,j+1,k)
                  sy_p(i,j,k) = sy_p(i,j,k) - 1.d0/dy * fsy(i,j+1,k)
                  sz_p(i,j,k) = sz_p(i,j,k) - 1.d0/dy * fsz(i,j+1,k)
                  e_p(i,j,k)  = e_p(i,j,k)  - 1.d0/dy * fe(i,j+1,k)
               else if (j==Ny) then
                  d_p(i,j,k)  = d_p(i,j,k)  + 1.d0/dy * fd(i,j,k)
                  sx_p(i,j,k) = sx_p(i,j,k) + 1.d0/dy * fsx(i,j,k)
                  sy_p(i,j,k) = sy_p(i,j,k) + 1.d0/dy * fsy(i,j,k)
                  sz_p(i,j,k) = sz_p(i,j,k) + 1.d0/dy * fsz(i,j,k)
                  e_p(i,j,k)  = e_p(i,j,k)  + 1.d0/dy * fe(i,j,k)
               else
                  d_p(i,j,k)  = d_p(i,j,k)  - 1.d0/dy * (fd(i,j+1,k)
     &                 -fd(i,j,k)) 
                  sx_p(i,j,k) = sx_p(i,j,k) - 1.d0/dy * (fsx(i,j+1,k)
     &                 -fsx(i,j,k)) 
                  sy_p(i,j,k) = sy_p(i,j,k) - 1.d0/dy * (fsy(i,j+1,k)
     &                 -fsy(i,j,k)) 
                  sz_p(i,j,k) = sz_p(i,j,k) - 1.d0/dy * (fsz(i,j+1,k)
     &                 -fsz(i,j,k)) 
                  e_p(i,j,k)  = e_p(i,j,k)  - 1.d0/dy * (fe(i,j+1,k)
     &                 -fe(i,j,k)) 
               end if

c     !        Save flux corrections if necessary.
               lfc_mask = fc_mask(i,j,k)
               
               lsign = -1.0
               if (iand(lfc_mask,SIGN_FLAG) .ne. 0) lsign = 1.0
               if ((iand(lfc_mask,PYFLAG).ne.0)
     &              .and. iter==2) then
c     !           saving corrections at + y
                  d_fcs(i,j,k)  = d_fcs(i,j,k)  - 
     &                 lsign * dt/dy * fd(i,j+1,k)
                  sx_fcs(i,j,k) = sx_fcs(i,j,k) - 
     &                 lsign * dt/dy * fsx(i,j+1,k)
                  sy_fcs(i,j,k) = sy_fcs(i,j,k) - 
     &                 lsign * dt/dy * fsy(i,j+1,k)
                  sz_fcs(i,j,k) = sz_fcs(i,j,k) - 
     &                 lsign * dt/dy * fsz(i,j+1,k)
                  e_fcs(i,j,k)  = e_fcs(i,j,k)  - 
     &                 lsign * dt/dy * fe(i,j+1,k)
               end if
               if ((iand(lfc_mask,MYFLAG).ne.0) 
     &              .and. iter==2) then
c     !           saving corrections at - y
                  d_fcs(i,j,k)  = d_fcs(i,j,k)  +
     &                 lsign * dt/dy * fd(i,j,k)
                  sx_fcs(i,j,k) = sx_fcs(i,j,k) + 
     &                 lsign * dt/dy * fsx(i,j,k)
                  sy_fcs(i,j,k) = sy_fcs(i,j,k) + 
     &                 lsign * dt/dy * fsy(i,j,k)
                  sz_fcs(i,j,k) = sz_fcs(i,j,k) + 
     &                 lsign * dt/dy * fsz(i,j,k)
                  e_fcs(i,j,k)  = e_fcs(i,j,k)  + 
     &                 lsign * dt/dy * fe(i,j,k)
               end if
            end do
         end do
      end do
         
      call hydro_advection(dir_z,limiter,rho,vx,vy,vz,P,u,cs,
     &     fd,fsx,fsy,fsz,fe,gamma,dx,dy,dz,Nx,Ny,Nz,eos_flag,
     &     lorentz_max)
      

c     ! add contribution
      do k=1,Nz
         do j=1,Ny
            do i=1,Nx
               if (k==1) then 
                  d_p(i,j,k)  = d_p(i,j,k)  - 1.d0/dz * fd(i,j,k+1)
                  sx_p(i,j,k) = sx_p(i,j,k) - 1.d0/dz * fsx(i,j,k+1)
                  sy_p(i,j,k) = sy_p(i,j,k) - 1.d0/dz * fsy(i,j,k+1)
                  sz_p(i,j,k) = sz_p(i,j,k) - 1.d0/dz * fsz(i,j,k+1)
                  e_p(i,j,k)  = e_p(i,j,k)  - 1.d0/dz * fe(i,j,k+1)
               else if (k==Nz) then
                  d_p(i,j,k)  = d_p(i,j,k)  + 1.d0/dz * fd(i,j,k)
                  sx_p(i,j,k) = sx_p(i,j,k) + 1.d0/dz * fsx(i,j,k)
                  sy_p(i,j,k) = sy_p(i,j,k) + 1.d0/dz * fsy(i,j,k)
                  sz_p(i,j,k) = sz_p(i,j,k) + 1.d0/dz * fsz(i,j,k)
                  e_p(i,j,k)  = e_p(i,j,k)  + 1.d0/dz * fe(i,j,k)
               else 
                  d_p(i,j,k)  = d_p(i,j,k)  - 1.d0/dz * (fd(i,j,k+1)
     &                 -fd(i,j,k)) 
                  sx_p(i,j,k) = sx_p(i,j,k) - 1.d0/dz * (fsx(i,j,k+1)
     &                 -fsx(i,j,k)) 
                  sy_p(i,j,k) = sy_p(i,j,k) - 1.d0/dz * (fsy(i,j,k+1)
     &                 -fsy(i,j,k)) 
                  sz_p(i,j,k) = sz_p(i,j,k) - 1.d0/dz * (fsz(i,j,k+1)
     &                 -fsz(i,j,k)) 
                  e_p(i,j,k)  = e_p(i,j,k)  - 1.d0/dz * (fe(i,j,k+1)
     &                 -fe(i,j,k)) 
               end if
               
c     !        Save flux corrections if necessary.
               lfc_mask = fc_mask(i,j,k) 
               lsign = -1.0
               if (iand(lfc_mask,SIGN_FLAG).ne.0) lsign = 1.0
               if ((iand(lfc_mask,PZFLAG).ne.0)
     &              .and. iter==2) then
c     !           saving corrections at + z
                  d_fcs(i,j,k)  = d_fcs(i,j,k)  - 
     &                 lsign * dt/dz * fd(i,j,k+1)
                  sx_fcs(i,j,k) = sx_fcs(i,j,k) - 
     &                 lsign * dt/dz * fsx(i,j,k+1)
                  sy_fcs(i,j,k) = sy_fcs(i,j,k) - 
     &                 lsign * dt/dz * fsy(i,j,k+1)
                  sz_fcs(i,j,k) = sz_fcs(i,j,k) - 
     &                 lsign * dt/dz * fsz(i,j,k+1)
                  e_fcs(i,j,k)  = e_fcs(i,j,k)  - 
     &                 lsign * dt/dz * fe(i,j,k+1)
                        
               end if
               if ((iand(lfc_mask,MZFLAG).ne.0)
     &              .and. iter==2) then
c     !           saving corrections at - z
                  d_fcs(i,j,k)  = d_fcs(i,j,k)  +
     &                 lsign * dt/dz * fd(i,j,k)
                  sx_fcs(i,j,k) = sx_fcs(i,j,k) + 
     &                 lsign * dt/dz * fsx(i,j,k)
                  sy_fcs(i,j,k) = sy_fcs(i,j,k) + 
     &                 lsign * dt/dz * fsy(i,j,k)
                  sz_fcs(i,j,k) = sz_fcs(i,j,k) + 
     &                 lsign * dt/dz * fsz(i,j,k)
                  e_fcs(i,j,k)  = e_fcs(i,j,k)  + 
     &                 lsign * dt/dz * fe(i,j,k)
               end if
            end do
         end do
      end do
      
      if (bc_type == 1) then
c     !  handle lower diagonal corner
         do k = 1, Nz
            do j = 1, 2*N_bound
               do i= 1, 2*N_bound
                  d_p(i,j,k)    = 0.0
                  e_p(i,j,k)    = 0.0
                  sx_p(i,j,k)   = 0.0
                  sy_p(i,j,k)   = 0.0
                  sz_p(i,j,k)   = 0.0
               end do
            end do
         end do

c     !  handle upper diagonal corner
         do k = 1, Nz
            do j = Ny-2*N_bound+1, Ny
               do i= Nx-2*N_bound+1, Nx
                  d_p(i,j,k)    = 0.0
                  e_p(i,j,k)    = 0.0
                  sx_p(i,j,k)   = 0.0
                  sy_p(i,j,k)   = 0.0
                  sz_p(i,j,k)   = 0.0
               end do
            end do
         end do
      end if

      return
      end

c---------------------------------------------------------------
c     Subroutine for doing advection in 1 dimension
c---------------------------------------------------------------
      subroutine hydro_advection(dir_flag,limiter,rho,vx,vy,vz,P,u,cs, 
     &     fd,fsx,fsy,fsz,fe,gamma,dx,dy,dz,Nx,Ny,Nz,eos_flag,
     &     lorentz_max)
      implicit none
      integer Nx,Ny,Nz
c     ! Conserved variables
      real*8 rho(Nx,Ny,Nz), P(Nx,Ny,Nz), u(Nx,Ny,Nz)
      real*8 vx(Nx,Ny,Nz),vy(Nx,Ny,Nz),vz(Nx,Ny,Nz)
      real*8 fd(Nx,Ny,Nz),fsx(Nx,Ny,Nz)
      real*8 fsy(Nx,Ny,Nz),fsz(Nx,Ny,Nz)
      real*8 fe(Nx,Ny,Nz), cs(Nx,Ny,Nz)
      real*8 gamma
      real*8 dx, dy, dz
      integer eos_flag
      
c     ! Local arrays
      real*8 rhol(Nx,Ny,Nz),rhor(Nx,Ny,Nz)
      real*8 vxl(Nx,Ny,Nz),vxr(Nx,Ny,Nz)
      real*8 vyl(Nx,Ny,Nz),vyr(Nx,Ny,Nz)
      real*8 vzl(Nx,Ny,Nz),vzr(Nx,Ny,Nz)
      real*8 Pl(Nx,Ny,Nz),Pr(Nx,Ny,Nz)
      real*8 ul(Nx,Ny,Nz),ur(Nx,Ny,Nz)
      real*8 csl(Nx,Ny,Nz),csr(Nx,Ny,Nz)
      
      real*8 rhol_l, rhor_l, Pl_l, Pr_l
      real*8 vxl_l, vxr_l, vyl_l, vyr_l, vzl_l, vzr_l
      real*8 fdl, fdr, fsxl, fsxr
      real*8 fsyl, fsyr, fszl, fszr
      real*8 fel, fer, vil_l, vir_l
      real*8 cmax_l, cmin_l
      real*8 g_tt, g_tx, g_ty, g_tz
      real*8 g_xx, g_xy, g_xz, g_yy, g_yz, g_zz
      real*8 gup_tt, gup_tx, gup_ty, gup_tz
      real*8 gup_xx, gup_xy, gup_xz, gup_yy, gup_yz, gup_zz
      real*8 u_tl, u_xl, u_yl, u_zl
      real*8 u_tr, u_xr, u_yr, u_zr
      real*8 utl, uxl, uyl, uzl
      real*8 utr, uxr, uyr, uzr
      real*8 hl_l,hr_l, ul_l, ur_l
      real*8 cpr,cpl,cmr,cml,csl_l,csr_l
      real*8 dl_l, el_l, sxl_l, syl_l, szl_l
      real*8 dr_l, er_l, sxr_l, syr_l, szr_l
      
      real*8 sg
      integer i,j,k
      integer dir_flag, limiter
      real*8 lorentz_max, factor

c     ! Set metric functions once and for all.
      g_tt = -1.d0
      g_xx = 1.d0
      g_yy = 1.d0
      g_zz = 1.d0
      
      g_tx = 0.d0
      g_ty = 0.d0
      g_tz = 0.d0
      g_xy = 0.d0
      g_xz = 0.d0
      g_yz = 0.d0
      
      gup_tt = -1.d0
      gup_xx = 1.d0
      gup_yy = 1.d0
      gup_zz = 1.d0
      
      gup_tx = 0.d0
      gup_ty = 0.d0
      gup_tz = 0.d0
      gup_xy = 0.d0
      gup_xz = 0.d0
      gup_yz = 0.d0
      sg = 1.d0

c     ! Calculate left and right states
      call find_ur_ul(dir_flag,limiter,rho,rhol,rhor,
     &     dx,dy,dz,Nx,Ny,Nz)
      call find_ur_ul(dir_flag,limiter,P,Pl,Pr,dx,dy,dz,Nx,Ny,Nz)
      call find_ur_ul(dir_flag,limiter,vx,vxl,vxr,dx,dy,dz,Nx,Ny,Nz)
      call find_ur_ul(dir_flag,limiter,vy,vyl,vyr,dx,dy,dz,Nx,Ny,Nz)
      call find_ur_ul(dir_flag,limiter,vz,vzl,vzr,dx,dy,dz,Nx,Ny,Nz)
      if (eos_flag == 1) then
         call find_ur_ul(dir_flag,limiter,u,ul,ur,dx,dy,dz,Nx,Ny,Nz)
         call find_ur_ul(dir_flag,limiter,cs,csl,csr,dx,dy,dz,Nx,Ny,Nz)
      end if

c     ! calculate and return fluxes
c     ! For now, simply apply the HLL formula.
      do k=1,Nz
         do j=1,Ny
            do i=1,Nx
               rhol_l = rhol(i,j,k)
               rhor_l = rhor(i,j,k)
               vxl_l  = vxl(i,j,k)
               vxr_l  = vxr(i,j,k)
               vyl_l  = vyl(i,j,k)
               vyr_l  = vyr(i,j,k)
               vzl_l  = vzl(i,j,k)
               vzr_l  = vzr(i,j,k)
               Pl_l   = Pl(i,j,k)
               Pr_l   = Pr(i,j,k)

               if (eos_flag == 0) then
                  hl_l = 1 + gamma * Pl_l / (gamma-1.d0) / rhol_l                 
                  hr_l = 1 + gamma * Pr_l / (gamma-1.d0) / rhor_l
               else if (eos_flag == 1) then ! using tabular eos.
                  ul_l = ul(i,j,k)
                  ur_l = ur(i,j,k)
                  hl_l = 1 + (Pl_l + ul_l) / rhol_l                 
                  hr_l = 1 + (Pr_l + ur_l) / rhor_l
               end if

c     !        Calculate u_\mu,L and u_\mu,R locally 
               utl = g_tt + 2.d0*g_tx*vxl_l + 
     &              2.d0*g_ty*vyl_l + 2.d0*g_tz*vzl_l +
     &              g_xx*vxl_l**2 + 2.d0*g_xy*vxl_l*vyl_l + 
     &              2.d0*g_xz*vxl_l*vzl_l +
     &              g_yy*vyl_l**2 + 2.d0*g_yz*vyl_l*vzl_l + 
     &              g_zz*vzl_l**2
               
               if (utl .ge. 0.d0) then
                  write(*,*) 'superluminality on the left! i,j,k = ', 
     &                 i, ' ', j, ' ', k
                  factor = (1.d0-1.d0/lorentz_max**2)/(utl + 1.d0)
                  factor = sqrt(factor)
                  vxl_l = factor*vxl_l
                  vyl_l = factor*vyl_l
                  vzl_l = factor*vzl_l
                  utl = -1.d0/lorentz_max**2
               end if
               
               utl = 1.d0/sqrt(-utl)
               uxl = utl*vxl_l
               uyl = utl*vyl_l
               uzl = utl*vzl_l
               
               u_tl = g_tt*utl + g_tx*uxl + g_ty*uyl + g_tz*uzl
               u_xl = g_tx*utl + g_xx*uxl + g_xy*uyl + g_xz*uzl
               u_yl = g_ty*utl + g_xy*uxl + g_yy*uyl + g_yz*uzl
               u_zl = g_tz*utl + g_xz*uxl + g_yz*uyl + g_zz*uzl
               
c     !        Compute conservatives on left side.
               dl_l  = sg * rhol_l * utl
               el_l  = sg * rhol_l * hl_l * utl * u_tl + sg * Pl_l
               sxl_l = sg * rhol_l * hl_l * utl * u_xl
               syl_l = sg * rhol_l * hl_l * utl * u_yl
               szl_l = sg * rhol_l * hl_l * utl * u_zl

               utr = g_tt + 2.d0*g_tx*vxr_l + 
     &              2.d0*g_ty*vyr_l + 2.d0*g_tz*vzr_l + 
     &              g_xx*vxr_l**2 + 2.d0*g_xy*vxr_l*vyr_l + 
     &              2.d0*g_xz*vxr_l*vzr_l + 
     &              g_yy*vyr_l**2 + 2.d0*g_yz*vyr_l*vzr_l + 
     &              g_zz*vzr_l**2
           
               if (utr .ge. 0.d0) then
                  write(*,*) 'superluminality on the right! i,j,k = ', 
     &                 i, ' ', j, ' ', k
                  factor = (1.d0-1.d0/lorentz_max**2)/(utr + 1.d0)
                  factor = sqrt(factor)
                  vxr_l = factor*vxr_l
                  vyr_l = factor*vyr_l
                  vzr_l = factor*vzr_l
                  utr = -1.d0/lorentz_max**2
               end if
           
               utr = 1.d0/sqrt(-utr)
               uxr = utr*vxr_l
               uyr = utr*vyr_l
               uzr = utr*vzr_l
               
               u_tr = g_tt*utr + g_tx*uxr + g_ty*uyr + g_tz*uzr
               u_xr = g_tx*utr + g_xx*uxr + g_xy*uyr + g_xz*uzr
               u_yr = g_ty*utr + g_xy*uxr + g_yy*uyr + g_yz*uzr
               u_zr = g_tz*utr + g_xz*uxr + g_yz*uyr + g_zz*uzr
               
c     !        Compute conservatives on right side.
               dr_l  = sg * rhor_l * utr
               er_l  = sg * rhor_l * hr_l * utr * u_tr + sg * Pr_l
               sxr_l = sg * rhor_l * hr_l * utr * u_xr
               syr_l = sg * rhor_l * hr_l * utr * u_yr
               szr_l = sg * rhor_l * hr_l * utr * u_zr
                 
               if (dir_flag == 1) then 
                  fdl = sg * rhol_l * utl * vxl_l
                  fdr = sg * rhor_l * utr * vxr_l
                  fel = fdl * hl_l * u_tl
                  fer = fdr * hr_l * u_tr
                  fsxl = fdl * hl_l * u_xl + sg * Pl_l
                  fsxr = fdr * hr_l * u_xr + sg * Pr_l
                  fsyl = fdl * hl_l * u_yl 
                  fsyr = fdr * hr_l * u_yr 
                  fszl = fdl * hl_l * u_zl 
                  fszr = fdr * hr_l * u_zr 
               else if (dir_flag == 2) then 
                  fdl = sg * rhol_l * utl * vyl_l
                  fdr = sg * rhor_l * utr * vyr_l
                  fel = fdl * hl_l * u_tl
                  fer = fdr * hr_l * u_tr
                  fsxl = fdl * hl_l * u_xl
                  fsxr = fdr * hr_l * u_xr
                  fsyl = fdl * hl_l * u_yl + sg * Pl_l 
                  fsyr = fdr * hr_l * u_yr + sg * Pr_l
                  fszl = fdl * hl_l * u_zl 
                  fszr = fdr * hr_l * u_zr 
               else 
                  fdl = sg * rhol_l * utl * vzl_l
                  fdr = sg * rhor_l * utr * vzr_l
                  fel = fdl * hl_l * u_tl
                  fer = fdr * hr_l * u_tr
                  fsxl = fdl * hl_l * u_xl
                  fsxr = fdr * hr_l * u_xr
                  fsyl = fdl * hl_l * u_yl
                  fsyr = fdr * hr_l * u_yr
                  fszl = fdl * hl_l * u_zl + sg * Pl_l  
                  fszr = fdr * hr_l * u_zr + sg * Pr_l 
               end if
               
c     !        Calculate local comoving sound speed.
               if (eos_flag == 0) then
                  csl_l = sqrt(gamma*Pl_l/(rhol_l*hl_l))
                  csr_l = sqrt(gamma*Pr_l/(rhor_l*hr_l))
               else if (eos_flag == 1) then
                  csl_l = csl(i,j,k)
                  csr_l = csr(i,j,k)
               end if
               
c     !        Calculate c_+L, c_+R, c_-L, and c_-R
               call find_cp_cm(cpl,cml,csl_l,gup_tt,gup_tx,gup_ty, 
     &              gup_tz,gup_xx,gup_yy,gup_zz,utl,uxl,uyl,uzl,
     &              dir_flag)
               call find_cp_cm(cpr,cmr,csr_l,gup_tt,gup_tx,gup_ty,
     &              gup_tz,gup_xx,gup_yy,gup_zz,utr,uxr,uyr,uzr,
     &              dir_flag)
               
c     !        Calculate c_max and c_min 
               cmax_l = max(0.d0,cpr,cpl)
               cmin_l = -min(0.d0,cmr,cml)

c     !        Apply HLL formula
               fd(i,j,k) = cmin_l * fdr + cmax_l * fdl - 
     &              cmax_l*cmin_l*(dr_l - dl_l)
               fd(i,j,k) = fd(i,j,k) / (cmax_l + cmin_l)

               fsx(i,j,k) = cmin_l * fsxr + cmax_l * fsxl - 
     &              cmax_l*cmin_l*(sxr_l - sxl_l)
               fsx(i,j,k) = fsx(i,j,k) / (cmax_l + cmin_l)
               
               fsy(i,j,k) = cmin_l * fsyr + cmax_l * fsyl - 
     &              cmax_l*cmin_l*(syr_l - syl_l)
               fsy(i,j,k) = fsy(i,j,k) / (cmax_l + cmin_l)
               
               fsz(i,j,k) = cmin_l * fszr + cmax_l * fszl - 
     &              cmax_l*cmin_l*(szr_l - szl_l)
               fsz(i,j,k) = fsz(i,j,k) / (cmax_l + cmin_l)
               
               fe(i,j,k) =  cmin_l * fer + cmax_l * fel - 
     &              cmax_l*cmin_l*(er_l - el_l)
               fe(i,j,k) = fe(i,j,k) / (cmax_l + cmin_l)
               
            end do
         end do 
      end do
      
      return
      end

c-------------------------------------------------------------
c      Reconstruction step.  
c-------------------------------------------------------------
      subroutine find_ur_ul(dir_flag,limiter,f,fl,fr,dx,dy,dz,
     &     Nx,Ny,Nz)
      implicit none
      integer Nx,Ny,Nz,i,j,k,dir_flag
c     ! Conserved variables
      real*8 f(Nx,Ny,Nz), fr(Nx,Ny,Nz)
      real*8 fl(Nx,Ny,Nz), f_slopes(Nx,Ny,Nz)
      real*8 dx, dy, dz, dxo2, dyo2, dzo2
      integer limiter           ! limiter = 0 for zero slope
                                ! limiter = 1 for MC limiter

c     ! The mc_slopes routine doesn't
c     ! divide by dx, so you shouldn't multiply by it here.
      dxo2 = 1.d0 / 2.d0
      dyo2 = 1.d0 / 2.d0
      dzo2 = 1.d0 / 2.d0
      
      if (limiter == 0) then
         do k = 1, Nz
            do j = 1, Ny
               do i = 1, Nx
                  f_slopes(i,j,k) = 0.d0
               end do
            end do 
         end do
      else if (limiter == 1) then
         call mc_slopes(dir_flag,f,f_slopes,Nx,Ny,Nz)
      end if
      
      if (dir_flag == 1) then
         do k = 1, Nz
            do j = 1, Ny
               do i = 1, Nx
                  if (i==1) then
                     fr(i,j,k) = f(i,j,k)   - dxo2 * f_slopes(i,j,k)
                     fl(i,j,k) = fr(i,j,k)  
                  else
                     fr(i,j,k) = f(i,j,k)   - dxo2 * f_slopes(i,j,k)
                     fl(i,j,k) = f(i-1,j,k) + dxo2 * f_slopes(i-1,j,k)
                  end if
               end do
            end do
         end do
      else if (dir_flag == 2) then
         do k = 1, Nz
            do j = 1, Ny
               do i = 1, Nx
                  if (j==1) then
                     fr(i,j,k) = f(i,j,k)   - dyo2 * f_slopes(i,j,k)
                     fl(i,j,k) = fr(i,j,k)
                  else
                     fr(i,j,k) = f(i,j,k)   - dyo2 * f_slopes(i,j,k)
                     fl(i,j,k) = f(i,j-1,k) + dyo2 * f_slopes(i,j-1,k)
                  end if
               end do
            end do
         end do
      else if (dir_flag == 3) then
         do k = 1, Nz
            do j = 1, Ny
               do i = 1, Nx
                  if (k==1) then
                     fr(i,j,k) = f(i,j,k)   - dzo2 * f_slopes(i,j,k)
                     fl(i,j,k) = fr(i,j,k)
                  else
                     fr(i,j,k) = f(i,j,k)   - dzo2 * f_slopes(i,j,k)
                     fl(i,j,k) = f(i,j,k-1) + dzo2 * f_slopes(i,j,k-1)
                  end if
               end do
            end do
         end do
      end if
      
      return
      end 

c-----------------------------------------------------------------------
c     Initial data 
c     Problem types (parameter riemann_prob)
c     0 -- linear shock tube 
c     1 -- 2d Riemann problem
c     2 -- 2d cylindrical blast wave
c     3 -- shock at 45 degrees.
c-----------------------------------------------------------------------
      subroutine init_riemann(gamma,dn,sxn,syn,szn,en,rho,
     &     vx,vy,vz,u,P,T,cs,x_c,y_c,z_c,Nx,Ny,Nz,
     &     x0,rho_l,rho_r,u_l,u_r,ut_l,ux_l,uy_l,uz_l,
     &     ut_r,ux_r,uy_r,uz_r,T_r,T_l,Ye,rho0,eos_flag,
     &     riemann_prob,lorentz_max)
      implicit none
      integer Nx,Ny,Nz, riemann_prob
c     ! Conserved variables
      real*8 rho(Nx,Ny,Nz)
      real*8 dn(Nx,Ny,Nz)
      real*8 Sxn(Nx,Ny,Nz)
      real*8 Syn(Nx,Ny,Nz)
      real*8 Szn(Nx,Ny,Nz)
      real*8 En(Nx,Ny,Nz)
c     ! Primitive variables
      real*8 x0, rho_l, rho_r, u_l, u_r, rho0
      real*8 T_r, T_l      ! guesses for the temperature if necessary.
      integer eos_flag
      real*8 Ye
      real*8 ut_l,ux_l,uy_l,uz_l
      real*8 ut_r,ux_r,uy_r,uz_r
      real*8 vx(Nx,Ny,Nz),vy(Nx,Ny,Nz),vz(Nx,Ny,Nz)
      real*8 u(Nx,Ny,Nz),cs(Nx,Ny,Nz)
      real*8 P(Nx,Ny,Nz), T(Nx,Ny,Nz)
      real*8 dt,gamma
      real*8 x_c(Nx),y_c(Ny),z_c(Nz)
      real*8 r, PI
      integer i,j,k
      real*8 rho_atm, p_deplete
      real*8 rholoc, uloc, utloc, uxloc, uyloc, uzloc,ploc
      real*8 tloc, csloc, epsloc
      real*8 vxloc, vyloc, vzloc
      real*8 dloc,eloc,sxloc,syloc,szloc
      real*8 x_comp
      real*8 lorentz_max

      do k=1,Nz
         do j=1,Ny
            do i=1,Nx
               
               if (riemann_prob .eq. 1) then

                  if (x_c(i) < 0.5 .and. y_c(j) < 0.5) then
c     !              You're in the southwest
                     ploc = 1.0
                     rholoc = 0.5
                     vxloc = 0.d0
                     vyloc = 0.d0
                  else if (x_c(i) < 0.5 .and.  y_c(j) > 0.5) then
c     !              You're in the northwest
                     ploc = 1.0
                     rholoc = 0.1
                     vxloc = 0.99
                     vyloc = 0.d0
                  else if (x_c(i) > 0.5 .and. y_c(j) < 0.5) then
c     !              You're in the southeast
                     ploc = 1.0
                     rholoc = 0.1
                     vxloc = 0.d0
                     vyloc = 0.99
                  else if (x_c(i) > 0.5 .and. y_c(j) > 0.5) then
c     !              You're in the northeast
                     ploc = 0.01
                     rholoc = 0.1
                     vxloc = 0.d0
                     vyloc = 0.d0
                  end if
                  
                  if (eos_flag == 1) then
                     write(*,*) 'Burrows eos not supported'
                     stop
                  end if
                  
                  vzloc = 0.d0
                  
                  call compute_conservatives(dloc,eloc,sxloc,syloc,
     &                 szloc,rholoc,ploc,uloc,vxloc,vyloc,vzloc,gamma,
     &                 eos_flag,lorentz_max)
                  
                  dn(i,j,k)  = dloc
                  sxn(i,j,k) = sxloc
                  syn(i,j,k) = syloc
                  szn(i,j,k) = szloc
                  en(i,j,k)  = eloc
                  rho(i,j,k) = rholoc
                  vx(i,j,k)  = vxloc
                  vy(i,j,k)  = vyloc
                  vz(i,j,k)  = vzloc
                  u(i,j,k)   = uloc
                  P(i,j,k)   = ploc
                  T(i,j,k)   = tloc
               else
                  if (riemann_prob .eq. 0) x_comp = x_c(i)
                  if (riemann_prob .eq. 2) 
     &                 x_comp = sqrt(x_c(i)**2 + y_c(j)**2)
                  if (riemann_prob .eq. 3) x_comp = x_c(i)+y_c(j)

                  if (x_comp < x0) then
                     rholoc = rho_l
                     uloc   = u_l
                     utloc  = ut_l
                     uxloc  = ux_l
                     uyloc  = uy_l
                     uzloc  = uz_l
                     tloc   = T_l
                  else
                     rholoc = rho_r
                     uloc   = u_r
                     utloc  = ut_r
                     uxloc  = ux_r
                     uyloc  = uy_r
                     uzloc  = uz_r
                     tloc   = T_r
                  end if

                  if (eos_flag == 0) then
                     ploc   = (gamma-1.d0)*uloc
                     csloc  = 0.d0
                  else if (eos_flag == 1) then
                     epsloc = uloc / rholoc
                     call eos_wrapper(ploc,rholoc,epsloc,Ye,tloc,rho0,
     &                    csloc)
                  else
                     write(*,*) 'problem with eos_flag in riemann_id!'
                     stop
                  end if
                  
                  vxloc = uxloc/utloc
                  vyloc = uyloc/utloc
                  vzloc = uzloc/utloc
                  
                  call compute_conservatives(dloc,eloc,sxloc,syloc,
     &                 szloc, rholoc, ploc, uloc, vxloc, vyloc, vzloc, 
     &                 gamma, eos_flag, lorentz_max)
                  
                  dn(i,j,k)  = dloc
                  sxn(i,j,k) = sxloc
                  syn(i,j,k) = syloc
                  szn(i,j,k) = szloc
                  en(i,j,k)  = eloc
                  rho(i,j,k) = rholoc
                  vx(i,j,k)  = vxloc
                  vy(i,j,k)  = vyloc
                  vz(i,j,k)  = vzloc
                  u(i,j,k)   = uloc
                  T(i,j,k)   = tloc
                  P(i,j,k)   = ploc
                  cs(i,j,k)  = csloc
               end if
            end do
         end do
      end do
      
      return
      end

c-----------------------------------------------------------------------
c     Find primitive variables.
c-----------------------------------------------------------------------
      subroutine find_primitives(dn,sxn,syn,
     &     szn,en,rho,vx,vy,vz,u,gamma, 
     &     Nx,Ny,Nz)
      implicit none
      integer Nx,Ny,Nz
c     !  Conserved variables
      real*8 dn(Nx,Ny,Nz)
      real*8 Sxn(Nx,Ny,Nz)
      real*8 Syn(Nx,Ny,Nz)
      real*8 Szn(Nx,Ny,Nz)
      real*8 En(Nx,Ny,Nz)
c     ! Primitive variables
      real*8 vx(Nx,Ny,Nz),vy(Nx,Ny,Nz),vz(Nx,Ny,Nz)
      real*8 u(Nx,Ny,Nz),rho(Nx,Ny,Nz)
      real*8 gamma
c     ! Local arrays
      real*8 P(Nx,Ny,Nz)

      real*8 g_tt,   g_tx,   g_ty,   g_tz,   g_xx
      real*8 g_xy,   g_xz,   g_yy,   g_yz,   g_zz
      real*8 gup_tt, gup_tx, gup_ty, gup_tz, gup_xx
      real*8 gup_xy, gup_xz, gup_yy, gup_yz, gup_zz
      real*8 detg
      real*8 rholoc,   ploc,   vxloc,  vyloc,  vzloc
      real*8 dloc, eloc, sxloc, syloc, szloc
      
      integer i,j,k, retval
      
c     ! Set metric variables once and for all
      g_tt = -1.d0
      g_xx = 1.d0
      g_yy = 1.d0
      g_zz = 1.d0
      
      g_tx = 0.d0
      g_ty = 0.d0
      g_tz = 0.d0
      g_xy = 0.d0
      g_xz = 0.d0
      g_yz = 0.d0
      
      gup_tt = -1.d0
      gup_xx = 1.d0
      gup_yy = 1.d0
      gup_zz = 1.d0
      
      gup_tx = 0.d0
      gup_ty = 0.d0
      gup_tz = 0.d0
      gup_xy = 0.d0
      gup_xz = 0.d0
      gup_yz = 0.d0
      
      detg = 1.d0
      
      do k=1,Nz
         do j=1,Ny
            do i=1,Nx
               
               dloc  = dn(i,j,k)
               eloc  = en(i,j,k)
               sxloc = sxn(i,j,k)
               syloc = syn(i,j,k)
               szloc = szn(i,j,k)
               
               rholoc = rho(i,j,k)
               ploc = u(i,j,k)*(gamma-1.d0)
               vxloc = vx(i,j,k)
               vyloc = vy(i,j,k)
               vzloc = vz(i,j,k)

c     !        This is the wrapper for the Noble et al. solver.
               call inversion_interface(retval,
     &              rholoc,ploc,vxloc,vyloc,vzloc, 
     &              dloc,eloc,sxloc,syloc,szloc, 
     &              g_tt,g_tx,g_ty,g_tz,g_xx,g_xy,g_xz,g_yy,g_yz,g_zz, 
     &              gup_tt,gup_tx,gup_ty,gup_tz,gup_xx,gup_xy,gup_xz,
     &              gup_yy,gup_yz,gup_zz,detg,gamma)
               
               dn(i,j,k) = dloc
               en(i,j,k) = eloc
               sxn(i,j,k) = sxloc
               syn(i,j,k) = syloc
               szn(i,j,k) = szloc
               
               rho(i,j,k) = rholoc
               u(i,j,k)   = ploc/(gamma-1.d0)
               vx(i,j,k)  = vxloc
               vy(i,j,k)  = vyloc
               vz(i,j,k)  = vzloc
               
            end do
         end do
      end do
      
      return
      end

c-----------------------------------------------------------------------
c     Find primitive variables.  For a tabulated EOS.
c-----------------------------------------------------------------------
      subroutine find_primitives_table(dn,sxn,syn,
     &     szn,en,rho,vx,vy,vz,u,P,T,cs,Ye,rho0,
     &     Nx,Ny,Nz)
      implicit none
      integer Nx,Ny,Nz
c     ! Conserved variables
      real*8 dn(Nx,Ny,Nz)
      real*8 Sxn(Nx,Ny,Nz)
      real*8 Syn(Nx,Ny,Nz)
      real*8 Szn(Nx,Ny,Nz)
      real*8 En(Nx,Ny,Nz)
c     ! Primitive variables
      real*8 vx(Nx,Ny,Nz),vy(Nx,Ny,Nz),vz(Nx,Ny,Nz)
      real*8 u(Nx,Ny,Nz),rho(Nx,Ny,Nz),T(Nx,Ny,Nz)
      real*8 P(Nx,Ny,Nz), cs(Nx,Ny,Nz)
      real*8 Ye, rho0

      real*8 gamup_xx, gamup_xy, gamup_xz
      real*8 gamup_yy, gamup_yz, gamup_zz
      real*8 gamxx,gamxy,gamxz,gamyy,gamyz,gamzz
      real*8 detgam, lapse
      real*8 shiftx, shifty, shiftz
      real*8 rholoc,   ploc,   vxloc,  vyloc,  vzloc
      real*8 dloc, eloc, sxloc, syloc, szloc, tloc, csloc
      real*8 uloc
      
      integer i,j,k, retval
      
c     ! Set metric variables once and for all
      gamup_xx = 1.d0
      gamup_yy = 1.d0
      gamup_zz = 1.d0
      gamup_xy = 0.d0
      gamup_xz = 0.d0
      gamup_yz = 0.d0
      shiftx   = 0.d0
      shifty   = 0.d0
      shiftz   = 0.d0
      gamxx = 1.d0
      gamyy = 1.d0
      gamzz = 1.d0
      gamxy = 0.d0
      gamxz = 0.d0
      gamyz = 0.d0
      lapse = 1.d0
      detgam = 1.d0
      
      do k=1,Nz
         do j=1,Ny
            do i=1,Nx
               
               dloc  = dn(i,j,k)
               eloc  = en(i,j,k)
               sxloc = sxn(i,j,k)
               syloc = syn(i,j,k)
               szloc = szn(i,j,k)
               
               rholoc = rho(i,j,k)
               ploc   = P(i,j,k)
               vxloc  = vx(i,j,k)
               vyloc  = vy(i,j,k)
               vzloc  = vz(i,j,k)
               tloc   = T(i,j,k)
               uloc   = u(i,j,k)

               call table_inversion_loc(retval,rholoc,ploc,uloc,tloc,
     &              vxloc,vyloc,vzloc,csloc,dloc,eloc,sxloc,syloc,
     &              szloc,Ye,rho0,gamxx,gamxy,gamxz,gamyy,gamyz,
     &              gamzz,gamup_xx,gamup_xy,gamup_xz,
     &              gamup_yy,gamup_yz,gamup_zz,lapse,detgam,shiftx,
     &              shifty,shiftz,i,j,k)
                  
c     !        Assign outgoing values.
               dn(i,j,k)  = dloc
               en(i,j,k)  = eloc
               sxn(i,j,k) = sxloc
               syn(i,j,k) = syloc
               szn(i,j,k) = szloc
               
               rho(i,j,k) = rholoc
               u(i,j,k)   = uloc
               vx(i,j,k)  = vxloc
               vy(i,j,k)  = vyloc
               vz(i,j,k)  = vzloc
               T(i,j,k)   = tloc
               cs(i,j,k)  = csloc
               
            end do
         end do
      end do
      
      return
      end

c-------------------------------------------------------------
c     Compute slopes using the Monotonized-Centered limiter
c-------------------------------------------------------------
      subroutine mc_slopes(dir_flag, f, f_slope, Nx, Ny, Nz)
      implicit none
      integer Nx, Ny, Nz        
      integer  i, j, k
      integer  dir_flag
      real*8 f(Nx,Ny,Nz), f_slope(Nx,Ny,Nz)
      real*8  temp
      real*8 dfp, dfm, dfc
      
      if (dir_flag == 1) then
         do k = 1, Nz
            do j = 1, Ny
               do i = 2, Nx - 1
              
                  dfp = f(i+1,j,k) - f(i,j,k)
                  dfm = f(i,j,k) - f(i-1,j,k)
                  dfc = 0.5d0 * (dfp + dfm)
                  
                  if (dfp*dfm < 0.d0) then 
                     f_slope(i,j,k) = 0.d0
                  else
                     if (abs(2.d0*dfp) .lt. abs(2.d0*dfm)) then
                        temp = 2.d0*dfp
                     else
                        temp = 2.d0*dfm
                     end if
                     
                     if (abs(temp) .lt. abs(dfc)) then
                        f_slope(i,j,k) = temp
                     else
                        f_slope(i,j,k) = dfc
                     end if
                  end if
               end do
            end do
         end do
c     !  Set slope values at outer boundaries.  These
c     !  slopes are obviously incorrect, but everything
c     !  is fine as long as you have two ghost zones.
         do k = 1, Nz
            do j = 1, Ny
               f_slope(1,j,k)  = 0.d0
               f_slope(Nx,j,k) = 0.d0
            end do
         end do

      else if (dir_flag == 2) then
         
         do k = 1, Nz
            do j = 2, Ny - 1
               do i = 1, Nx 
                 
                  dfp = f(i,j+1,k) - f(i,j,k)
                  dfm = f(i,j,k) - f(i,j-1,k)
                  dfc = 0.5d0 * (dfp + dfm)
                 
                  if (dfp*dfm < 0.d0) then 
                     f_slope(i,j,k) = 0.d0
                  else
                     if (abs(2.d0*dfp) .lt. abs(2.d0*dfm)) then
                        temp = 2.d0*dfp
                     else
                        temp = 2.d0*dfm
                     end if
                     
                     if (abs(temp) .lt. abs(dfc)) then
                        f_slope(i,j,k) = temp
                     else
                        f_slope(i,j,k) = dfc
                     end if
                  end if
               end do
            end do
         end do
c     !  Set slope values at outer boundaries.  These
c     !  slopes are obviously incorrect, but everything
c     !  is fine as long as you have two ghost zones.
         do k = 1, Nz
            do i = 1, Nx
               f_slope(i,1,k)   = 0.d0
               f_slope(i,Ny,k)  = 0.d0
            end do
         end do

      else if (dir_flag == 3) then 
         
         do k = 2, Nz - 1
            do j = 1, Ny 
               do i = 1, Nx 
                  
                  dfp = f(i,j,k+1) - f(i,j,k)
                  dfm = f(i,j,k) - f(i,j,k-1)
                  dfc = 0.5d0 * (dfp + dfm)
                  
                  if (dfp*dfm < 0.d0) then 
                     f_slope(i,j,k) = 0.d0
                  else
                     if (abs(2.d0*dfp) .lt. abs(2.d0*dfm)) then
                        temp = 2.d0*dfp
                     else
                        temp = 2.d0*dfm
                     end if
                     
                     if (abs(temp) .lt. abs(dfc)) then
                        f_slope(i,j,k) = temp
                     else
                        f_slope(i,j,k) = dfc
                     end if
                  end if
               end do
            end do
         end do
c     !  Set slope values at outer boundaries.  These
c     !  slopes are obviously incorrect, but everything
c     !  is fine as long as you have two ghost zones.
         do j = 1, Ny
            do i = 1, Nx
               f_slope(i,j,1)  = 0.d0
               f_slope(i,j,Nz) = 0.d0
            end do
         end do

      end if

      end subroutine mc_slopes

c-----------------------------------------------------------------
c     I'm not sure whether it's actually necessary to initialize
c     the advanced time level like this.
c-----------------------------------------------------------------
      subroutine init_rest(dn,dnp1,sxn,sxnp1,syn,synp1,szn,sznp1,
     &     en,enp1,Nx,Ny,Nz)
      implicit none
      integer Nx,Ny,Nz
c     ! Conserved variables
      real*8 dn(Nx,Ny,Nz),dnp1(Nx,Ny,Nz)
      real*8 Sxn(Nx,Ny,Nz),Sxnp1(Nx,Ny,Nz)
      real*8 Syn(Nx,Ny,Nz),Synp1(Nx,Ny,Nz)
      real*8 Szn(Nx,Ny,Nz),Sznp1(Nx,Ny,Nz)
      real*8 En(Nx,Ny,Nz),Enp1(Nx,Ny,Nz)
c     ! Primitive variables
      integer i,j,k

      do k=1,Nz
         do j=1,Ny
            do i=1,Nx
               dnp1(i,j,k)  = dn(i,j,k)
               sxnp1(i,j,k) = sxn(i,j,k)
               synp1(i,j,k) = syn(i,j,k)
               sznp1(i,j,k) = szn(i,j,k)
               enp1(i,j,k)  = en(i,j,k)
            end do
         end do
      end do
      
      return
      end

c-----------------------------------------------------------------
c     Subroutine for zeroing fcs variables.  The flag determines
c     whether type a or type b cells should be zeroed.     
c-----------------------------------------------------------------
      subroutine zero_fcs_vars(type_flag,fc_mask,
     &     d_fcs,sx_fcs,sy_fcs,sz_fcs,e_fcs,Nx,Ny,Nz)
      implicit none
      integer Nx,Ny,Nz,type_flag
      integer fc_mask(Nx,Ny,Nz)
      integer lfc_mask
      real*8 d_fcs(Nx,Ny,Nz), e_fcs(Nx,Ny,Nz)
      real*8 sx_fcs(Nx,Ny,Nz), sy_fcs(Nx,Ny,Nz), sz_fcs(Nx,Ny,Nz)
      include 'fc_mask.inc'
      integer i,j,k

      do k=1,Nz
         do j=1,Ny
            do i=1,Nx
               
               lfc_mask = fc_mask(i,j,k)
               if (iand(lfc_mask,SIGN_FLAG) == type_flag) then 
                  if ((iand(lfc_mask,PXFLAG) .ne. 0) .or. 
     &                 (iand(lfc_mask,MXFLAG) .ne. 0) .or. 
     &                 (iand(lfc_mask,PYFLAG) .ne. 0) .or. 
     &                 (iand(lfc_mask,MYFLAG) .ne. 0) .or. 
     &                 (iand(lfc_mask,PZFLAG) .ne. 0) .or. 
     &                 (iand(lfc_mask,MZFLAG) .ne. 0)) then 
                     d_fcs(i,j,k)  = 0.d0
                     sx_fcs(i,j,k) = 0.d0
                     sy_fcs(i,j,k) = 0.d0
                     sz_fcs(i,j,k) = 0.d0
                     e_fcs(i,j,k)  = 0.d0
                  end if
               end if
            end do
         end do
      end do
      
      end subroutine

c-----------------------------------------------------------------
c     Subroutine for applying the flux correction.  This
c     probably could have been done in C.     
c-----------------------------------------------------------------
      subroutine apply_fc(dn,sxn,syn,szn,en,
     &     d_fcs,sx_fcs,sy_fcs,sz_fcs,e_fcs,Nx,Ny,Nz)
      implicit none
      integer Nx,Ny,Nz
      real*8 d_fcs(Nx,Ny,Nz), e_fcs(Nx,Ny,Nz)
      real*8 sx_fcs(Nx,Ny,Nz), sy_fcs(Nx,Ny,Nz), sz_fcs(Nx,Ny,Nz)
      real*8 dn(Nx,Ny,Nz), en(Nx,Ny,Nz)
      real*8 sxn(Nx,Ny,Nz), syn(Nx,Ny,Nz), szn(Nx,Ny,Nz)
      integer i,j,k

      do k=1,Nz
         do j=1,Ny
            do i=1,Nx
               dn(i,j,k)  = dn(i,j,k)  + d_fcs(i,j,k)
               sxn(i,j,k) = sxn(i,j,k) + sx_fcs(i,j,k)
               syn(i,j,k) = syn(i,j,k) + sy_fcs(i,j,k)
               szn(i,j,k) = szn(i,j,k) + sz_fcs(i,j,k)
               en(i,j,k)  = en(i,j,k)  + e_fcs(i,j,k)
               
            end do
         end do
      end do
      
      end subroutine

c--------------------------------------------------------------------
c     Compute the conserved (evolution) variables given the 
c     primitive variables.
c--------------------------------------------------------------------
      subroutine compute_conservatives(rho_s,s_t,s_x,s_y,s_z,rho,P,u,
     &     vx,vy,vz,GAMMA,eos_flag,lorentz_max)
      implicit none

      real*8 rho_s, s_t, s_x, s_y, s_z ! The conservative variables
      real*8 rho, P, vx, vy, vz, u     ! The primitive variables
      
      real*8 g_tt, g_tx, g_ty, g_tz
      real*8 g_xx, g_xy, g_xz, g_yy, g_yz, g_zz
      real*8 sg             ! sqrt(-g)
      real*8 GAMMA
      integer eos_flag
      
c     ! Variables for intermediate quantities.
      real*8  :: ut, ux, uy, uz, u_x, u_y, u_z, u_t, rhoh
      real*8  :: lorentz_max, factor

c     ! Set metric variables once and for all
      g_tt = -1.d0
      g_xx = 1.d0
      g_yy = 1.d0
      g_zz = 1.d0
      
      g_tx = 0.d0
      g_ty = 0.d0
      g_tz = 0.d0
      g_xy = 0.d0
      g_xz = 0.d0
      g_yz = 0.d0
      sg = 1.d0
      
      ut = g_tt + 2.d0*g_tx*vx + 2.d0*g_ty*vy + 2.d0*g_tz*vz + 
     &     g_xx*vx**2 + 2.d0*g_xy*vx*vy + 2.d0*g_xz*vx*vz + 
     &     g_yy*vy**2 + 2.d0*g_yz*vy*vz + g_zz*vz**2
      
      if (ut .ge. 0.d0) then
         write(*,*) 'superluminality!'
         factor = (1.d0-1.d0/lorentz_max**2)/(ut + 1.d0)
         factor = sqrt(factor)
         vx = factor*vx
         vy = factor*vy
         vz = factor*vz
         ut = -1.d0/lorentz_max**2
      end if

      ut = 1.d0/sqrt(-ut)
      ux = ut * vx
      uy = ut * vy
      uz = ut * vz
      
      u_t = g_tt * ut + g_tx * ux + g_ty * uy + g_tz * uz
      u_x = g_tx * ut + g_xx * ux + g_xy * uy + g_xz * uz
      u_y = g_ty * ut + g_xy * ux + g_yy * uy + g_yz * uz
      u_z = g_tz * ut + g_xz * ux + g_yz * uy + g_zz * uz
      
      if (eos_flag == 0) then
         rhoh = rho + GAMMA * P / (GAMMA - 1.d0)
      else if (eos_flag == 1) then
         rhoh = rho + P + u
      end if
      
      rho_s = sg * rho * ut
      s_t = sg * rhoh * ut * u_t + sg * P
      s_x = sg * rhoh * ut * u_x 
      s_y = sg * rhoh * ut * u_y 
      s_z = sg * rhoh * ut * u_z 

      end subroutine compute_conservatives

c--------------------------------------------------------------------
c     For a given fluid state, find the maximum left- and right-going
c     wave speeds.  (These are only equal in the comoving frame 
c     of the fluid.  We are calculating grid speeds here.)
c--------------------------------------------------------------------
      subroutine find_cp_cm(cp,cm,cs,gup_tt,gup_tx,gup_ty,gup_tz, 
     &     gup_xx,gup_yy,gup_zz,ut,ux,uy,uz,dir)
      implicit none

      real*8 cp,cm,cs
      real*8 gup_tt, gup_tx, gup_ty, gup_tz
      real*8 gup_xx, gup_yy, gup_zz
      real*8 ut, ux, uy, uz
      integer dir
      real*8 a, b, c, d, e, temp, desc

      d = ut
      a = d**2*(1.d0-cs**2) - cs**2*gup_tt
  
      if (dir == 1) then
         e = ux
         b = -2.d0*(d*e*(1.d0-cs**2) - cs**2*gup_tx)
         c = -gup_xx*cs**2 + e**2*(1.d0-cs**2)
      else if (dir == 2) then
         e = uy
         b = -2.d0*(d*e*(1.d0-cs**2) - cs**2*gup_ty)
         c = -gup_yy*cs**2 + e**2*(1.d0-cs**2)
      else if (dir == 3) then
         e = uz
         b = -2.d0*(d*e*(1.d0-cs**2) - cs**2*gup_tz)
         c = -gup_zz*cs**2 + e**2*(1.d0-cs**2)
      end if
      desc = b**2 - 4.d0 * a * c

      if (desc .lt. 0.d0) then
         write(*,*) 'descriminant less than zero! complexity!'
         stop
      end if

      cp = (-b + sqrt(desc))/2.d0/a
      cm = (-b - sqrt(desc))/2.d0/a

      end subroutine find_cp_cm

      subroutine find_cp_cm2(cp,cm,cs,gup_tt,gup_tx,gup_ty,gup_tz, 
     &     gup_xx,gup_yy,gup_zz,ut,ux,uy,uz,dir)
      implicit none

      real*8 cp,cm,cs
      real*8 gup_tt, gup_tx, gup_ty, gup_tz
      real*8 gup_xx, gup_yy, gup_zz
      real*8 ut, ux, uy, uz
      integer dir
      real*8 a, b, c, d, e, temp, desc
      real*8 vx, vy, vz

      if (dir==1) then
         vx = ux/ut
         cp = (vx + cs)/(1+cs*vx)
         cm = (vx - cs)/(1-cs*vx)
      else if (dir==2) then
         vy = uy/ut
         cp = (vy + cs)/(1+cs*vy)
         cm = (vy - cs)/(1-cs*vy)
      else
         vx = uz/ut
         cp = (vz + cs)/(1+cs*vz)
         cm = (vz - cs)/(1-cs*vz)
      end if

      end subroutine find_cp_cm2

c-------------------------------------------------------------
c      The refinement grid function.
c-------------------------------------------------------------
      subroutine find_tre_hydro(f_tre,rho,Nx,Ny,Nz)
      implicit none
      integer Nx,Ny,Nz,i,j,k
      real*8 f_tre(Nx,Ny,Nz), rho(Nx,Ny,Nz)
      
         do k = 1, Nz
            do j = 1, Ny
               do i = 1, Nx
                  if (i==1) then
                     f_tre(i,j,k) = 0.d0
                  else if (i==Nx) then
                     f_tre(i,j,k) = 0.d0
                  else
c     !              2nd derivative in the x-direction.
                     f_tre(i,j,k) = rho(i+1,j,k) - 2.d0*rho(i,j,k) + 
     &                    rho(i-1,j,k)
                  end if
               end do
            end do 
         end do
      
      return
      end 

        subroutine total_mass(m,rho,x,y,z,Nx,Ny,Nz,ghost_width,my_rank)
        implicit none
        integer Nx,Ny,Nz
        real*8 m
        real*8 rho(Nx,Ny,Nz)
        real*8 x(Nx+1),y(Ny+1),z(Nz+1)
        integer i, j, k
        integer is,ie,js,je,ks,ke
        real*8 dx, dy, dz
        integer ghost_width(6)
        integer my_rank
        
        dx = x(2)-x(1)
        dy = y(2)-y(1)
        dz = z(2)-z(1)

        is = ghost_width(1) + 1
        ie = Nx - ghost_width(2) 
        js = ghost_width(3) + 1
        je = Ny - ghost_width(4)
        ks = ghost_width(5) + 1
        ke = Nz - ghost_width(6)

        do i=is,ie
           do j=js,je
              do k=ks,ke
                 m = m + rho(i,j,k)
              end do
           end do
        end do
        
        m = m * dx * dy * dz

        return
        end

      subroutine cons_bc(bc_type,d,e,sx,sy,sz,Nx,Ny,Nz,N_bound,phys_bdy)
      implicit none
      integer Nx, Ny, Nz
      real*8 d(Nx,Ny,Nz)
      real*8 sx(Nx,Ny,Nz)
      real*8 sy(Nx,Ny,Nz)
      real*8 sz(Nx,Ny,Nz)
      real*8 e(Nx,Ny,Nz)
      integer phys_bdy(6)
      integer bc_type,N_bound
      integer i,j,k

      if (bc_type == 0) then
c     ! Apply copy boundary conditions.
c     ! Copy in x-direction
         if (phys_bdy(1) .ne. 0) then
            do k = 1, Nz
               do j = 1, Ny
                  do i= 1, N_bound
                     d(i,j,k)    = d(N_bound+1,j,k)
                     e(i,j,k)    = e(N_bound+1,j,k)
                     sx(i,j,k)   = sx(N_bound+1,j,k)
                     sy(i,j,k)   = sy(N_bound+1,j,k)
                     sz(i,j,k)   = sz(N_bound+1,j,k)
                  end do
               end do
            end do
         end if
         
         if (phys_bdy(2) .ne. 0) then
            do k = 1, Nz
               do j = 1, Ny
                  do i= Nx-N_bound+1, Nx
                     d(i,j,k)    = d(Nx-N_bound,j,k)
                     e(i,j,k)    = e(Nx-N_bound,j,k)
                     sx(i,j,k)   = sx(Nx-N_bound,j,k)
                     sy(i,j,k)   = sy(Nx-N_bound,j,k)
                     sz(i,j,k)   = sz(Nx-N_bound,j,k)
                  end do
               end do
            end do
         end if
      
c     ! Copy in y-direction
         if (phys_bdy(3) .ne. 0) then
            do k = 1, Nz
               do j = 1, N_bound
                  do i= 1, Nx
                     d(i,j,k)    = d(i,N_bound+1,k)
                     e(i,j,k)    = e(i,N_bound+1,k)
                     sx(i,j,k)   = sx(i,N_bound+1,k)
                     sy(i,j,k)   = sy(i,N_bound+1,k)
                     sz(i,j,k)   = sz(i,N_bound+1,k)
                  end do
               end do
            end do
         end if
         
         if (phys_bdy(4) .ne. 0) then
            do k = 1, Nz
               do j = Ny-N_bound+1, Ny
                  do i= 1, Nx
                     d(i,j,k)    = d(i,Ny-N_bound,k)
                     e(i,j,k)    = e(i,Ny-N_bound,k)
                     sx(i,j,k)   = sx(i,Ny-N_bound,k)
                     sy(i,j,k)   = sy(i,Ny-N_bound,k)
                     sz(i,j,k)   = sz(i,Ny-N_bound,k)
                  end do
               end do
            end do
         end if

c     ! Copy in z-direction
         if (phys_bdy(5) .ne. 0) then
            do k = 1, N_bound
               do j = 1, Ny
                  do i= 1, Nx
                     d(i,j,k)    = d(i,j,N_bound+1)
                     e(i,j,k)    = e(i,j,N_bound+1)
                     sx(i,j,k)   = sx(i,j,N_bound+1)
                     sy(i,j,k)   = sy(i,j,N_bound+1)
                     sz(i,j,k)   = sz(i,j,N_bound+1)
                  end do
               end do
            end do
         end if         

         if (phys_bdy(6) .ne. 0) then
            do k = Nz-N_bound+1, Nz
               do j = 1,Ny
                  do i= 1, Nx
                     d(i,j,k)    = d(i,j,Nz-N_bound)
                     e(i,j,k)    = e(i,j,Nz-N_bound)
                     sx(i,j,k)   = sx(i,j,Nz-N_bound)
                     sy(i,j,k)   = sy(i,j,Nz-N_bound)
                     sz(i,j,k)   = sz(i,j,Nz-N_bound)
                  end do
               end do
            end do
         end if
      end if

      if (bc_type == 1) then
c     ! Apply slanty boundary conditions.
c     ! handle lower x boundary
         if (phys_bdy(1) .ne. 0 .and. phys_bdy(3) .ne. 0) then
            do k = 1, Nz
               do j = N_bound+1, Ny
                  do i= 1,N_bound
                     d(i,j,k)    = d(i+N_bound,j-N_bound,k)
                     e(i,j,k)    = e(i+N_bound,j-N_bound,k)
                     sx(i,j,k)   = sx(i+N_bound,j-N_bound,k)
                     sy(i,j,k)   = sy(i+N_bound,j-N_bound,k)
                     sz(i,j,k)   = sz(i+N_bound,j-N_bound,k)
                  end do
               end do
            end do
         end if

         
c     ! handle the upper x boundary
         if (phys_bdy(2) .ne. 0 .and. phys_bdy(4) .ne. 0) then
            do k = 1, Nz
               do j = 1, Ny-N_bound
                  do i= Nx-N_bound+1, Nx
                     d(i,j,k)    = d(i-N_bound,j+N_bound,k)
                     e(i,j,k)    = e(i-N_bound,j+N_bound,k)
                     sx(i,j,k)   = sx(i-N_bound,j+N_bound,k)
                     sy(i,j,k)   = sy(i-N_bound,j+N_bound,k)
                     sz(i,j,k)   = sz(i-N_bound,j+N_bound,k)
                  end do
               end do
            end do
         end if
      
c     ! handle the lower y boundary
         if (phys_bdy(2) .ne. 0 .and. phys_bdy(3) .ne. 0) then
            do k = 1, Nz
               do j = 1, N_bound
                  do i= N_bound+1,Nx
                     d(i,j,k)    = d(i-N_bound,j+N_bound,k)
                     e(i,j,k)    = e(i-N_bound,j+N_bound,k)
                     sx(i,j,k)   = sx(i-N_bound,j+N_bound,k)
                     sy(i,j,k)   = sy(i-N_bound,j+N_bound,k)
                     sz(i,j,k)   = sz(i-N_bound,j+N_bound,k)
                  end do
               end do
            end do
         end if
         
c     ! handle the upper y boundary
         if (phys_bdy(1) .ne. 0 .and. phys_bdy(4) .ne. 0) then 
            do k = 1, Nz
               do j = Ny-N_bound+1, Ny
                  do i= 1, Nx-N_bound
                     d(i,j,k)    = d(i+N_bound,j-N_bound,k)
                     e(i,j,k)    = e(i+N_bound,j-N_bound,k)
                     sx(i,j,k)   = sx(i+N_bound,j-N_bound,k)
                     sy(i,j,k)   = sy(i+N_bound,j-N_bound,k)
                     sz(i,j,k)   = sz(i+N_bound,j-N_bound,k)
                  end do
               end do
            end do
         end if

c     ! Copy in z-direction
         if (phys_bdy(5) .ne. 0) then
            do k = 1, N_bound
               do j = 1, Ny
                  do i= 1, Nx
                     d(i,j,k)    = d(i,j,N_bound+1)
                     e(i,j,k)    = e(i,j,N_bound+1)
                     sx(i,j,k)   = sx(i,j,N_bound+1)
                     sy(i,j,k)   = sy(i,j,N_bound+1)
                     sz(i,j,k)   = sz(i,j,N_bound+1)
                  end do
               end do
            end do
         end if

         if (phys_bdy(6) .ne. 0) then
            do k = Nz-N_bound+1, Nz
               do j = 1,Ny
                  do i= 1, Nx
                     d(i,j,k)    = d(i,j,Nz-N_bound)
                     e(i,j,k)    = e(i,j,Nz-N_bound)
                     sx(i,j,k)   = sx(i,j,Nz-N_bound)
                     sy(i,j,k)   = sy(i,j,Nz-N_bound)
                     sz(i,j,k)   = sz(i,j,Nz-N_bound)
                  end do
               end do
            end do
         end if
      end if

      end subroutine cons_bc

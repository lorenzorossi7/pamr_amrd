c----------------------------------------------------------------------
c numerical routines for the Newtonian boson star code
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c    Assign initial data
c----------------------------------------------------------------------

c-----------------------------------------------------------------------
c differential operator L, computed where cmask=CMASK_ON (saved in LV)
c
c L(V) = 0
c 
c =>   V,xx + V,yy + V,zz - rho = 0
c
c with V=0 on the physical boundary
c Notice that you need to plug the vertex centered rho into here.
c-----------------------------------------------------------------------
        subroutine lop(LV,V,rho,cmask,x,y,z,Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz
        real*8 V(Nx+1,Ny+1,Nz+1),rho(Nx+1,Ny+1,Nz+1)
        real*8 cmask(Nx+1,Ny+1,Nz+1),LV(Nx+1,Ny+1,Nz+1)
        real*8 x(Nx+1),y(Ny+1),z(Nz+1)

        include 'cmask.inc'

        integer i,j,k
        real*8 dx,dy,dz,dx2,dy2,dz2

        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

        dx2=dx**2
        dy2=dy**2
        dz2=dz**2

        do k=2,Nz
           do j=2,Ny
              do i=2,Nx
                 if (cmask(i,j,k).eq.CMASK_ON) then
                    LV(i,j,k)=(V(i+1,j,k)-2*V(i,j,k)+V(i-1,j,k))/dx2
     &                   + (V(i,j+1,k)-2*V(i,j,k)+V(i,j-1,k))/dy2
     &                   + (V(i,j,k+1)-2*V(i,j,k)+V(i,j,k-1))/dz2
     &                   - rho(i,j,k)
                 end if
              end do
           end do
        end do
        
        return
        end
      

c-----------------------------------------------------------------------
c computes the residual L[V]-rhs L, computed where cmask=CMASK_ON 
c (saved in res) ... the L2 norm of the residual is stored in norm
c-----------------------------------------------------------------------
        subroutine residual(res,rhs,V,rho,cmask,x,y,z,norm,
     &                      Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz
        real*8 V(Nx+1,Ny+1,Nz+1),rho(Nx+1,Ny+1,Nz+1)
        real*8 cmask(Nx+1,Ny+1,Nz+1),res(Nx+1,Ny+1,Nz+1)
        real*8 rhs(Nx+1,Ny+1,Nz+1)
        real*8 x(Nx+1),y(Ny+1),z(Nz+1),norm

        integer i,j,k,sum
        include 'cmask.inc'

        call  lop(res,V,rho,cmask,x,y,z,Nx,Ny,Nz)

        norm=0
        sum=0

        do k=2,Nz
           do j=2,Ny
              do i=2,Nx
                 if (cmask(i,j,k).eq.CMASK_ON) then
                    res(i,j,k)=res(i,j,k)-rhs(i,j,k)
                    norm=norm+res(i,j,k)**2
                    sum=sum+1
                 end if
              end do
           end do
        end do

        norm=sqrt(norm/sum)
        return
        end

c-----------------------------------------------------------------------
c Applies a single red-black, gauss-seidel relaxation sweep to V
c cmask=CMASK_ON ... see Lop() for the equations
c
c the 'approximate' L2 norm of the residual prior to the sweep is
c returned
c
c phys_bdy=1 when corresponding boundary is physical
c-----------------------------------------------------------------------
        subroutine relax(V,rhs,rho,cmask,phys_bdy,x,y,z,norm,
     &                   Nx,Ny,Nz,ghost_width)
        implicit none
        integer Nx,Ny,Nz
        real*8 V(Nx+1,Ny+1,Nz+1),rho(Nx+1,Ny+1,Nz+1)
        real*8 cmask(Nx+1,Ny+1,Nz+1),rhs(Nx+1,Ny+1,Nz+1)
        real*8 x(Nx+1),y(Ny+1),z(Nz+1),norm
        integer phys_bdy(6)
        integer ghost_width(6)

        integer i,j,k,N,pass,red,sum
        real*8 PI
        real*8 dx,dy,dz,dx2,dy2,dz2
        real*8 res,Jac
        integer have_point

        include 'cmask.inc'

        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

        dx2=dx**2
        dy2=dy**2
        dz2=dz**2

        norm=0
        sum=0

        do pass=0,1
         do i=2,Nx
          do j=2,Ny
           do k=2,Nz
            if (mod(i+j+k+pass,2).eq.0.and.
     &          (cmask(i,j,k).eq.CMASK_ON)) then
             res=(V(i+1,j,k)-2*V(i,j,k)+V(i-1,j,k))/dx2
     &         + (V(i,j+1,k)-2*V(i,j,k)+V(i,j-1,k))/dy2
     &         + (V(i,j,k+1)-2*V(i,j,k)+V(i,j,k-1))/dz2
     &         - rho(i,j,k)
     &         - rhs(i,j,k)

             Jac=-2/dx2 -2/dy2 -2/dz2
      
             V(i,j,k)=V(i,j,k)-res/Jac

             norm=norm+res**2
             sum=sum+1
            end if
           end do
          end do
         end do
        end do

        norm=sqrt(norm/sum)

        return
        end

c-----------------------------------------------------------------------
c  Do one hydro step
c-----------------------------------------------------------------------
        subroutine hydro_1step(step_flag,d,sx,sy,sz,e,phi,
     &     d_p,sx_p,sy_p,sz_p,e_p,
     &     rho,vx,vy,vz,u,x,y,z,gamma, 
     &     limiter,rho_atm,dt,Nx,Ny,Nz,p_deplete,
     &     fc_mask,d_fcs,sx_fcs,sy_fcs,sz_fcs,e_fcs)
        implicit none
        integer Nx,Ny,Nz,step_flag
c     Conserved variables
        real*8 d(Nx,Ny,Nz), d_p(Nx,Ny,Nz)
        real*8 Sx(Nx,Ny,Nz), Sx_p(Nx,Ny,Nz)
        real*8 Sy(Nx,Ny,Nz), Sy_p(Nx,Ny,Nz)
        real*8 Sz(Nx,Ny,Nz), Sz_p(Nx,Ny,Nz)
        real*8 E(Nx,Ny,Nz), E_p(Nx,Ny,Nz)
        real*8 phi(Nx+1,Ny+1,Nz+1)
c     Primitive variables
        real*8 vx(Nx,Ny,Nz),vy(Nx,Ny,Nz),vz(Nx,Ny,Nz)
        real*8 u(Nx,Ny,Nz),rho(Nx,Ny,Nz)

        real*8 x(Nx),y(Ny),z(Nz),dt,gamma
        real*8 rho_atm, u_thresh, phil

c     Flux correction stuff.
        integer fc_mask(Nx,Ny,Nz)
        integer lfc_mask
        real*8 d_fcs(Nx,Ny,Nz), e_fcs(Nx,Ny,Nz)
        real*8 sx_fcs(Nx,Ny,Nz), sy_fcs(Nx,Ny,Nz), sz_fcs(Nx,Ny,Nz)
        real*8 lsign  ! the local sign with which fluxes enter the 
                      ! flux correction.  
        include 'fc_mask.inc'

c     Local arrays
        real*8 P(Nx,Ny,Nz)
        real*8 rhol(Nx,Ny,Nz),rhor(Nx,Ny,Nz)
        real*8 vxl(Nx,Ny,Nz),vxr(Nx,Ny,Nz)
        real*8 vyl(Nx,Ny,Nz),vyr(Nx,Ny,Nz)
        real*8 vzl(Nx,Ny,Nz),vzr(Nx,Ny,Nz)
        real*8 Pl(Nx,Ny,Nz),Pr(Nx,Ny,Nz)
        real*8 fd(Nx,Ny,Nz),fsx(Nx,Ny,Nz)
        real*8 fsy(Nx,Ny,Nz),fsz(Nx,Ny,Nz)
        real*8 fe(Nx,Ny,Nz)

        integer i,j,k, limiter, dir_x, dir_y, dir_z
        real*8 dx,dy,dz,dx2,dy2,dz2,p_deplete
        real*8 v2
        real*8 mass
        integer il,jl,kl
        integer have_point, have_point0, have_point1
        integer have_pointc0, have_pointc1
        real*8 temp_cont
        real*8 gamma1

        gamma1 = gamma

        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))
        dir_x = 1
        dir_y = 2
        dir_z = 3

c     First call primitive variable inversion
        call find_primitives(d,sx,sy,
     &       sz,e,phi,rho,vx,vy,vz,u,gamma, 
     &       Nx,Ny,Nz,rho_atm,p_deplete)

c     Use EOS to calculate pressure
        do k=1,Nz
           do j=1,Ny
              do i=1,Nx
                 P(i,j,k) = (gamma-1.d0)*u(i,j,k)
              end do
           end do 
        end do

        if (step_flag < 3) then 
c     call subroutine to advect in x-direction
           call hydro_advection(dir_x,limiter,rho,vx,vy,vz,P, 
     &          fd,fsx,fsy,fsz,fe,phi,gamma1,dx,dy,dz,Nx,Ny,Nz)

c     add contribution to conservatives at 
c     Assume zero flux through outer boundaries
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
                       d_p(i,j,k)  =  1.d0/dx * fd(i,j,k)
                       sx_p(i,j,k) =  1.d0/dx * fsx(i,j,k)
                       sy_p(i,j,k) =  1.d0/dx * fsy(i,j,k)
                       sz_p(i,j,k) =  1.d0/dx * fsz(i,j,k)
                       e_p(i,j,k)  =  1.d0/dx * fe(i,j,k)
                    else
                       d_p(i,j,k)  = -1.d0/dx*(fd(i+1,j,k)-fd(i,j,k)) 
                       sx_p(i,j,k) = -1.d0/dx*(fsx(i+1,j,k)-fsx(i,j,k)) 
                       sy_p(i,j,k) = -1.d0/dx*(fsy(i+1,j,k)-fsy(i,j,k)) 
                       sz_p(i,j,k) = -1.d0/dx*(fsz(i+1,j,k)-fsz(i,j,k)) 
                       e_p(i,j,k)  = -1.d0/dx*(fe(i+1,j,k)-fe(i,j,k)) 
                    end if

c     Save flux corrections if necessary.
c     We first cast the real value of the mask to an integer.
                    lfc_mask = fc_mask(i,j,k)
                    lsign = -1.0
                    if (iand(lfc_mask,SIGN_FLAG) .ne. 0) lsign = 1.0
c     Note that lsign = 1 corresponds to the
c     natural sign, i.e., the sign with which
c     the fluxes contribute for evolution.
c     This natural sign should be used for type 
c     B cells.  
                    if ((iand(lfc_mask,PXFLAG).ne.0) .and. 
     &                   step_flag==2) then
c     saving corrections at + x
                       d_fcs(i,j,k)  = d_fcs(i,j,k)  - 
     &                      lsign * dt/dx * fd(i+1,j,k)
                       sx_fcs(i,j,k) = sx_fcs(i,j,k) - 
     &                      lsign * dt/dx * fsx(i+1,j,k)
                       sy_fcs(i,j,k) = sy_fcs(i,j,k) - 
     &                      lsign * dt/dx * fsy(i+1,j,k)
                       sz_fcs(i,j,k) = sz_fcs(i,j,k) - 
     &                      lsign * dt/dx * fsz(i+1,j,k)
                       e_fcs(i,j,k)  = e_fcs(i,j,k)  - 
     &                      lsign * dt/dx * fe(i+1,j,k)

                    end if
                    if ((iand(lfc_mask,MXFLAG).ne.0) .and. 
     &                   step_flag==2) then
c     saving corrections at - x
                       d_fcs(i,j,k)  = d_fcs(i,j,k)  +
     &                      lsign * dt/dx * fd(i,j,k)
                       sx_fcs(i,j,k) = sx_fcs(i,j,k) + 
     &                      lsign * dt/dx * fsx(i,j,k)
                       sy_fcs(i,j,k) = sy_fcs(i,j,k) + 
     &                      lsign * dt/dx * fsy(i,j,k)
                       sz_fcs(i,j,k) = sz_fcs(i,j,k) + 
     &                      lsign * dt/dx * fsz(i,j,k)
                       e_fcs(i,j,k)  = e_fcs(i,j,k)  + 
     &                      lsign * dt/dx * fe(i,j,k)
                    end if
                 end do
              end do
           end do
           
           call hydro_advection(dir_y,limiter,rho,vx,vy,vz,P, 
     &          fd,fsx,fsy,fsz,fe,phi,gamma1,dx,dy,dz,Nx,Ny,Nz)
           
c     add contribution
           do k=1,Nz
              do j=1,Ny
                 do i=1,Nx
                    if (j==1) then
                       d_p(i,j,k)  = d_p(i,j,k)  - 1.d0/dy *fd(i,j+1,k)
                       sx_p(i,j,k) = sx_p(i,j,k) - 1.d0/dy *fsx(i,j+1,k)
                       sy_p(i,j,k) = sy_p(i,j,k) - 1.d0/dy *fsy(i,j+1,k)
                       sz_p(i,j,k) = sz_p(i,j,k) - 1.d0/dy *fsz(i,j+1,k)
                       e_p(i,j,k)  = e_p(i,j,k)  - 1.d0/dy *fe(i,j+1,k)
                    else if (j==Ny) then
                       d_p(i,j,k)  = d_p(i,j,k)  + 1.d0/dy * fd(i,j,k)
                       sx_p(i,j,k) = sx_p(i,j,k) + 1.d0/dy * fsx(i,j,k)
                       sy_p(i,j,k) = sy_p(i,j,k) + 1.d0/dy * fsy(i,j,k)
                       sz_p(i,j,k) = sz_p(i,j,k) + 1.d0/dy * fsz(i,j,k)
                       e_p(i,j,k)  = e_p(i,j,k)  + 1.d0/dy * fe(i,j,k)
                    else
                       d_p(i,j,k)  = d_p(i,j,k) -1.d0/dy*(fd(i,j+1,k) - 
     &                      fd(i,j,k)) 
                       sx_p(i,j,k) = sx_p(i,j,k) -1.d0/dy*(fsx(i,j+1,k)-
     &                      fsx(i,j,k)) 
                       sy_p(i,j,k) = sy_p(i,j,k) -1.d0/dy*(fsy(i,j+1,k)-
     &                      fsy(i,j,k)) 
                       sz_p(i,j,k) = sz_p(i,j,k) -1.d0/dy*(fsz(i,j+1,k)-
     &                      fsz(i,j,k)) 
                       e_p(i,j,k)  = e_p(i,j,k) -1.d0/dy*(fe(i,j+1,k)-
     &                      fe(i,j,k)) 
                    end if

c     Save flux corrections if necessary.
                    lfc_mask = fc_mask(i,j,k)
                    
                    lsign = -1.0
                    if (iand(lfc_mask,SIGN_FLAG) .ne. 0) lsign = 1.0
                    if ((iand(lfc_mask,PYFLAG).ne.0)
     &                   .and. step_flag==2) then
c     saving corrections at + y
                       d_fcs(i,j,k)  = d_fcs(i,j,k)  - 
     &                      lsign * dt/dy * fd(i,j+1,k)
                       sx_fcs(i,j,k) = sx_fcs(i,j,k) - 
     &                      lsign * dt/dy * fsx(i,j+1,k)
                       sy_fcs(i,j,k) = sy_fcs(i,j,k) - 
     &                      lsign * dt/dy * fsy(i,j+1,k)
                       sz_fcs(i,j,k) = sz_fcs(i,j,k) - 
     &                      lsign * dt/dx * fsz(i,j+1,k)
                       e_fcs(i,j,k)  = e_fcs(i,j,k)  - 
     &                      lsign * dt/dy * fe(i,j+1,k)
                    end if
                    if ((iand(lfc_mask,MYFLAG).ne.0) 
     &                   .and. step_flag==2) then
c     saving corrections at - y
                       d_fcs(i,j,k)  = d_fcs(i,j,k)  +
     &                      lsign * dt/dy * fd(i,j,k)
                       sx_fcs(i,j,k) = sx_fcs(i,j,k) + 
     &                      lsign * dt/dy * fsx(i,j,k)
                       sy_fcs(i,j,k) = sy_fcs(i,j,k) + 
     &                      lsign * dt/dy * fsy(i,j,k)
                       sz_fcs(i,j,k) = sz_fcs(i,j,k) + 
     &                      lsign * dt/dy * fsz(i,j,k)
                       e_fcs(i,j,k)  = e_fcs(i,j,k)  + 
     &                      lsign * dt/dy * fe(i,j,k)
                    end if
                 end do
              end do
           end do
           
           call hydro_advection(dir_z,limiter,rho,vx,vy,vz,P, 
     &          fd,fsx,fsy,fsz,fe,phi,gamma1,dx,dy,dz,Nx,Ny,Nz)
          
c     add contribution
           do k=1,Nz
              do j=1,Ny
                 do i=1,Nx
                    if (k==1) then 
                       d_p(i,j,k)  = d_p(i,j,k) - 1.d0/dz *fd(i,j,k+1)
                       sx_p(i,j,k) = sx_p(i,j,k) - 1.d0/dz *fsx(i,j,k+1)
                       sy_p(i,j,k) = sy_p(i,j,k) - 1.d0/dz *fsy(i,j,k+1)
                       sz_p(i,j,k) = sz_p(i,j,k) - 1.d0/dz *fsz(i,j,k+1)
                       e_p(i,j,k)  = e_p(i,j,k) - 1.d0/dz *fe(i,j,k+1)
                    else if (k==Nz) then
                       d_p(i,j,k)  = d_p(i,j,k) + 1.d0/dz * fd(i,j,k)
                       sx_p(i,j,k) = sx_p(i,j,k) + 1.d0/dz * fsx(i,j,k)
                       sy_p(i,j,k) = sy_p(i,j,k) + 1.d0/dz * fsy(i,j,k)
                       sz_p(i,j,k) = sz_p(i,j,k) + 1.d0/dz * fsz(i,j,k)
                       e_p(i,j,k)  = e_p(i,j,k) + 1.d0/dz * fe(i,j,k)
                    else 
                       d_p(i,j,k)  = d_p(i,j,k) -1.d0/dz*(fd(i,j,k+1)-
     &                      fd(i,j,k)) 
                       sx_p(i,j,k) = sx_p(i,j,k) -1.d0/dz*(fsx(i,j,k+1)-
     &                      fsx(i,j,k)) 
                       sy_p(i,j,k) = sy_p(i,j,k) -1.d0/dz*(fsy(i,j,k+1)-
     &                      fsy(i,j,k)) 
                       sz_p(i,j,k) = sz_p(i,j,k) -1.d0/dz*(fsz(i,j,k+1)-
     &                      fsz(i,j,k)) 
                       e_p(i,j,k)  = e_p(i,j,k) -1.d0/dz*(fe(i,j,k+1)-
     &                      fe(i,j,k)) 
                    end if

c     Save flux corrections if necessary.
                    lfc_mask = fc_mask(i,j,k) 
                    lsign = -1.0
                    if (iand(lfc_mask,SIGN_FLAG).ne.0) lsign = 1.0
                    if ((iand(lfc_mask,PZFLAG).ne.0)
     &                   .and. step_flag==2) then
c     saving corrections at + z
                       d_fcs(i,j,k)  = d_fcs(i,j,k)  - 
     &                      lsign * dt/dz * fd(i,j,k+1)
                       sx_fcs(i,j,k) = sx_fcs(i,j,k) - 
     &                      lsign * dt/dz * fsx(i,j,k+1)
                       sy_fcs(i,j,k) = sy_fcs(i,j,k) - 
     &                      lsign * dt/dz * fsy(i,j,k+1)
                       sz_fcs(i,j,k) = sz_fcs(i,j,k) - 
     &                      lsign * dt/dz * fsz(i,j,k+1)
                       e_fcs(i,j,k)  = e_fcs(i,j,k)  - 
     &                      lsign * dt/dz * fe(i,j,k+1)
                       
                    end if
                    if ((iand(lfc_mask,MZFLAG).ne.0)
     &                   .and. step_flag==2) then
c     saving corrections at - z
                       d_fcs(i,j,k)  = d_fcs(i,j,k)  +
     &                      lsign * dt/dz * fd(i,j,k)
                       sx_fcs(i,j,k) = sx_fcs(i,j,k) + 
     &                      lsign * dt/dz * fsx(i,j,k)
                       sy_fcs(i,j,k) = sy_fcs(i,j,k) + 
     &                      lsign * dt/dz * fsy(i,j,k)
                       sz_fcs(i,j,k) = sz_fcs(i,j,k) + 
     &                      lsign * dt/dz * fsz(i,j,k)
                       e_fcs(i,j,k)  = e_fcs(i,j,k)  + 
     &                      lsign * dt/dz * fe(i,j,k)
                    end if
                 end do
              end do
           end do

c     Add gravitational source terms to the momentum equation.
c     This is assuming that the phi array is vertex-centered.
c     But you have to average in the transverse directions.
           do k=1,Nz
              do j=1,Ny
                 do i=1,Nx
                    sx_p(i,j,k) = sx_p(i,j,k) - rho(i,j,k) * 
     &                   ((phi(i+1,j,k)-phi(i,j,k))/dx +
     &                   (phi(i+1,j+1,k)-phi(i,j+1,k))/dx +
     &                   (phi(i+1,j,k+1)-phi(i,j,k+1))/dx +
     &                   (phi(i+1,j+1,k+1)-phi(i,j+1,k+1))/dx)/4.d0
                    sy_p(i,j,k) = sy_p(i,j,k) - rho(i,j,k) * 
     &                   ((phi(i,j+1,k)-phi(i,j,k))/dy +
     &                   (phi(i+1,j+1,k)-phi(i+1,j,k))/dy +
     &                   (phi(i,j+1,k+1)-phi(i,j,k+1))/dy +
     &                   (phi(i+1,j+1,k+1)-phi(i+1,j,k+1))/dy)/4.d0
                    sz_p(i,j,k) = sz_p(i,j,k) - rho(i,j,k) * 
     &                   ((phi(i,j,k+1)-phi(i,j,k))/dz +
     &                   (phi(i+1,j,k+1)-phi(i+1,j,k))/dz +
     &                   (phi(i,j+1,k+1)-phi(i,j+1,k))/dz +
     &                   (phi(i+1,j+1,k+1)-phi(i+1,j+1,k))/dz)/4.d0
                 end do
              end do
           end do
           
        end if

        return
        end

      ! Subroutine for doing advection in 1 dimension
      subroutine hydro_advection(dir_flag,limiter,rho,vx,vy,vz,P, 
     &     fd,fsx,fsy,fsz,fe,phi,gamma,dx,dy,dz,Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz
        ! Conserved variables
        real*8 rho(Nx,Ny,Nz), P(Nx,Ny,Nz)
        real*8 vx(Nx,Ny,Nz),vy(Nx,Ny,Nz),vz(Nx,Ny,Nz)
        real*8 fd(Nx,Ny,Nz),fsx(Nx,Ny,Nz)
        real*8 fsy(Nx,Ny,Nz),fsz(Nx,Ny,Nz)
        real*8 fe(Nx,Ny,Nz)
        real*8 phi(Nx+1,Ny+1,Nz+1)
        real*8 gamma
        real*8 dx, dy, dz

        ! Local arrays
        real*8 rhol(Nx,Ny,Nz),rhor(Nx,Ny,Nz)
        real*8 vxl(Nx,Ny,Nz),vxr(Nx,Ny,Nz)
        real*8 vyl(Nx,Ny,Nz),vyr(Nx,Ny,Nz)
        real*8 vzl(Nx,Ny,Nz),vzr(Nx,Ny,Nz)
        real*8 Pl(Nx,Ny,Nz),Pr(Nx,Ny,Nz)
        real*8 cmax(Nx,Ny,Nz), cmin(Nx,Ny,Nz)

        real*8 rhol_l, rhor_l, Pl_l, Pr_l
        real*8 vxl_l, vxr_l, vyl_l, vyr_l, vzl_l, vzr_l
        real*8 fdl, fdr, fsxl, fsxr
        real*8 fsyl, fsyr, fszl, fszr
        real*8 fel, fer, vil_l, vir_l
        real*8 epsl_l, epsr_l, v2r, v2l
        real*8 phil_l, phir_l, el_l, er_l
        real*8 cmax_l, cmin_l

        integer i,j,k
        integer dir_flag, limiter
        real*8 gamma1 ! gamma governing adiabatic perturbations
       
        gamma1 = gamma
!        gamma1 = 5.d0/3.d0

        ! Calculate left and right states
        call find_ur_ul(dir_flag,limiter,rho,rhol,rhor,
     &       dx,dy,dz,Nx,Ny,Nz)
        call find_ur_ul(dir_flag,limiter,P,Pl,Pr,dx,dy,dz,Nx,Ny,Nz)
        call find_ur_ul(dir_flag,limiter,vx,vxl,vxr,dx,dy,dz,Nx,Ny,Nz)
        call find_ur_ul(dir_flag,limiter,vy,vyl,vyr,dx,dy,dz,Nx,Ny,Nz)
        call find_ur_ul(dir_flag,limiter,vz,vzl,vzr,dx,dy,dz,Nx,Ny,Nz)

        ! determine signal speeds 
        call find_cmax_cmin(dir_flag,cmax,cmin,vxl,vxr,vyl,vyr,vzl,vzr,
     &       Pl,Pr,rhol,rhor,gamma1,Nx,Ny,Nz)

        ! calculate and return fluxes
        ! For now, simply apply the HLL formula.
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
                 cmax_l = cmax(i,j,k)
                 cmin_l = cmin(i,j,k)

                 ! find phi on the left and right hand sides
                 ! Remember that you've passed in the vertex-centered phi
                 if (dir_flag == 1) then
                    phil_l = (phi(i,j,k) + 
     &                        phi(i,j+1,k) + 
     &                        phi(i,j,k+1) + 
     &                        phi(i,j+1,k+1))/4.d0 
                    phir_l = phil_l
                 else if (dir_flag == 2) then
                    phil_l = (phi(i,j,k) + 
     &                        phi(i+1,j,k) + 
     &                        phi(i,j,k+1) + 
     &                        phi(i+1,j,k+1))/4.d0 
                    phir_l = phil_l
                 else if (dir_flag == 3) then
                    phil_l = (phi(i,j,k) + 
     &                        phi(i+1,j,k) + 
     &                        phi(i,j+1,k) + 
     &                        phi(i+1,j+1,k))/4.d0 
                    phir_l = phil_l
                 end if   
                 
                 ! Note assumption of eos.
                 epsl_l = Pl_l/((gamma-1)*rhol_l)
                 epsr_l = Pr_l/((gamma-1)*rhor_l)
                 
                 if (dir_flag == 1) then
                    vil_l = vxl_l
                    vir_l = vxr_l
                 else if (dir_flag == 2) then
                    vil_l = vyl_l
                    vir_l = vyr_l
                 else if (dir_flag == 3) then
                    vil_l = vzl_l
                    vir_l = vzr_l
                 end if
                 
                 fdl = rhol_l*vil_l
                 fdr = rhor_l*vir_l
                 fsxl  = fdl*vxl_l
                 fsxr  = fdr*vxr_l
                 fsyl  = fdl*vyl_l
                 fsyr  = fdr*vyr_l
                 fszl  = fdl*vzl_l
                 fszr  = fdr*vzr_l


                 if (dir_flag == 1) then
                    fsxl = fsxl + Pl_l
                    fsxr = fsxr + Pr_l
                 else if (dir_flag == 2) then
                    fsyl = fsyl + Pl_l
                    fsyr = fsyr + Pr_l
                 else if (dir_flag == 3) then
                    fszl = fszl + Pl_l
                    fszr = fszr + Pr_l
                 end if         

                 v2l = vxl_l*vxl_l + vyl_l*vyl_l + vzl_l*vzl_l
                 v2r = vxr_l*vxr_l + vyr_l*vyr_l + vzr_l*vzr_l
                 fel   = fdl * (0.5d0*v2l + epsl_l + Pl_l/rhol_l 
     &                + phil_l)
                 fer   = fdr * (0.5d0*v2r + epsr_l + Pr_l/rhor_l 
     &                + phir_l)

                 el_l = rhol_l*epsl_l + 0.5d0*rhol_l*v2l +  
     &                0.5d0*rhol_l*phil_l
                 er_l = rhor_l*epsr_l + 0.5d0*rhor_l*v2r +  
     &                0.5d0*rhor_l*phir_l

                 ! Apply HLL formula
                 fd(i,j,k) = cmin_l * fdr + cmax_l * fdl - 
     &                cmax_l*cmin_l*(rhor_l - rhol_l)
                 fd(i,j,k) = fd(i,j,k) / (cmax_l + cmin_l)

                 fsx(i,j,k) = cmin_l * fsxr + cmax_l * fsxl - 
     &                cmax_l*cmin_l*(rhor_l*vxr_l - rhol_l*vxl_l)
                 fsx(i,j,k) = fsx(i,j,k) / (cmax_l + cmin_l)

                 fsy(i,j,k) = cmin_l * fsyr + cmax_l * fsyl - 
     &                cmax_l*cmin_l*(rhor_l*vyr_l - rhol_l*vyl_l)
                 fsy(i,j,k) = fsy(i,j,k) / (cmax_l + cmin_l)

                 fsz(i,j,k) = cmin_l * fszr + cmax_l * fszl - 
     &                cmax_l*cmin_l*(rhor_l*vzr_l - rhol_l*vzl_l)
                 fsz(i,j,k) = fsz(i,j,k) / (cmax_l + cmin_l)

                 fe(i,j,k) =  cmin_l * fer + cmax_l * fel - 
     &                cmax_l*cmin_l*(er_l - el_l)
                 fe(i,j,k) = fe(i,j,k) / (cmax_l + cmin_l)

              end do
           end do 
        end do
     
        return
        end

! Reconstruction step.  Zero slope for now.
        subroutine find_ur_ul(dir_flag,limiter,f,fl,fr,dx,dy,dz,
     &     Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz,i,j,k,dir_flag
        ! Conserved variables
        real*8 f(Nx,Ny,Nz), fr(Nx,Ny,Nz)
        real*8 fl(Nx,Ny,Nz), f_slopes(Nx,Ny,Nz)
        real*8 dx, dy, dz, dxo2, dyo2, dzo2
        integer limiter   ! limiter = 0 for zero slope
                          ! limiter = 1 for MC limiter

        ! You silly goose!  The mc_slopes routine doesn't
        ! divide by dx, so you shouldn't multiply by it here.
c$$$        dxo2 = dx / 2.d0
c$$$        dyo2 = dy / 2.d0
c$$$        dzo2 = dz / 2.d0
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
        
        end 

! Find maximum left and right propagating wave speeds.
! Use galilean transformations!
        subroutine find_cmax_cmin(dir_flag,cmax,cmin,vxl,vxr,vyl,vyr, 
     &     vzl,vzr,Pl,Pr,rhol,rhor,gamma,Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz,i,j,k,dir_flag
        ! Conserved variables
        real*8 cmax(Nx,Ny,Nz),cmin(Nx,Ny,Nz)
        real*8 vxl(Nx,Ny,Nz),vxr(Nx,Ny,Nz)
        real*8 vyl(Nx,Ny,Nz),vyr(Nx,Ny,Nz)
        real*8 vzl(Nx,Ny,Nz),vzr(Nx,Ny,Nz)
        real*8 rhol(Nx,Ny,Nz),rhor(Nx,Ny,Nz)
        real*8 Pl(Nx,Ny,Nz), Pr(Nx,Ny,Nz)
        real*8 gamma

        real*8 csl, csr

        do k=1,Nz
           do j=1,Ny
              do i=1,Nx
                 csl = sqrt(Pl(i,j,k)*gamma/rhol(i,j,k))
                 csr = sqrt(Pr(i,j,k)*gamma/rhor(i,j,k))

                 if (dir_flag==1) then
                    cmax(i,j,k) = max(0.d0,csl+vxl(i,j,k),
     &                   csr+vxr(i,j,k))
                    cmin(i,j,k) = -min(0.d0,-csl+vxl(i,j,k), 
     &                   -csr+vxr(i,j,k))

                 else if (dir_flag==2) then
                    cmax(i,j,k) = max(0.d0,csl+vyl(i,j,k), 
     &                   csr+vyr(i,j,k))
                    cmin(i,j,k) = -  
     &                   min(0.d0,-csl+vyl(i,j,k),-csr+vyr(i,j,k))
                 else
                    cmax(i,j,k) = max(0.d0,csl+vzl(i,j,k),
     &                   csr+vzr(i,j,k))
                    cmin(i,j,k) = -  
     &                   min(0.d0,-csl+vzl(i,j,k),-csr+vzr(i,j,k))
                 end if
              end do
           end do
        end do

        end

c-----------------------------------------------------------------------
c  Initial data for the polytrope, n=1 case
c-----------------------------------------------------------------------
        subroutine polytrope_id(gamma,dn,sxn,syn,szn,en,rho,rhov,
     &     vx,vy,vz,u,phiv,p_deplete,x,y,z,x_c,y_c,z_c,
     &     Nx,Ny,Nz,rho_atm)
        implicit none
        integer Nx,Ny,Nz
        ! Conserved variables
        real*8 rho(Nx,Ny,Nz)
        real*8 rhov(Nx+1,Ny+1,Nz+1)
        real*8 dn(Nx,Ny,Nz)
        real*8 Sxn(Nx,Ny,Nz)
        real*8 Syn(Nx,Ny,Nz)
        real*8 Szn(Nx,Ny,Nz)
        real*8 En(Nx,Ny,Nz)
        real*8 phiv(Nx+1,Ny+1,Nz+1)
        ! Primitive variables
        real*8 vx(Nx,Ny,Nz),vy(Nx,Ny,Nz),vz(Nx,Ny,Nz)
        real*8 u(Nx,Ny,Nz)
        real*8 P(Nx,Ny,Nz)
        real*8 x(Nx+1),y(Ny+1),z(Nz+1),dt,gamma
        real*8 x_c(Nx),y_c(Ny),z_c(Nz)
        real*8 r, PI
        integer i,j,k
        real*8 dx,dy,dz,dx2,dy2,dz2
        real*8 rho_atm, p_deplete
        real*8 v2,phil

        ! temporary.
        real*8 res(Nx+1,Ny+1,Nz+1), rhs(Nx+1,Ny+1,Nz+1)
        real*8 cmask(Nx+1,Ny+1,Nz+1)
        integer il, jl, kl

        real*8 y20l,costh
        real*8 deltarho

        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

        il = Nx / 2
        jl = Ny / 2
        kl = Nz / 2

        PI = acos(-1.d0)

        deltarho = 0.0

        ! Calculate the potential.
        do k=1,Nz+1
           do j=1,Ny+1
              do i=1,Nx+1
                 r = x(i)**2 + y(j)**2 + z(k)**2
                 r = sqrt(r)
                 
                 if (r < PI) then
                    if (r > 0.d0) then 
                       phiv(i,j,k) = -sin(r)/r
                       rhov(i,j,k) = sin(r)/r
                    else
                       phiv(i,j,k) = -1.0
                       rhov(i,j,k) = 1.0
                    end if
                 else
                    phiv(i,j,k) = 1.d0-PI/r
                    rhov(i,j,k) = rho_atm
                 end if
              end do
           end do
        end do

        ! Now do the other stuff.
        do k=1,Nz
           do j=1,Ny
              do i=1,Nx
                 r = x_c(i)**2 + y_c(j)**2 + z_c(k)**2
                 r = sqrt(r)

                 !modify r by the perturbation
                 costh = z_c(k)/r
                 y20l = 0.25d0*sqrt(5.d0/3.d0)*(3.d0*costh**2 - 1.d0)

                 if (r < PI) then
                    rho(i,j,k) = sin(r)/r
                    phil       = -sin(r)/r
                    rho(i,j,k) = rho(i,j,k) + deltarho*y20l*r/PI
                 else
                    rho(i,j,k) = rho_atm
                    phil       = 1.d0-PI/r 
                 end if

                 ! Assumes n=1
                 P(i,j,k) = p_deplete*rho(i,j,k)**2/2.d0
                 u(i,j,k) = P(i,j,k)/(gamma-1.d0)
                 en(i,j,k) = u(i,j,k) + 0.5d0*rho(i,j,k)*phil
                 sxn(i,j,k) = 0.d0
                 syn(i,j,k) = 0.d0
                 szn(i,j,k) = 0.d0
                 dn(i,j,k)  = rho(i,j,k)
                 vx(i,j,k) = 0.d0
                 vy(i,j,k) = 0.d0
                 vz(i,j,k) = 0.d0
                 
              end do
           end do
        end do

        do k=1,Nz+1
           do j=1,Ny+1
              do i=1,Nx+1
                 res(i,j,k) = 0.d0
                 rhs(i,j,k) = 0.d0
                 cmask(i,j,k) = 1.d0
              end do
           end do
        end do

        return
        end

c-----------------------------------------------------------------------
c  Find primitive variables.
c-----------------------------------------------------------------------
      subroutine find_primitives(dn,sxn,syn,
     &     szn,en,phin,rho,vx,vy,vz,u,gamma, 
     &     Nx,Ny,Nz,rho_atm,p_deplete)
        implicit none
        integer Nx,Ny,Nz
        ! Conserved variables
        real*8 dn(Nx,Ny,Nz)
        real*8 Sxn(Nx,Ny,Nz)
        real*8 Syn(Nx,Ny,Nz)
        real*8 Szn(Nx,Ny,Nz)
        real*8 En(Nx,Ny,Nz)
        real*8 phin(Nx+1,Ny+1,Nz+1)
        ! Primitive variables
        real*8 vx(Nx,Ny,Nz),vy(Nx,Ny,Nz),vz(Nx,Ny,Nz)
        real*8 u(Nx,Ny,Nz),rho(Nx,Ny,Nz)

        real*8 gamma
        real*8 u_thresh, phinl, rho_atm, p_deplete

        ! Local arrays
        real*8 P(Nx,Ny,Nz)

        integer i,j,k
        real*8 v2

        ! Achtung--coupling to gravity should be done more carefully!
        ! In particular, should you pass in a CC representation of phi
        ! which could, in general, be interpolated to higher than 
        ! second order?

        do k=1,Nz
           do j=1,Ny
              do i=1,Nx
                 ! impose density floor
                 if (dn(i,j,k) .lt. rho_atm) then
                    dn(i,j,k) = rho_atm
                    sxn(i,j,k) = 0.d0
                    syn(i,j,k) = 0.d0
                    szn(i,j,k) = 0.d0
                 end if

                 rho(i,j,k) = dn(i,j,k)
                 vx(i,j,k)  = sxn(i,j,k)/dn(i,j,k)
                 vy(i,j,k)  = syn(i,j,k)/dn(i,j,k)
                 vz(i,j,k)  = szn(i,j,k)/dn(i,j,k)
                 v2 = vx(i,j,k)*vx(i,j,k) + vy(i,j,k)*vy(i,j,k) + 
     &                vz(i,j,k)*vz(i,j,k)

                 ! This is a second order estimate of phi from the 
                 ! vertex centered function.  This is a bit stupid,
                 ! given that you have a cell centered representation
                 ! of phi (but only at certain time levels.
                 phinl = 1.d0/8.d0*(phin(i,j,k) + phin(i+1,j,k) + 
     &                phin(i,j+1,k) + phin(i,j,k+1) + 
     &                phin(i+1,j+1,k) + phin(i+1,j,k+1) + 
     &                phin(i,j+1,k+1) + phin(i+1,j+1,k+1) )

                 u(i,j,k) = en(i,j,k)-0.5d0*dn(i,j,k)*v2 -
     &                0.5d0*dn(i,j,k)*phinl

                 u_thresh = p_deplete*rho(i,j,k)**2/(gamma-1.d0)/2.d0
                 ! Impose font fix
                 if (u(i,j,k) .lt. u_thresh) then
                    u(i,j,k) = u_thresh
                    en(i,j,k) = 0.5d0*dn(i,j,k)*v2 + u(i,j,k) +
     &                   0.5d0*dn(i,j,k)*phinl
                 end if

                 ! also impose a limit on spurious heating in the
                 ! atmosphere
                 if ((gamma-1.d0)*u(i,j,k)/(rho(i,j,k))**2.gt.100) then
                    u(i,j,k) = 100.d0*(gamma-1.d0)*rho(i,j,k)**2
                    en(i,j,k) = 0.5d0*dn(i,j,k)*v2 + u(i,j,k) +
     &                   0.5d0*dn(i,j,k)*phinl
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

        m = 0.d0

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

        subroutine get_rho_c(has_origin,rho_c,rho,x,y,z,Nx,Ny,Nz,
     &     ghost_width)
        implicit none
        integer Nx,Ny,Nz,has_origin
        real*8 rho_c
        real*8 rho(Nx+1,Ny+1,Nz+1)
        real*8 x(Nx+1),y(Ny+1),z(Nz+1)
        integer ghost_width(6)
        integer i, j, k
        integer io, jo, ko
        integer origin_x 
        integer origin_y 
        integer origin_z 

        origin_x = 0
        origin_y = 0
        origin_z = 0

        do i=1+ghost_width(1),Nx+1-ghost_width(2)
           if (x(i)==0) then
              origin_x = 1
              io = i
           end if
        end do

        do i=1+ghost_width(3),Ny+1-ghost_width(4)
           if (y(i)==0) then
              origin_y = 1
              jo = i
           end if
        end do

        do i=1+ghost_width(5),Nz+1-ghost_width(6)
           if (z(i)==0) then
              origin_z = 1
              ko = i
           end if
        end do
        
        has_origin = 0
        if (origin_x>0 .and. origin_y>0 .and. origin_z>0) has_origin = 1

        if (has_origin>0) then
           rho_c = rho(io,jo,ko)
        else
           rho_c = 0.d0
        end if

        return
        end


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
         ! Set slope values at outer boundaries.  These
         ! slopes are obviously incorrect, but everything
         ! is fine as long as you have two ghost zones.
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
         ! Set slope values at outer boundaries.  These
         ! slopes are obviously incorrect, but everything
         ! is fine as long as you have two ghost zones.
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
         ! Set slope values at outer boundaries.  These
         ! slopes are obviously incorrect, but everything
         ! is fine as long as you have two ghost zones.
         do j = 1, Ny
            do i = 1, Nx
               f_slope(i,j,1)  = 0.d0
               f_slope(i,j,Nz) = 0.d0
            end do
         end do

      end if

      end subroutine mc_slopes

        subroutine init_rest(dn,dnp1,sxn,sxnp1,syn,synp1,szn,sznp1,
     &     en,enp1,phiv,phiv_np1,phi_extrap_tm1,phi_extrap_tm2,
     &     Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz
        ! Conserved variables
        real*8 dn(Nx,Ny,Nz),dnp1(Nx,Ny,Nz)
        real*8 Sxn(Nx,Ny,Nz),Sxnp1(Nx,Ny,Nz)
        real*8 Syn(Nx,Ny,Nz),Synp1(Nx,Ny,Nz)
        real*8 Szn(Nx,Ny,Nz),Sznp1(Nx,Ny,Nz)
        real*8 En(Nx,Ny,Nz),Enp1(Nx,Ny,Nz)
        real*8 phiv(Nx+1,Ny+1,Nz+1)
        real*8 phiv_np1(Nx+1,Ny+1,Nz+1)
        real*8 phi_extrap_tm1(Nx+1,Ny+1,Nz+1)
        real*8 phi_extrap_tm2(Nx+1,Ny+1,Nz+1)
        ! Primitive variables
        integer i,j,k
        integer il,jl,kl

        il = Nx / 2
        jl = Ny / 2
        kl = Nz / 2

        ! Calculate the potential.
        do k=1,Nz+1
           do j=1,Ny+1
              do i=1,Nx+1
                 phiv_np1(i,j,k) = phiv(i,j,k)
                 phi_extrap_tm1(i,j,k) = phiv(i,j,k)
                 phi_extrap_tm2(i,j,k) = phiv(i,j,k)
              end do
           end do
        end do

        ! Now do the other stuff.
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

c       Subroutine for zeroing fcs variables.  The flag determines
c       whether type a or type b cells should be zeroed.     
        subroutine zero_fcs_vars(type_flag,fc_mask,
     &     d_fcs,sx_fcs,sy_fcs,sz_fcs,e_fcs,Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz,type_flag
        ! Flux correction stuff.
        integer fc_mask(Nx,Ny,Nz)
        integer lfc_mask
        real*8 d_fcs(Nx,Ny,Nz), e_fcs(Nx,Ny,Nz)
        real*8 sx_fcs(Nx,Ny,Nz), sy_fcs(Nx,Ny,Nz), sz_fcs(Nx,Ny,Nz)
        include 'fc_mask.inc'
        real*8 lsign  ! the local sign with which fluxes enter the 
                      ! flux correction.  
        integer i,j,k

        do k=1,Nz
           do j=1,Ny
              do i=1,Nx

                 lfc_mask = fc_mask(i,j,k)
                 if (iand(lfc_mask,SIGN_FLAG) == type_flag) then 
                    if ((iand(lfc_mask,PXFLAG) .ne. 0) .or. 
     &                  (iand(lfc_mask,MXFLAG) .ne. 0) .or. 
     &                  (iand(lfc_mask,PYFLAG) .ne. 0) .or. 
     &                  (iand(lfc_mask,MYFLAG) .ne. 0) .or. 
     &                  (iand(lfc_mask,PZFLAG) .ne. 0) .or. 
     &                  (iand(lfc_mask,MZFLAG) .ne. 0)) then 
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

c       Subroutine for zeroing fcs variables.  The flag determines
c       whether type a or type b cells should be zeroed.     
      subroutine apply_fc(dn,sxn,syn,szn,en,
     &     d_fcs,sx_fcs,sy_fcs,sz_fcs,e_fcs,Nx,Ny,Nz)
      implicit none
      integer Nx,Ny,Nz
      real*8 d_fcs(Nx,Ny,Nz), e_fcs(Nx,Ny,Nz)
      real*8 sx_fcs(Nx,Ny,Nz), sy_fcs(Nx,Ny,Nz), sz_fcs(Nx,Ny,Nz)
      real*8 dn(Nx,Ny,Nz), en(Nx,Ny,Nz)
      real*8 sxn(Nx,Ny,Nz), syn(Nx,Ny,Nz), szn(Nx,Ny,Nz)
     
      ! Locals
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



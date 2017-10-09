c----------------------------------------------------------------------
c numerical routines for the Newtonian boson star code
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c initializes f with a generalized gaussian:
c
c f = A * exp (- (r-r0)^2/delta^2)
c
c where r = sqrt ( (1-ex^2)*(xb)^2 + (-ey^2)*(yb)^2 + (1-ez^2)*(zb)^2 )
c
c and xb = x-xu0, yb = ...
c----------------------------------------------------------------------
        subroutine gauss3d(f,A,r0,delta,xu0,yu0,zu0,ex,ey,ez,
     &                     x,y,z,Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz
        real*8 f(Nx,Ny,Nz),x(Nx),y(Ny),z(Nz)
        real*8 A,r0,delta,ex,ey,ez,xu0,yu0,zu0
        !TEST
        !real*8 trunc
        !parameter (trunc=0.1)

        integer i,j,k
        real*8 r,xb,yb,zb

        do i=1,Nx
           xb=x(i)-xu0
           do j=1,Ny
              yb=y(j)-yu0
              do k=1,Nz
                 zb=z(k)-zu0
                 if (abs(x(i)).eq.1.or.abs(y(j)).eq.1.or.abs(z(k)).eq.1)
     &           then
                    f(i,j,k)=0
                 else
                    r=sqrt((1-ex**2)*xb**2+
     &                     (1-ey**2)*yb**2+
     &                     (1-ez**2)*zb**2)
                    f(i,j,k)=A*exp(-((r-r0)/delta)**2) 
                    ! TEST
                    if ((r-r0).gt.delta) then
                       f(i,j,k)=0
                    else
                       f(i,j,k)=A
                    end if
                 end if
              end do
           end do
        end do

        return
        end

        !! TEST ... initial guess for V for TEST above
        subroutine initv(V,A,delta,x,y,z,Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz
        real*8 V(Nx,Ny,Nz),x(Nx),y(Ny),z(Nz)
        real*8 A,r,delta
        integer i,j,k

        do i=1,Nx
           do j=1,Ny
              do k=1,Nz
                 r=sqrt(x(i)**2+y(j)**2+z(k)**2)
                 if (r.lt.delta) r=delta
                 V(i,j,k)=-(2.0d0*delta**3*A**2)/3.0d0/r
              end do
           end do
        end do

        return 
        end


c-----------------------------------------------------------------------
c differential operator L, computed where cmask=CMASK_ON (saved in LV)
c
c L(V) = 0
c 
c =>   V,xx + V,yy + V,zz - [phi_r^2+phi_i^2] = 0
c
c with V=0 on the physical boundary
c-----------------------------------------------------------------------
        subroutine lop(LV,V,phi_r,phi_i,cmask,x,y,z,Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz
        real*8 V(Nx,Ny,Nz),phi_r(Nx,Ny,Nz),phi_i(Nx,Ny,Nz)
        real*8 cmask(Nx,Ny,Nz),LV(Nx,Ny,Nz)
        real*8 x(Nx),y(Ny),z(Nz)

        integer i,j,k
        real*8 dx,dy,dz,dx2,dy2,dz2
        real*8 v_x,v_y,v_z

        include 'cmask.inc'

        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

        dx2=dx**2
        dy2=dy**2
        dz2=dz**2

        do i=2,Nx-1
         do j=2,Ny-1
          do k=2,Nz-1
           if (cmask(i,j,k).eq.CMASK_ON) then
            v_x=(V(i+1,j,k)-V(i-1,j,k))/2/dx
            v_y=(V(i,j+1,k)-V(i,j-1,k))/2/dy
            v_z=(V(i,j,k+1)-V(i,j,k-1))/2/dz
            LV(i,j,k)=(V(i+1,j,k)-2*V(i,j,k)+V(i-1,j,k))/dx2
     &              + (V(i,j+1,k)-2*V(i,j,k)+V(i,j-1,k))/dy2
     &              + (V(i,j,k+1)-2*V(i,j,k)+V(i,j,k-1))/dz2
     &          - (phi_r(i,j,k)**2+phi_i(i,j,k)**2)
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
        subroutine residual(res,rhs,V,phi_r,phi_i,cmask,x,y,z,norm,
     &                      Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz
        real*8 V(Nx,Ny,Nz),phi_r(Nx,Ny,Nz),phi_i(Nx,Ny,Nz)
        real*8 cmask(Nx,Ny,Nz),res(Nx,Ny,Nz),rhs(Nx,Ny,Nz)
        real*8 x(Nx),y(Ny),z(Nz),norm

        integer i,j,k,sum

        include 'cmask.inc'

        call  lop(res,V,phi_r,phi_i,cmask,x,y,z,Nx,Ny,Nz)

        norm=0
        sum=0

        do i=2,Nx-1
           do j=2,Ny-1
              do k=2,Nz-1
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
        subroutine relax(V,rhs,phi_r,phi_i,cmask,phys_bdy,x,y,z,norm,
     &                   Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz
        real*8 V(Nx,Ny,Nz),phi_r(Nx,Ny,Nz),phi_i(Nx,Ny,Nz)
        real*8 cmask(Nx,Ny,Nz),rhs(Nx,Ny,Nz)
        real*8 x(Nx),y(Ny),z(Nz),norm
        integer phys_bdy(6)

        integer i,j,k,N,pass,red,sum
        integer iis,iie,jis,jie,kis,kie
        real*8 PI
        real*8 dx,dy,dz,dx2,dy2,dz2
        real*8 v_x,v_y,v_z,r

        real*8 res,Jac

        include 'cmask.inc'

        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

        dx2=dx**2
        dy2=dy**2
        dz2=dz**2

        norm=0
        sum=0

        iis=2
        iie=Nx-1
        jis=2
        jie=Ny-1
        kis=2
        kie=Nz-1

        if (phys_bdy(1).eq.1) iis=iis+1
        if (phys_bdy(2).eq.1) iie=iie-1
        if (phys_bdy(3).eq.1) jis=jis+1
        if (phys_bdy(4).eq.1) jie=jie-1
        if (phys_bdy(5).eq.1) kis=kis+1
        if (phys_bdy(6).eq.1) kie=kie-1

        do pass=0,1
         do i=2,Nx-1
          do j=2,Ny-1
           do k=2,Nz-1
            if (mod(i+j+k+pass,2).eq.0.and.
     &          (cmask(i,j,k).eq.CMASK_ON)) then
             v_x=(V(i+1,j,k)-V(i-1,j,k))/2/dx
             v_y=(V(i,j+1,k)-V(i,j-1,k))/2/dy
             v_z=(V(i,j,k+1)-V(i,j,k-1))/2/dz
             res=(V(i+1,j,k)-2*V(i,j,k)+V(i-1,j,k))/dx2
     &         + (V(i,j+1,k)-2*V(i,j,k)+V(i,j-1,k))/dy2
     &         + (V(i,j,k+1)-2*V(i,j,k)+V(i,j,k-1))/dz2
     &         - (phi_r(i,j,k)**2+phi_i(i,j,k)**2)
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
c Solve for phi via a CN discretization of the s.e., but solve
c 'simultaneously' at each point by using complex arithmetic
c-----------------------------------------------------------------------
        subroutine phi_1step_cnc(phi_rn,phi_rnp1,phi_in,phi_inp1,
     &             Vn,Vnp1,cmask,x,y,z,dt,Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz
        real*8 phi_rn(Nx,Ny,Nz),phi_rnp1(Nx,Ny,Nz)
        real*8 phi_in(Nx,Ny,Nz),phi_inp1(Nx,Ny,Nz)
        real*8 Vn(Nx,Ny,Nz),Vnp1(Nx,Ny,Nz),cmask(Nx,Ny,Nz)
        real*8 x(Nx),y(Ny),z(Nz),dt

        integer i,j,k
        real*8 dx,dy,dz,dx2,dy2,dz2

        complex *16 phi_np1

        include 'cmask.inc'

        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

        dx2=dx**2
        dy2=dy**2
        dz2=dz**2

        !--------------------------------------------------------------
        !
        ! i d(phi)/dt = -D^2(phi)/2 + V phi
        !
        !--------------------------------------------------------------

        do i=2,Nx-1
           do j=2,Ny-1
              do k=2,Nz-1
               if (cmask(i,j,k).eq.CMASK_ON) then
                  phi_np1=
     &        ( dcmplx(0,1/dt)*dcmplx(phi_rn(i,j,k),phi_in(i,j,k)) 
     &         -dcmplx(0.25d0/dx2,0)*
     &          dcmplx(phi_rn(i+1,j,k)-2*phi_rn(i,j,k)+phi_rn(i-1,j,k)+
     &                 phi_rnp1(i+1,j,k)+phi_rnp1(i-1,j,k),
     &                 phi_in(i+1,j,k)-2*phi_in(i,j,k)+phi_in(i-1,j,k)+
     &                 phi_inp1(i+1,j,k)+phi_inp1(i-1,j,k))
     &          
     &         -dcmplx(0.25d0/dy2,0)*
     &          dcmplx(phi_rn(i,j+1,k)-2*phi_rn(i,j,k)+phi_rn(i,j-1,k)+
     &                 phi_rnp1(i,j+1,k)+phi_rnp1(i,j-1,k),
     &                 phi_in(i,j+1,k)-2*phi_in(i,j,k)+phi_in(i,j-1,k)+
     &                 phi_inp1(i,j+1,k)+phi_inp1(i,j-1,k))
     &          
     &         -dcmplx(0.25d0/dz2,0)*
     &          dcmplx(phi_rn(i,j,k+1)-2*phi_rn(i,j,k)+phi_rn(i,j,k-1)+
     &                 phi_rnp1(i,j,k+1)+phi_rnp1(i,j,k-1),
     &                 phi_in(i,j,k+1)-2*phi_in(i,j,k)+phi_in(i,j,k-1)+
     &                 phi_inp1(i,j,k+1)+phi_inp1(i,j,k-1))
     &           
     &         +dcmplx((Vnp1(i,j,k)+Vn(i,j,k))/4,0)*
     &          dcmplx(phi_rn(i,j,k),phi_in(i,j,k)))
     &         /
     &          dcmplx(-0.5d0/dx2- 0.5d0/dy2 - 0.5d0/dz2
     &                 -(Vnp1(i,j,k)+Vn(i,j,k))/4,1/dt)
                  phi_rnp1(i,j,k)=dreal(phi_np1)
                  phi_inp1(i,j,k)=dimag(phi_np1)
               end if
              end do
           end do
        end do

        return
        end

c-----------------------------------------------------------------------
c the following computes a residual of the Schroedinger equation
c (based upon the CN scheme) and returns the L2 norm of it 
c-----------------------------------------------------------------------
        subroutine se_ires(res,l2norm,phi_rn,phi_rnp1,phi_in,phi_inp1,
     &             Vn,Vnp1,cmask,x,y,z,dt,Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz
        real*8 phi_rn(Nx,Ny,Nz),phi_rnp1(Nx,Ny,Nz)
        real*8 phi_in(Nx,Ny,Nz),phi_inp1(Nx,Ny,Nz)
        real*8 Vn(Nx,Ny,Nz),Vnp1(Nx,Ny,Nz),res(Nx,Ny,Nz)
        real*8 x(Nx),y(Ny),z(Nz),dt,cmask(Nx,Ny,Nz),l2norm

        integer i,j,k,sum,N
        real*8 norm,resr,resi
        real*8 dx,dy,dz,dx2,dy2,dz2

        include 'cmask.inc'

        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

        dx2=dx**2
        dy2=dy**2
        dz2=dz**2

        norm=0
        sum=0

        do i=2,Nx-1
           do j=2,Ny-1
              do k=2,Nz-1
                 if (cmask(i,j,k).eq.CMASK_on) then
                    resr=(phi_rnp1(i,j,k)-phi_rn(i,j,k))/dt+
     &         (0.25d0/dx2*
     &          (phi_in(i+1,j,k)-2*phi_in(i,j,k)+phi_in(i-1,j,k)+
     &           phi_inp1(i+1,j,k)-2*phi_inp1(i,j,k)+phi_inp1(i-1,j,k)))
     &         +
     &         (0.25d0/dy2*
     &          (phi_in(i,j+1,k)-2*phi_in(i,j,k)+phi_in(i,j-1,k)+
     &           phi_inp1(i,j+1,k)-2*phi_inp1(i,j,k)+phi_inp1(i,j-1,k)))
     &         +
     &         (0.25d0/dz2*
     &          (phi_in(i,j,k+1)-2*phi_in(i,j,k)+phi_in(i,j,k-1)+
     &           phi_inp1(i,j,k+1)-2*phi_inp1(i,j,k)+phi_inp1(i,j,k-1)))
     &         - 
     &           (Vnp1(i,j,k)+Vn(i,j,k))/2*
     &           (phi_in(i,j,k)+phi_inp1(i,j,k))/2

                    resi=(phi_inp1(i,j,k)-phi_in(i,j,k))/dt-
     &         (0.25d0/dx2*
     &          (phi_rn(i+1,j,k)-2*phi_rn(i,j,k)+phi_rn(i-1,j,k)+
     &           phi_rnp1(i+1,j,k)-2*phi_rnp1(i,j,k)+phi_rnp1(i-1,j,k)))
     &         -
     &         (0.25d0/dy2*
     &          (phi_rn(i,j+1,k)-2*phi_rn(i,j,k)+phi_rn(i,j-1,k)+
     &           phi_rnp1(i,j+1,k)-2*phi_rnp1(i,j,k)+phi_rnp1(i,j-1,k)))
     &         -
     &         (0.25d0/dz2*
     &          (phi_rn(i,j,k+1)-2*phi_rn(i,j,k)+phi_rn(i,j,k-1)+
     &           phi_rnp1(i,j,k+1)-2*phi_rnp1(i,j,k)+phi_rnp1(i,j,k-1)))
     &         + 
     &           (Vnp1(i,j,k)+Vn(i,j,k))/2*
     &           (phi_rn(i,j,k)+phi_rnp1(i,j,k))/2

                    res(i,j,k)=sqrt(resi**2+resr**2)
                  
                    norm=norm+res(i,j,k)**2
                    sum=sum+1
                 end if
              end do
           end do
        end do

        l2norm=sqrt(norm/sum)

        return 
        end

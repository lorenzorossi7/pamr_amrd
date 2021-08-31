c---------------------------------------------------------------------------
c 'real *8' Math utility routines ... for copying, interpolation, 
c restrictions, etc.
c
c some functions borrowed from dmatlib.f
c
c For versions of the routines supporting excision (_ex suffix),
c the characteristic masks (chr1,chr2,etc) at different resolutions are 
c expected to be discrete representations of the same anayltic excision 
c function, i.e. matching coarse-fine points must agree on being 
c excised (as defined by values of chr==ex)
c---------------------------------------------------------------------------

c---------------------------------------------------------------------------
c dmcopy3d copies a segment of grid u1 into grid u2.
c offset in u1 = [i1,j1,k1] with stride [si1,sj1,sk1],
c offset in u2 = [i2,j2,k2] with stride [si2,sj2,sk2]
c size of segment=ni x nj x nk
c---------------------------------------------------------------------------
        subroutine dmcopy3d(u1,u2,nx1,ny1,nz1,nx2,ny2,nz2,
     &                      i1,j1,k1,i2,j2,k2,ni,nj,nk,
     &                      si1,sj1,sk1,si2,sj2,sk2)
        implicit none
        integer nx1,ny1,nz1,nx2,ny2,nz2
        integer i1,j1,k1,i2,j2,k2,ni,nj,nk
        integer si1,sj1,sk1,si2,sj2,sk2
        real*8 u1(nx1,ny1,nz1),u2(nx2,ny2,nz2)

        integer i,j,k
 
        if (i1.lt.1.or.(i1+(ni-1)*si1).gt.nx1) return
        if (i2.lt.1.or.(i2+(ni-1)*si2).gt.nx2) return

        if (j1.lt.1.or.(j1+(nj-1)*sj1).gt.ny1) return
        if (j2.lt.1.or.(j2+(nj-1)*sj2).gt.ny2) return

        if (k1.lt.1.or.(k1+(nk-1)*sk1).gt.nz1) return
        if (k2.lt.1.or.(k2+(nk-1)*sk2).gt.nz2) return

        do i=1,ni
           do j=1,nj
              do k=1,nk
                 u2(i2+(i-1)*si2,j2+(j-1)*sj2,k2+(k-1)*sk2)=
     &           u1(i1+(i-1)*si1,j1+(j-1)*sj1,k1+(k-1)*sk1)
              end do
           end do
        end do

        return
        end

c-----------------------------------------------------------------------
c Performs a irho:1 cubic interpolation of v1(n1) into v2(n2) 
c (based on similar routine from dveclib).
c
c interpolates v1[i1,i1+s1,...,i1+n*s1] into
c              v2[i2,i2+s2,...,i2+n*rho*s2]
c
c if 1<=i1<=s1, then an interior stencil is used at the beginning,
c similarly if i1+n*s1<=(n1-s1)
c
c NOTES: 1. no bounds checking done here
c        2. if n<4 then dvi2q1d is called.
c-----------------------------------------------------------------------
        subroutine dvi4q1d(v1,v2,i1,i2,n1,n2,n,irho,s1,s2)
        implicit none
        integer i1,i2,n1,n2,n,irho,s1,s2
        real*8 v1(n1),v2(n2)
 
        real*8 c0,c1,c2,c3,sg
        integer i,j,is,ie
 
        real*8 half,mhalf,sixth,msixth
        parameter  ( half   =  0.5000 0000 0000 0000 d0,
     &               mhalf  = -0.5000 0000 0000 0000 d0,
     &               sixth  =  0.1666 6666 6666 6667 d0,
     &               msixth = -0.1666 6666 6666 6667 d0  )

        if (n.lt.4) then
           call dvi2q1d(v1,v2,i1,i2,n1,n2,n,irho,s1,s2)
           return
        end if
 
        !--------------------------------------------------------------
        ! common points
        !--------------------------------------------------------------
        do i=1,n
           v2(irho*(i-1)*s2+i2)=v1((i-1)*s1+i1)
        end do

        !--------------------------------------------------------------
        ! left boundary
        !--------------------------------------------------------------
        if (i1.le.s1) then
           do j=1,irho-1
              sg = (j*1.0d0)/irho
              c0 = msixth*(sg-1.0d0)*(sg-2.0d0)*(sg-3.0d0)
              c1 = half*sg*(sg-2.0d0)*(sg-3.0d0)
              c2 = mhalf*sg*(sg-1.0d0)*(sg-3.0d0)
              c3 = sixth*sg*(sg-1.0d0)*(sg-2.0d0)
              v2(j*s2+i2)=c0*v1(i1)+c1*v1(s1+i1)+
     &                    c2*v1(2*s1+i1)+c3*v1(3*s1+i1)
           end do
           is=1
        else
           is=0
        end if
         
        !--------------------------------------------------------------
        ! right boundary
        !--------------------------------------------------------------
        if ((i1+n*s1).gt.(n1-2*s1)) then
           ie=n-3
           do j=1,irho-1
              sg = 2.0d0+(j*1.0d0)/irho
              c0 = msixth*(sg-1.0d0)*(sg-2.0d0)*(sg-3.0d0)
              c1 = half*sg*(sg-2.0d0)*(sg-3.0d0)
              c2 = mhalf*sg*(sg-1.0d0)*(sg-3.0d0)
              c3 = sixth*sg*(sg-1.0d0)*(sg-2.0d0)
              v2((irho*(n-2)+j)*s2+i2)=
     &           c0*v1((n-4)*s1+i1)+c1*v1((n-3)*s1+i1)+
     &           c2*v1((n-2)*s1+i1)+c3*v1((n-1)*s1+i1)
           end do
        else
           ie=n-2
        end if

        !--------------------------------------------------------------
        ! interior
        !--------------------------------------------------------------
        do j=1,irho-1
           sg = 1.0d0+(j*1.0d0)/irho
           c0 = msixth*(sg-1.0d0)*(sg-2.0d0)*(sg-3.0d0)
           c1 = half*sg*(sg-2.0d0)*(sg-3.0d0)
           c2 = mhalf*sg*(sg-1.0d0)*(sg-3.0d0)
           c3 = sixth*sg*(sg-1.0d0)*(sg-2.0d0)
           do i=is,ie
              v2((irho*i+j)*s2+i2)=c0*v1((i-1)*s1+i1)+c1*v1(i*s1+i1)+
     &                             c2*v1((i+1)*s1+i1)+c3*v1((i+2)*s1+i1)
           end do
        end do

        return
        end

c-----------------------------------------------------------------------
c version of dvi4q1d with excision ... if *any* non-excised segment
c of v1 (coarse) is less that 4 in width then passed to linear routine
c
c NOTE: this function overwrite's portions of v1 and v2's excised zone!
c the [3] vars are associated with chr1, v1's mask
c-----------------------------------------------------------------------
        subroutine dvi4q1d_ex(v1,v2,i1,i2,n1,n2,n,irho,s1,s2,
     &                        chr1,i3,n3,s3,ex)
        implicit none
        integer i1,i2,i3,n1,n2,n3,n,irho,s1,s2,s3
        real*8 v1(n1),v2(n2)
        real*8 chr1(n3),ex
 
        real*8 c0(20),c1(20),c2(20),c3(20),sg
        integer i,j,is,ie,cn,j0,i1s,i1e,i10,i30
        logical init_2
 
        real*8 half,mhalf,sixth,msixth
        parameter  ( half   =  0.5000 0000 0000 0000 d0,
     &               mhalf  = -0.5000 0000 0000 0000 d0,
     &               sixth  =  0.1666 6666 6666 6667 d0,
     &               msixth = -0.1666 6666 6666 6667 d0  )

        if (irho.gt.20) then
           write(*,*) 'dvi4q1d_ex : irho expected to be less than 20'
           stop
        end if

        !--------------------------------------------------------------
        ! because we can interpolate from outside the coarse grid range,
        ! we also check wether these points on the coarse grid are
        ! unexcised
        !--------------------------------------------------------------
        i=3
        i1s=i1
        do while(i.gt.0.and.i1s.gt.s1)
           i1s=i1s-s1
           i=i-1
        end do
        i=3
        i1e=i1+(n-1)*s1
        do while(i.gt.0.and.i1e.lt.(n1-s1))
           i1e=i1e+s1
           i=i-1
        end do

        cn=0 
        init_2=.false.
        do i=i1s,i1e,s1
           if (chr1((i-i1)*s3/s1+i3).eq.ex.or.i.eq.i1e) then
              if (cn.gt.0.and.cn.lt.4) then
                 init_2=.true.
              else
                 cn=0
              end if
           else
              cn=cn+1
           end if
        end do
        if (init_2) call dvi2q1d_ex(v1,v2,i1,i2,n1,n2,n,irho,s1,s2,
     &                              chr1,i3,n3,s3,ex)

        !--------------------------------------------------------------
        ! first, common points, and extrapolate all coarse grid points 
        ! next to excision boundaries by 1 point.
        !--------------------------------------------------------------
        do i=i1s,i1e,s1
!LR: equivalent to the commented out piece below, but it does not go out of bounds
           if (chr1(i3+s3*(i-i1)/s1).eq.ex) then
              if (i.ge.(i1s+3*s1)) then
               if (chr1(i3+s3*(i-s1-i1)/s1).ne.ex.and.
     &            chr1(i3+s3*(i-2*s1-i1)/s1).ne.ex.and.
     &            chr1(i3+s3*(i-3*s1-i1)/s1).ne.ex) then
                 if (i.ge.(i1s+4*s1)) then
                  if (chr1(i3+s3*(i-4*s1-i1)/s1).ne.ex) then
                   v1(i)=4*(v1(i-s1)+v1(i-3*s1))-6*v1(i-2*s1)-v1(i-4*s1)
                  end if
                 else
                  v1(i)=3*(v1(i-s1)-v1(i-2*s1))+v1(i-3*s1)
                 end if
               end if
              else if (i.le.(i1e-3*s1)) then
                if (chr1(i3+s3*(i+s1-i1)/s1).ne.ex.and.
     &                 chr1(i3+s3*(i+2*s1-i1)/s1).ne.ex.and.
     &                 chr1(i3+s3*(i+3*s1-i1)/s1).ne.ex) then
                 if (i.le.(i1s-4*s1)) then
                  if (chr1(i3+s3*(i+4*s1-i1)/s1).ne.ex) then
                   v1(i)=4*(v1(i+s1)+v1(i+3*s1))-6*v1(i+2*s1)-v1(i+4*s1)
                  end if
                 else
                  v1(i)=3*(v1(i+s1)-v1(i+2*s1))+v1(i+3*s1)
                 end if
                end if
              end if
           end if
!LR commented this out
!           if (chr1(i3+s3*(i-i1)/s1).eq.ex) then
!              if (i.ge.(i1s+3*s1).and.
!     &            chr1(i3+s3*(i-s1-i1)/s1).ne.ex.and.
!     &            chr1(i3+s3*(i-2*s1-i1)/s1).ne.ex.and.
!     &            chr1(i3+s3*(i-3*s1-i1)/s1).ne.ex) then
!                 if (i.ge.(i1s+4*s1).and.
!     &              chr1(i3+s3*(i-4*s1-i1)/s1).ne.ex) then
!                  v1(i)=4*(v1(i-s1)+v1(i-3*s1))-6*v1(i-2*s1)-v1(i-4*s1)
!                 else
!                  v1(i)=3*(v1(i-s1)-v1(i-2*s1))+v1(i-3*s1)
!                 end if
!              else if (i.le.(i1e-3*s1).and.
!     &                 chr1(i3+s3*(i+s1-i1)/s1).ne.ex.and.
!     &                 chr1(i3+s3*(i+2*s1-i1)/s1).ne.ex.and.
!     &                 chr1(i3+s3*(i+3*s1-i1)/s1).ne.ex) then
!                 if (i.le.(i1s-4*s1).and.
!     &               chr1(i3+s3*(i+4*s1-i1)/s1).ne.ex) then
!                  v1(i)=4*(v1(i+s1)+v1(i+3*s1))-6*v1(i+2*s1)-v1(i+4*s1)
!                 else
!                  v1(i)=3*(v1(i+s1)-v1(i+2*s1))+v1(i+3*s1)
!                 end if
!              end if
!           end if
           if (i.ge.i1.and.i.le.(i1+(n-1)*s1)) 
     &        v2(irho*((i-i1)/s1)*s2+i2)=v1(i)
        end do

        !--------------------------------------------------------------
        ! left boundaries cells
        !--------------------------------------------------------------
        do j=1,irho-1
           sg = (j*1.0d0)/irho
           c0(j) = msixth*(sg-1.0d0)*(sg-2.0d0)*(sg-3.0d0)
           c1(j) = half*sg*(sg-2.0d0)*(sg-3.0d0)
           c2(j) = mhalf*sg*(sg-1.0d0)*(sg-3.0d0)
           c3(j) = sixth*sg*(sg-1.0d0)*(sg-2.0d0)
        end do
        do i=1,n-1
           i10=(i-1)*s1+i1
           i30=(i-1)*s3+i3
           j0=(i-1)*irho*s2+i2
!LR: equivalent to the commented out piece below, but it does not go out of bounds      
           if (i10.le.(i1e-3*s1)) then
            if (chr1(i30+s3).ne.ex.and.
     &          chr1(i30+2*s3).ne.ex) then
             if (i.eq.1.and.i10.eq.i1s) then
              do j=1,irho-1
                 v2(j*s2+j0)=c0(j)*v1(i10)+c1(j)*v1(s1+i10)+
     &                       c2(j)*v1(2*s1+i10)+c3(j)*v1(3*s1+i10)
              end do
             else if (i10.ge.(i1s+s1).and.
     &           chr1(i30).ne.ex.and.
     &           chr1(i30-s3).eq.ex) then
              do j=1,irho-1
                 v2(j*s2+j0)=c0(j)*v1(i10)+c1(j)*v1(s1+i10)+
     &                       c2(j)*v1(2*s1+i10)+c3(j)*v1(3*s1+i10)
              end do     
             else if (chr1(i30).eq.ex) then
              do j=1,irho-1
                 v2(j*s2+j0)=c0(j)*v1(i10)+c1(j)*v1(s1+i10)+
     &                       c2(j)*v1(2*s1+i10)+c3(j)*v1(3*s1+i10)
              end do
             end if
            end if
           end if
!LR commented this out
!           if (i10.le.(i1e-3*s1).and.chr1(i30+s3).ne.ex.and.
!     &         chr1(i30+2*s3).ne.ex.and.
!     &         ((i.eq.1.and.i10.eq.i1s).or.
!     &          (i10.ge.(i1s+s1).and.chr1(i30).ne.ex.and.
!     &           chr1(i30-s3).eq.ex).or.
!     &           chr1(i30).eq.ex)) then
!              do j=1,irho-1
!                 v2(j*s2+j0)=c0(j)*v1(i10)+c1(j)*v1(s1+i10)+
!     &                       c2(j)*v1(2*s1+i10)+c3(j)*v1(3*s1+i10)
!              end do
!           end if
        end do
         
        !--------------------------------------------------------------
        ! right boundary cells
        !--------------------------------------------------------------
        do j=1,irho-1
           sg = 2.0d0+(j*1.0d0)/irho
           c0(j) = msixth*(sg-1.0d0)*(sg-2.0d0)*(sg-3.0d0)
           c1(j) = half*sg*(sg-2.0d0)*(sg-3.0d0)
           c2(j) = mhalf*sg*(sg-1.0d0)*(sg-3.0d0)
           c3(j) = sixth*sg*(sg-1.0d0)*(sg-2.0d0)
        end do
        do i=1,n-1
           i10=(i-1)*s1+i1
           i30=(i-1)*s3+i3
           j0=(i-1)*irho*s2+i2
!LR: equivalent to the commented out piece below, but it does not go out of bounds
           if (i10.ge.(i1s+2*s1)) then
            if (chr1(i30-s3).ne.ex) then
             if (i.eq.(n-1).and.i10.eq.(i1e-s1)) then
              do j=1,irho-1
                 v2(j*s2+j0)=
     &              c0(j)*v1(i10-2*s1)+c1(j)*v1(i10-s1)+
     &              c2(j)*v1(i10)+c3(j)*v1(i10+s1)
              end do
             else if (i10.le.(i1e-s1).and.
     &          chr1(i30).ne.ex.and.
     &          chr1(i30+s3).eq.ex) then
              do j=1,irho-1
                 v2(j*s2+j0)=
     &              c0(j)*v1(i10-2*s1)+c1(j)*v1(i10-s1)+
     &              c2(j)*v1(i10)+c3(j)*v1(i10+s1)
              end do
             end if
            end if
           end if
!LR commented this out
!           if (i10.ge.(i1s+2*s1).and.chr1(i30-s3).ne.ex.and.
!     &        ((i.eq.(n-1).and.i10.eq.(i1e-s1)).or.
!     &         (i10.le.(i1e-s1).and.
!     &          chr1(i30).ne.ex.and.chr1(i30+s3).eq.ex))) then
!              do j=1,irho-1
!                 v2(j*s2+j0)=
!     &              c0(j)*v1(i10-2*s1)+c1(j)*v1(i10-s1)+
!     &              c2(j)*v1(i10)+c3(j)*v1(i10+s1)
!              end do
!           end if
        end do

        !--------------------------------------------------------------
        ! interior
        !--------------------------------------------------------------
        do j=1,irho-1
           sg = 1.0d0+(j*1.0d0)/irho
           c0(j) = msixth*(sg-1.0d0)*(sg-2.0d0)*(sg-3.0d0)
           c1(j) = half*sg*(sg-2.0d0)*(sg-3.0d0)
           c2(j) = mhalf*sg*(sg-1.0d0)*(sg-3.0d0)
           c3(j) = sixth*sg*(sg-1.0d0)*(sg-2.0d0)
        end do
        do i=1,n
           i10=(i-1)*s1+i1
           i30=(i-1)*s3+i3
           j0=(i-1)*irho*s2+i2
!LR: equivalent to the commented out piece below, but it does not go out of bounds
           if (i10.ge.(i1s+s1).and.i10.le.(i1e-2*s1)) then
            if (chr1(i30).ne.ex.and.
     &          chr1(i30-s3).ne.ex.and.
     &          chr1(i30+s3).ne.ex) then
              do j=1,irho-1
                 v2(j0+j*s2)=c0(j)*v1(i10-s1)+c1(j)*v1(i10)+
     &                       c2(j)*v1(i10+s1)+c3(j)*v1(i10+2*s1)
              end do
            end if
           end if
!LR commented this out 
!           if (chr1(i30).ne.ex.and.
!     &         i10.ge.(i1s+s1).and.chr1(i30-s3).ne.ex.and.
!     &         i10.le.(i1e-2*s1).and.chr1(i30+s3).ne.ex) then
!              do j=1,irho-1
!                 v2(j0+j*s2)=c0(j)*v1(i10-s1)+c1(j)*v1(i10)+
!     &                       c2(j)*v1(i10+s1)+c3(j)*v1(i10+2*s1)
!              end do
!           end if
        end do

        return
        end


c-----------------------------------------------------------------------
c 2nd order version of dvi4q1d
c Performs a irho:1 linear interpolation of v1(n1) into v2(n2) 
c
c interpolates v1[i1,i1+s1,...,i1+n*s1] into
c              v2[i2,i2+s2,...,i2+n*rho*s2]
c
c if 1<=i1<=s1, then an interior stencil is used at the beginning,
c similarly if i1+n*s1<=(n1-s1)
c-----------------------------------------------------------------------
        subroutine dvi2q1d(v1,v2,i1,i2,n1,n2,n,irho,s1,s2)
        implicit none
        integer i1,i2,n1,n2,n,irho,s1,s2
        real*8 v1(n1),v2(n2)
 
        real*8 c0,c1,sg
        integer i,j,is,ie
 
        !--------------------------------------------------------------
        ! common points
        !--------------------------------------------------------------
        do i=1,n
           v2(irho*(i-1)*s2+i2)=v1((i-1)*s1+i1)
        end do

        !--------------------------------------------------------------
        ! interior
        !--------------------------------------------------------------
        do j=1,irho-1
           sg = (j*1.0d0)/irho
           c0 = 1-sg
           c1 = sg
           do i=0,n-2
              v2((irho*i+j)*s2+i2)=c0*v1(i*s1+i1)+c1*v1((i+1)*s1+i1)
           end do
        end do

        return
        end

c-----------------------------------------------------------------------
c version of dvi2q1d with excision 
c
c NOTE: this function overwrites some of v1 and v2's excised zone!!
c the [3] vars are associated with chr1, v1's mask.
c-----------------------------------------------------------------------
        subroutine dvi2q1d_ex(v1,v2,i1,i2,n1,n2,n,irho,s1,s2,
     &                        chr1,i3,n3,s3,ex)
        implicit none
        integer i1,i2,i3,n1,n2,n3,n,irho,s1,s2,s3
        real*8 v1(n1),v2(n2)
        real*8 chr1(n3),ex
 
        real*8 c0,c1,sg
        integer i,j,is,ie,i10,i30
 
        !--------------------------------------------------------------
        ! common points, and extrapolate all coarse points by 1 cell
        ! into the excised zone
        !--------------------------------------------------------------
        do i=1,n
           i10=(i-1)*s1+i1
           i30=(i-1)*s3+i3
!LR: equivalent to the commented out piece below, but it does not go out of bounds
           if (chr1(i30).eq.ex) then
              if (i10.gt.3*s1) then
               if (chr1(i30-s3).ne.ex.and.
     &            chr1(i30-2*s3).ne.ex.and.chr1(i30-3*s3).ne.ex) then
                 v1(i10)=3*(v1(i10-s1)-v1(i10-2*s1))+v1(i10-3*s1)
               end if
              else if (i10.le.(n1-3*s1)) then
               if (chr1(i30+s3).ne.ex.and.
     &           chr1(i30+2*s3).ne.ex.and.chr1(i30+3*s3).ne.ex) then
                 v1(i10)=3*(v1(i10+s1)-v1(i10+2*s1))+v1(i10+3*s1)
               end if
              else if (i10.gt.2*s1) then
               if (chr1(i30-s3).ne.ex.and.
     &             chr1(i30-2*s3).ne.ex) then
                 v1(i10)=2*v1(i10-s1)-v1(i10-2*s1)
               end if
              else if (i10.le.(n1-2*s1)) then
               if (chr1(i30+s3).ne.ex.and.
     &             chr1(i30+2*s3).ne.ex) then
                 v1(i10)=2*v1(i10+s1)-v1(i10+2*s1)
               end if
              else if (i10.gt.s1) then
               if (chr1(i30-s3).ne.ex) then
                 v1(i10)=v1(i10-s1)
               end if
              else if (i10.le.(n1-s1)) then
               if (chr1(i30+s3).ne.ex) then
                 v1(i10)=v1(i10+s1)
               end if
              end if
           end if
!LR commented this out
!           if (chr1(i30).eq.ex) then
!              if (i10.gt.3*s1.and.chr1(i30-s3).ne.ex.and.
!     &            chr1(i30-2*s3).ne.ex.and.chr1(i30-3*s3).ne.ex) then
!                 v1(i10)=3*(v1(i10-s1)-v1(i10-2*s1))+v1(i10-3*s1)
!              else if (i10.le.(n1-3*s1).and.chr1(i30+s3).ne.ex.and.
!     &           chr1(i30+2*s3).ne.ex.and.chr1(i30+3*s3).ne.ex) then
!                 v1(i10)=3*(v1(i10+s1)-v1(i10+2*s1))+v1(i10+3*s1)
!              else if (i10.gt.2*s1.and.chr1(i30-s3).ne.ex.and.
!     &                            chr1(i30-2*s3).ne.ex) then
!                 v1(i10)=2*v1(i10-s1)-v1(i10-2*s1)
!              else if (i10.le.(n1-2*s1).and.chr1(i30+s3).ne.ex.and.
!     &                                      chr1(i30+2*s3).ne.ex) then
!                 v1(i10)=2*v1(i10+s1)-v1(i10+2*s1)
!              else if (i10.gt.s1.and.chr1(i30-s3).ne.ex) then
!                 v1(i10)=v1(i10-s1)
!              else if (i10.le.(n1-s1).and.chr1(i30+s3).ne.ex) then
!                 v1(i10)=v1(i10+s1)
!              end if
!           end if
           v2(irho*(i-1)*s2+i2)=v1(i10)
        end do

        !--------------------------------------------------------------
        ! interior
        !--------------------------------------------------------------
        do j=1,irho-1
           sg = (j*1.0d0)/irho
           c0 = 1-sg
           c1 = sg
           do i=0,n-2
              i10=(irho*i+j)*s2+i2
              v2(i10)=c0*v1(i*s1+i1)+c1*v1((i+1)*s1+i1)
           end do
        end do

        return
        end

c---------------------------------------------------------------------------
c utility routines called by dminterp3d
c---------------------------------------------------------------------------
        subroutine init_chrt(chrt,chr1,n,n3,s3,ex)
        implicit none
        integer n3,s3,n
        real*8 chr1(n3),chrt(n),ex
 
        integer i

        do i=1,n
           chrt(i)=chr1((i-1)*s3+1)
        end do

        return
        end

        subroutine and_chrt(chrt,chr1,n,n3,s3,ex)
        implicit none
        integer n3,s3,n
        real*8 chr1(n3),chrt(n),ex
 
        integer i
         
        do i=1,n
           if ((chrt(i).eq.ex).and.(chr1((i-1)*s3+1).eq.ex)) then
              chrt(i)=ex
           else
              chrt(i)=ex+1
           end if
        end do

        return
        end

c---------------------------------------------------------------------------
c dminterp3d interpolates a segment of grid u1 into grid u2
c
c offset in u1 = [i1,j1,k1] with stride [si1,sj1,sk1] (nominally [1,1,1])
c offset in u2 = [i2,j2,k2] with rho    [irho,jrho,krho] 
c size of segment=ni x nj x nk in u1
c
c ord is the interpolation order
c
c 2 and 4 currently supported
c
c if (ex==1) then the _ex versions of the interpolation functions are called
c
c NOTE: portions of u1 and u2's excised zones will be overwritten!!
c---------------------------------------------------------------------------
        subroutine dminterp3d(u1,u2,nx1,ny1,nz1,nx2,ny2,nz2,
     &                       i1,j1,k1,i2,j2,k2,ni,nj,nk,
     &                       si1,sj1,sk1,irho,jrho,krho,ord,
     &                       chr1,ex,do_ex)
        implicit none
        integer nx1,ny1,nz1,nx2,ny2,nz2
        integer i1,j1,k1,i2,j2,k2,ni,nj,nk
        integer si1,sj1,sk1,irho,jrho,krho,ord,do_ex
        real*8 u1(nx1,ny1,nz1),u2(nx2,ny2,nz2)
        real*8 chr1(nx1,ny1,nz1),ex

        integer i,j,k,nmax
        parameter (nmax=3000)
        real*8 chrt(nmax)

        if (i1.lt.1.or.(i1+ni-1)*si1.gt.nx1) return
        if (i2.lt.1.or.(i2+(ni-1)*irho).gt.nx2) return

        if (j1.lt.1.or.(j1+nj-1)*sj1.gt.ny1) return
        if (j2.lt.1.or.(j2+(nj-1)*jrho).gt.ny2) return

        if (k1.lt.1.or.(k1+nk-1)*sk1.gt.nz1) return
        if (k2.lt.1.or.(k2+(nk-1)*krho).gt.nz2) return

        if (.not.(ord.eq.2.or.ord.eq.4)) then
           write (*,*) 'dminterp3d: only ord=2,4 currently supported' 
           stop
        end if 

        if (do_ex.eq.1.and.(nx1.gt.nmax.or.ny1.gt.nmax.or.nz1.gt.nmax))
     &     then
           write(*,*) 'dminterp3d: nmax too small'
           stop
        end if

        !--------------------------------------------------------------
        ! copy common points
        !--------------------------------------------------------------
        call dmcopy3d(u1,u2,nx1,ny1,nz1,nx2,ny2,nz2,
     &                i1,j1,k1,i2,j2,k2,ni,nj,nk,
     &                si1,sj1,sk1,irho,jrho,krho)

        !--------------------------------------------------------------
        ! interpolate i line segments with anchor points
        !--------------------------------------------------------------
        do j=j2,j2+(nj-1)*jrho,jrho
           do k=k2,k2+(nk-1)*krho,krho
              if (ord.eq.2) then
                 if (do_ex.eq.1) then
                    call dvi2q1d_ex(u2(1,j,k),u2(1,j,k),i2,i2,
     &                   nx2,nx2,ni,irho,irho,1,
     &                   chr1(1,j1+sj1*(j-j2)/jrho,k1+sk1*(k-k2)/krho),
     &                   i1,nx1,si1,ex)
                 else
                    call dvi2q1d(u2(1,j,k),u2(1,j,k),i2,i2,
     &                   nx2,nx2,ni,irho,irho,1)
                 end if
              else
                 if (do_ex.eq.1) then
                    call dvi4q1d_ex(u2(1,j,k),u2(1,j,k),i2,i2,
     &                   nx2,nx2,ni,irho,irho,1,
     &                   chr1(1,j1+sj1*(j-j2)/jrho,k1+sk1*(k-k2)/krho),
     &                   i1,nx1,si1,ex)
                 else
                    call dvi4q1d(u2(1,j,k),u2(1,j,k),i2,i2,
     &                   nx2,nx2,ni,irho,irho,1)
                 end if
              end if
           end do
        end do

        if (nj.eq.1.and.nk.eq.1) return

        !--------------------------------------------------------------
        ! now all j line segments, using points filled in by i above
        !--------------------------------------------------------------
        do i=i2,i2+(ni-1)*irho
           do k=k2,k2+(nk-1)*krho,krho
              if (do_ex.eq.1) then
                 call init_chrt(chrt,
     &              chr1(i1+si1*((i-i2)/irho),1,k1+sk1*((k-k2)/krho)),
     &              ny1,nx1*(ny1-1)+1,nx1,ex)
                 if (mod(i-i2,irho).ne.0) 
     &              call and_chrt(chrt,
     &              chr1(i1+si1*((i-i2)/irho+1),1,k1+sk1*((k-k2)/krho)),
     &              ny1,nx1*(ny1-1)+1,nx1,ex)
              end if
              if (ord.eq.2) then
                 if (do_ex.eq.1) then
                    call dvi2q1d_ex(u2(i,1,k),u2(i,1,k),1+(j2-1)*nx2,
     &                   1+(j2-1)*nx2,nx2*(ny2-1)+1,nx2*(ny2-1)+1,
     &                   nj,jrho,jrho*nx2,nx2,chrt,j1,ny1,sj1,ex)
                 else
                    call dvi2q1d(u2(i,1,k),u2(i,1,k),1+(j2-1)*nx2,
     &                   1+(j2-1)*nx2,nx2*(ny2-1)+1,nx2*(ny2-1)+1,
     &                   nj,jrho,jrho*nx2,nx2)
                 end if
              else
                 if (do_ex.eq.1) then
                    call dvi4q1d_ex(u2(i,1,k),u2(i,1,k),1+(j2-1)*nx2,
     &                   1+(j2-1)*nx2,nx2*(ny2-1)+1,nx2*(ny2-1)+1,
     &                   nj,jrho,jrho*nx2,nx2,chrt,j1,ny1,sj1,ex)
                 else
                    call dvi4q1d(u2(i,1,k),u2(i,1,k),1+(j2-1)*nx2,
     &                   1+(j2-1)*nx2,nx2*(ny2-1)+1,nx2*(ny2-1)+1,
     &                   nj,jrho,jrho*nx2,nx2)
                 end if
              end if
           end do
        end do

        if (nk.eq.1) return

        !--------------------------------------------------------------
        ! now all k line segments, using points filled in above
        !--------------------------------------------------------------
        do i=i2,i2+(ni-1)*irho
           do j=j2,j2+(nj-1)*jrho
              if (do_ex.eq.1) then
                 call init_chrt(chrt,
     &              chr1(i1+si1*((i-i2)/irho),j1+sj1*((j-j2)/jrho),1),
     &              nz1,nx1*ny1*(nz1-1)+1,nx1*ny1,ex)
                 if (mod(i-i2,irho).ne.0) 
     &              call and_chrt(chrt,
     &              chr1(i1+si1*((i-i2)/irho+1),j1+sj1*((j-j2)/jrho),1),
     &              nz1,nx1*ny1*(nz1-1)+1,nx1*ny1,ex)
                 if (mod(j-j2,jrho).ne.0) 
     &              call and_chrt(chrt,
     &              chr1(i1+si1*((i-i2)/irho),j1+sj1*((j-j2)/jrho+1),1),
     &              nz1,nx1*ny1*(nz1-1)+1,nx1*ny1,ex)
                 if (mod(j-j2,jrho).ne.0.and.mod(i-i2,irho).ne.0) 
     &              call and_chrt(chrt,
     &            chr1(i1+si1*((i-i2)/irho+1),j1+sj1*((j-j2)/jrho+1),1),
     &              nz1,nx1*ny1*(nz1-1)+1,nx1*ny1,ex)
              end if
              if (ord.eq.2) then
                 if (do_ex.eq.1) then
                    call dvi2q1d_ex(u2(i,j,1),u2(i,j,1),
     &                 (k2-1)*nx2*ny2+1,(k2-1)*nx2*ny2+1,
     &                 nx2*ny2*(nz2-1)+1,nx2*ny2*(nz2-1)+1,nk,
     &                 krho,krho*nx2*ny2,nx2*ny2,chrt,k1,nz1,sk1,ex)
                 else
                    call dvi2q1d(u2(i,j,1),u2(i,j,1),
     &                 (k2-1)*nx2*ny2+1,(k2-1)*nx2*ny2+1,
     &                 nx2*ny2*(nz2-1)+1,nx2*ny2*(nz2-1)+1,nk,
     &                 krho,krho*nx2*ny2,nx2*ny2)
                 end if
              else
                 if (do_ex.eq.1) then
                    call dvi4q1d_ex(u2(i,j,1),u2(i,j,1),
     &                 (k2-1)*nx2*ny2+1,(k2-1)*nx2*ny2+1,
     &                 nx2*ny2*(nz2-1)+1,nx2*ny2*(nz2-1)+1,nk,
     &                 krho,krho*nx2*ny2,nx2*ny2,chrt,k1,nz1,sk1,ex)
                 else
                    call dvi4q1d(u2(i,j,1),u2(i,j,1),
     &                 (k2-1)*nx2*ny2+1,(k2-1)*nx2*ny2+1,
     &                 nx2*ny2*(nz2-1)+1,nx2*ny2*(nz2-1)+1,nk,
     &                 krho,krho*nx2*ny2,nx2*ny2)
                 end if
              end if
           end do
        end do

        return
        end

c---------------------------------------------------------------------------
c half-weight restriction (if s1=[2,2,2]): 
c dmhwr3d copies a segment of grid u1, with half-weight averaging, into grid u2.
c
c offset in u1 = [i1,j1,k1] with stride [si1,sj1,sk1],
c offset in u2 = [i2,j2,k2] with stride [si2,sj2,sk2]
c size of segment=ni x nj x nk in u2
c
c if (do_ex==1), then excise
c---------------------------------------------------------------------------
        subroutine dmhwr3d(u1,u2,nx1,ny1,nz1,nx2,ny2,nz2,
     &                     i1,j1,k1,i2,j2,k2,ni,nj,nk,
     &                     si1,sj1,sk1,si2,sj2,sk2,
     &                     chr1,ex,do_ex)
        implicit none
        integer nx1,ny1,nz1,nx2,ny2,nz2
        integer i1,j1,k1,i2,j2,k2,ni,nj,nk,do_ex
        integer si1,sj1,sk1,si2,sj2,sk2
        real*8 u1(nx1,ny1,nz1),u2(nx2,ny2,nz2)
        real*8 chr1(nx1,ny1,nz1),ex

        integer i,j,k,w,i0,j0,k0
        real *8 sum

        if (i1.lt.1.or.(i1+(ni-1)*si1).gt.nx1) return
        if (i2.lt.1.or.(i2+(ni-1)*si2).gt.nx2) return
         
        if (j1.lt.1.or.(j1+(nj-1)*sj1).gt.ny1) return
        if (j2.lt.1.or.(j2+(nj-1)*sj2).gt.ny2) return

        if (k1.lt.1.or.(k1+(nk-1)*sk1).gt.nz1) return
        if (k2.lt.1.or.(k2+(nk-1)*sk2).gt.nz2) return
        
        do i=1,ni
           do j=1,nj
              do k=1,nk
                 i0=i1+(i-1)*si1
                 j0=j1+(j-1)*sj1
                 k0=k1+(k-1)*sk1
                 w=0
                 sum=0
                 if (i0.gt.1.and.i0.le.(nx1-1)) then
                    if (do_ex.eq.0.or.
     &                  (chr1(i0-1,j0,k0).ne.ex.and.
     &                   chr1(i0,j0,k0).ne.ex.and.
     &                   chr1(i0+1,j0,k0).ne.ex)) then
                       w=w+4
                       sum=sum+
     &                    u1(i0-1,j0,k0)+2*u1(i0,j0,k0)+u1(i0+1,j0,k0)
                    end if
                 end if
                 if (j0.gt.1.and.j0.le.(ny1-1)) then
                    if (do_ex.eq.0.or.
     &                  (chr1(i0,j0-1,k0).ne.ex.and.
     &                   chr1(i0,j0,k0).ne.ex.and.
     &                   chr1(i0,j0+1,k0).ne.ex)) then
                       w=w+4
                       sum=sum+
     &                    u1(i0,j0-1,k0)+2*u1(i0,j0,k0)+u1(i0,j0+1,k0)
                    end if
                 end if
                 if (k0.gt.1.and.k0.le.(nz1-1)) then
                    if (do_ex.eq.0.or.
     &                  (chr1(i0,j0,k0-1).ne.ex.and.
     &                   chr1(i0,j0,k0).ne.ex.and.
     &                   chr1(i0,j0,k0+1).ne.ex)) then
                       w=w+4
                       sum=sum+
     &                    u1(i0,j0,k0-1)+2*u1(i0,j0,k0)+u1(i0,j0,k0+1)
                    end if
                 end if
                 if (do_ex.eq.1.and.chr1(i0,j0,k0).eq.ex) then
                    u2(i2+(i-1)*si2,j2+(j-1)*sj2,k2+(k-1)*sk2)=0
                 else if (w.gt.0) then
                    u2(i2+(i-1)*si2,j2+(j-1)*sj2,k2+(k-1)*sk2)=sum/w
                 else
                    u2(i2+(i-1)*si2,j2+(j-1)*sj2,k2+(k-1)*sk2)=
     &              u1(i0,j0,k0)
                 end if
              end do
           end do
        end do

        return
        end

c---------------------------------------------------------------------------
c full-weight restriction (is s1=[2,2,2]): 
c dmhwr3d copies a segment of grid u1, with full-weight averaging, into grid u2.
c
c offset in u1 = [i1,j1,k1] with stride [si1,sj1,sk1],
c offset in u2 = [i2,j2,k2] with stride [si2,sj2,sk2]
c size of segment=ni x nj x nk in u2
c---------------------------------------------------------------------------
        subroutine dmfwr3d(u1,u2,nx1,ny1,nz1,nx2,ny2,nz2,
     &                     i1,j1,k1,i2,j2,k2,ni,nj,nk,
     &                     si1,sj1,sk1,si2,sj2,sk2,
     &                     chr1,ex,do_ex)
        implicit none
        integer nx1,ny1,nz1,nx2,ny2,nz2,do_ex
        integer i1,j1,k1,i2,j2,k2,ni,nj,nk
        integer si1,sj1,sk1,si2,sj2,sk2
        real*8 u1(nx1,ny1,nz1),u2(nx2,ny2,nz2)
        real*8 chr1(nx1,ny1,nz1),ex

        integer i,j,k,i0,j0,k0,im,jm,km,ip,jp,kp

        if (i1.lt.1.or.(i1+(ni-1)*si1).gt.nx1) return
        if (i2.lt.1.or.(i2+(ni-1)*si2).gt.nx2) return
         
        if (j1.lt.1.or.(j1+(nj-1)*sj1).gt.ny1) return
        if (j2.lt.1.or.(j2+(nj-1)*sj2).gt.ny2) return

        if (k1.lt.1.or.(k1+(nk-1)*sk1).gt.nz1) return
        if (k2.lt.1.or.(k2+(nk-1)*sk2).gt.nz2) return

        if (do_ex.eq.1) then
           call dmfwr3d_ex(u1,u2,nx1,ny1,nz1,nx2,ny2,nz2,
     &                     i1,j1,k1,i2,j2,k2,ni,nj,nk,
     &                     si1,sj1,sk1,si2,sj2,sk2,
     &                     chr1,ex)
           return
        end if
        
        do i=1,ni
           do j=1,nj
              do k=1,nk
                 i0=i1+(i-1)*si1
                 j0=j1+(j-1)*sj1
                 k0=k1+(k-1)*sk1
                 im=i0-1
                 jm=j0-1
                 km=k0-1
                 ip=i0+1
                 jp=j0+1
                 kp=k0+1
                 if (i0.gt.1.and.i0.le.(nx1-1).and.
     &               j0.gt.1.and.j0.le.(ny1-1).and.
     &               k0.gt.1.and.k0.le.(nz1-1)) then
                  u2(i2+(i-1)*si2,j2+(j-1)*sj2,k2+(k-1)*sk2)=(
     &                u1(im,jm,km)+2*u1(i0,jm,km)+  u1(ip,jm,km)+
     &              2*u1(im,j0,km)+4*u1(i0,j0,km)+2*u1(ip,j0,km)+
     &                u1(im,jp,km)+2*u1(i0,jp,km)+  u1(ip,jp,km)+
     &
     &              2*u1(im,jm,k0)+4*u1(i0,jm,k0)+2*u1(ip,jm,k0)+
     &              4*u1(im,j0,k0)+8*u1(i0,j0,k0)+4*u1(ip,j0,k0)+
     &              2*u1(im,jp,k0)+4*u1(i0,jp,k0)+2*u1(ip,jp,k0)+
     &
     &                u1(im,jm,kp)+2*u1(i0,jm,kp)+  u1(ip,jm,kp)+
     &              2*u1(im,j0,kp)+4*u1(i0,j0,kp)+2*u1(ip,j0,kp)+
     &                u1(im,jp,kp)+2*u1(i0,jp,kp)+  u1(ip,jp,kp)
     &               )/64.0d0
                 else if (i0.gt.1.and.i0.le.(nx1-1).and.
     &                    j0.gt.1.and.j0.le.(ny1-1)) then
                  u2(i2+(i-1)*si2,j2+(j-1)*sj2,k2+(k-1)*sk2)=(
     &              1*u1(im,jm,k0)+2*u1(i0,jm,k0)+1*u1(ip,jm,k0)+
     &              2*u1(im,j0,k0)+4*u1(i0,j0,k0)+2*u1(ip,j0,k0)+
     &              1*u1(im,jp,k0)+2*u1(i0,jp,k0)+1*u1(ip,jp,k0)
     &               )/16.0d0
                 else if (j0.gt.1.and.j0.le.(ny1-1).and.
     &                    k0.gt.1.and.k0.le.(nz1-1)) then
                  u2(i2+(i-1)*si2,j2+(j-1)*sj2,k2+(k-1)*sk2)=(
     &              1*u1(i0,jm,km)+2*u1(i0,jm,k0)+1*u1(i0,jm,kp)+
     &              2*u1(i0,j0,km)+4*u1(i0,j0,k0)+2*u1(i0,j0,kp)+
     &              1*u1(i0,jp,km)+2*u1(i0,jp,k0)+1*u1(i0,jp,kp)
     &               )/16.0d0
                 else if (i0.gt.1.and.i0.le.(nx1-1).and.
     &                    k0.gt.1.and.k0.le.(nz1-1)) then
                  u2(i2+(i-1)*si2,j2+(j-1)*sj2,k2+(k-1)*sk2)=(
     &              1*u1(im,j0,km)+2*u1(im,j0,k0)+1*u1(im,j0,kp)+
     &              2*u1(i0,j0,km)+4*u1(i0,j0,k0)+2*u1(i0,j0,kp)+
     &              1*u1(ip,j0,km)+2*u1(ip,j0,k0)+1*u1(ip,j0,kp)
     &               )/16.0d0
                 else if (i0.gt.1.and.i0.le.(nx1-1)) then
                  u2(i2+(i-1)*si2,j2+(j-1)*sj2,k2+(k-1)*sk2)=
     &              (u1(im,j0,k0)+2*u1(i0,j0,k0)+u1(ip,j0,k0))/4.0d0
                 else if (j0.gt.1.and.j0.le.(ny1-1)) then
                  u2(i2+(i-1)*si2,j2+(j-1)*sj2,k2+(k-1)*sk2)=
     &              (u1(i0,jm,k0)+2*u1(i0,j0,k0)+u1(i0,jp,k0))/4.0d0
                 else if (k0.gt.1.and.k0.le.(nz1-1)) then
                  u2(i2+(i-1)*si2,j2+(j-1)*sj2,k2+(k-1)*sk2)=
     &              (u1(i0,j0,km)+2*u1(i0,j0,k0)+u1(i0,j0,kp))/4.0d0
                 else
                  u2(i2+(i-1)*si2,j2+(j-1)*sj2,k2+(k-1)*sk2)=
     &               u1(i0,j0,k0)
                 end if
              end do
           end do
        end do

        return
        end

c---------------------------------------------------------------------------
c called by dmfwr3d for excision ... bounds checking already done
c---------------------------------------------------------------------------
        subroutine dmfwr3d_ex(u1,u2,nx1,ny1,nz1,nx2,ny2,nz2,
     &                     i1,j1,k1,i2,j2,k2,ni,nj,nk,
     &                     si1,sj1,sk1,si2,sj2,sk2,
     &                     chr1,ex)
        implicit none
        integer nx1,ny1,nz1,nx2,ny2,nz2
        integer i1,j1,k1,i2,j2,k2,ni,nj,nk
        integer si1,sj1,sk1,si2,sj2,sk2
        real*8 u1(nx1,ny1,nz1),u2(nx2,ny2,nz2)
        real*8 chr1(nx1,ny1,nz1),ex

        integer i,j,k,i0,j0,k0,im,jm,km,ip,jp,kp
        integer i01,j01,k01

        integer w3d(3,3,3),w
        real*8 sum
        data w3d/1,2,1,2,4,2,1,2,1,
     &           2,4,2,4,8,4,2,4,2,
     &           1,2,1,2,4,2,1,2,1/

        do i=1,ni
           do j=1,nj
              do k=1,nk
                 i0=i1+(i-1)*si1
                 j0=j1+(j-1)*sj1
                 k0=k1+(k-1)*sk1
                 im=i0-1
                 jm=j0-1
                 km=k0-1
                 ip=i0+1
                 jp=j0+1
                 kp=k0+1
                 sum=0
                 w=0
                 if (i0.gt.1.and.i0.le.(nx1-1).and.
     &               j0.gt.1.and.j0.le.(ny1-1).and.
     &               k0.gt.1.and.k0.le.(nz1-1)) then
                    do i01=im,ip
                       do j01=jm,jp
                          do k01=km,kp
                             if (chr1(i01,j01,k01).ne.ex) then  
                                sum=sum+w3d(i01-im+1,j01-jm+1,k01-km+1)*
     &                                  u1(i01,j01,k01)
                                w=w+w3d(i01-im+1,j01-jm+1,k01-km+1)
                             end if
                          end do
                       end do
                    end do
                 else if (i0.gt.1.and.i0.le.(nx1-1).and.
     &                    j0.gt.1.and.j0.le.(ny1-1)) then
                    do i01=im,ip
                       do j01=jm,jp
                        if (chr1(i01,j01,k0).ne.ex) then  
                         sum=sum+w3d(i01-im+1,j01-jm+1,1)*u1(i01,j01,k0)
                         w=w+w3d(i01-im+1,j01-jm+1,1)
                        end if
                       end do
                    end do
                 else if (j0.gt.1.and.j0.le.(ny1-1).and.
     &                    k0.gt.1.and.k0.le.(nz1-1)) then
                    do j01=jm,jp
                       do k01=km,kp
                        if (chr1(i0,j01,k01).ne.ex) then  
                         sum=sum+w3d(1,j01-jm+1,k01-km+1)*u1(i0,j01,k01)
                         w=w+w3d(1,j01-jm+1,k01-km+1)
                        end if
                       end do
                    end do
                 else if (i0.gt.1.and.i0.le.(nx1-1).and.
     &                    k0.gt.1.and.k0.le.(nz1-1)) then
                    do i01=im,ip
                       do k01=km,kp
                        if (chr1(i01,j0,k01).ne.ex) then  
                         sum=sum+w3d(i01-im+1,1,k01-km+1)*u1(i01,j0,k01)
                         w=w+w3d(i01-im+1,1,k01-km+1)
                        end if
                       end do
                    end do
                 else if (i0.gt.1.and.i0.le.(nx1-1)) then
                    do i01=im,ip
                       if (chr1(i01,j0,k0).ne.ex) then  
                          sum=sum+w3d(i01-im+1,1,1)*u1(i01,j0,k0)
                          w=w+w3d(i01-im+1,1,1)
                       end if
                    end do
                 else if (j0.gt.1.and.j0.le.(ny1-1)) then
                    do j01=jm,jp
                       if (chr1(i0,j01,k0).ne.ex) then  
                          sum=sum+w3d(1,j01-jm+1,1)*u1(i0,j01,k0)
                          w=w+w3d(1,j01-jm+1,1)
                       end if
                    end do
                 else if (k0.gt.1.and.k0.le.(nz1-1)) then
                    do k01=km,kp
                       if (chr1(i0,j0,k01).ne.ex) then  
                          sum=sum+w3d(1,1,k01-km+1)*u1(i0,j0,k01)
                          w=w+w3d(i,1,k01-km+1)
                       end if
                    end do
                 else
                    sum=u1(i0,j0,k0)
                    w=1
                 end if

                 if (w.eq.0.or.chr1(i0,j0,k0).eq.ex) then
                    u2(i2+(i-1)*si2,j2+(j-1)*sj2,k2+(k-1)*sk2)=0
                 else 
                    u2(i2+(i-1)*si2,j2+(j-1)*sj2,k2+(k-1)*sk2)=sum/w
                 end if
              end do
           end do
        end do

        return
        end

c---------------------------------------------------------------------------
c applies a simple KO filter to f, where mask != mask_off, and modified
c near excision boundaries as follows:
c (or all boundaries if do_bdy=1, 
c  or neither if do_bdy=0 and internal flag no_uc_ex_bdy=.true.)
c
c with undivided operators Up(f) = f(i+1)-f(i)
c                          Um(f) = f(i)-f(i-1)
c
c Interior: eps*w4*(Up Um)^2
c left+1  : eps*w3*(Up Up Um)
c left    : eps*w2*(Up Um)
c right-1 : eps*w3*(Um Um Up)
c right   : eps*w2*(Up Um)
c
c NOTE: update DV version too
c
c
c NOTE: SPECIAL FLAG ... if do_ex>0, eps is treated as a scalar, 
c       else eps is treated as an array.
c------------------------------------------------------------------------
        subroutine dmdiss3d_ex(f,work,eps,do_bdy,phys_bdy_type,even,
     &                         odd,mask,mask_off,nx,ny,
     &                         nz,chr,ex,do_ex,ind_sweeps)
        implicit none
        integer nx,ny,nz,do_bdy,phys_bdy_type(6),even,odd,do_ex
        real*8 f(nx,ny,nz),work(nx,ny,nz),mask(nx,ny,nz),eps(nx,ny,nz)
        real*8 mask_off
        real*8 chr(nx,ny,nz),ex
        logical ind_sweeps

        integer i,j,k,bo1,bo2
        real*8 eps_eff,f_hf,norm_f
        real*8 w4,w3,w2
c        parameter (w4=1.0d0/16.0d0,w3=1.0d0/16.0d0,w2=1.0d0/16.0d0)
        parameter (w4=1.0d0/16.0d0,w3=1.0d0/8.0d0,w2=1.0d0/4.0d0)
        integer pass,npass
        logical no_uc_ex_diss
        parameter (no_uc_ex_diss=.true.)
        
        eps_eff=eps(1,1,1)

        if (do_bdy.eq.0) then
           bo1=0
           bo2=0
           if (no_uc_ex_diss) then
              bo1=-max(nx,ny,nz)
              bo2=bo1
           end if
        else
           bo1=1
           bo2=2
        end if

        npass=1
        if (ny.gt.1.and.ind_sweeps) npass=2
        if (nz.gt.1.and.ind_sweeps) npass=3

        do pass=1,npass

         do i=1,nx
           do j=1,ny
              do k=1,nz
                 work(i,j,k)=f(i,j,k)
              end do
           end do
         end do

         do i=1,nx
           do j=1,ny
              do k=1,nz
                 if (mask(i,j,k).ne.mask_off.and.
     &               (chr(i,j,k).ne.ex)) then

                   if (do_ex.lt.0) eps_eff=eps(i,j,k)

!LR: equivalent to the commented out piece below, but it does not go out of bounds
                   if (.not.ind_sweeps.or.pass.eq.1) then
                    f_hf=0
                    if (i.gt.2.and.i.lt.(nx-1)) then
                     if ((chr(i-2,j,k).ne.ex.and.chr(i-1,j,k).ne.ex.and.
     &                    chr(i+2,j,k).ne.ex.and.chr(i+1,j,k).ne.ex))
     &               then
                       f_hf=w4*(
     &                     work(i-2,j,k)+work(i+2,j,k)
     &                 -4*(work(i-1,j,k)+work(i+1,j,k))+6*work(i,j,k))
                     end if
                    else if (i.eq.2.and.phys_bdy_type(1).eq.odd) then
                     if ((chr(i-1,j,k).ne.ex.and.
     &                    chr(i+2,j,k).ne.ex.and.chr(i+1,j,k).ne.ex))
     &               then
                       f_hf=w4*(
     &                    (-work(i,j,k))+work(i+2,j,k)
     &                 -4*(work(i-1,j,k)+work(i+1,j,k))+6*work(i,j,k))
                     end if
                    else if (i.eq.1.and.phys_bdy_type(1).eq.odd) then
                     if ((chr(i+2,j,k).ne.ex.and.chr(i+1,j,k).ne.ex))
     &               then
                       f_hf=w4*(
     &                    (-work(i+2,j,k))+work(i+2,j,k)
     &                -4*((-work(i+1,j,k))+work(i+1,j,k))+6*work(i,j,k))
                     end if
                    else if (i.eq.Nx-1.and.phys_bdy_type(2).eq.odd) then
                     if ((chr(i+1,j,k).ne.ex.and.
     &                    chr(i-2,j,k).ne.ex.and.chr(i-1,j,k).ne.ex))
     &               then
                       f_hf=w4*(
     &                     work(i-2,j,k)+(-work(i,j,k))
     &                 -4*(work(i-1,j,k)+work(i+1,j,k))+6*work(i,j,k))
                     end if
                    else if (i.eq.Nx.and.phys_bdy_type(2).eq.odd) then
                     if ((chr(i-2,j,k).ne.ex.and.chr(i-1,j,k).ne.ex))
     &               then
                       f_hf=w4*(
     &                    work(i-2,j,k)+(-work(i-2,j,k))
     &                -4*(work(i-1,j,k)+(-work(i-1,j,k)))+6*work(i,j,k))
                     end if
                    else if (i.eq.2.and.
     &               phys_bdy_type(1).eq.even) then
                     if ((chr(i-1,j,k).ne.ex.and.
     &                    chr(i+2,j,k).ne.ex.and.chr(i+1,j,k).ne.ex))
     &               then
                       f_hf=w4*(
     &                    (work(i,j,k))+work(i+2,j,k)
     &                 -4*(work(i-1,j,k)+work(i+1,j,k))+6*work(i,j,k))
                     end if
                    else if (i.eq.1.and.phys_bdy_type(1).eq.even) then
                     if ((chr(i+2,j,k).ne.ex.and.chr(i+1,j,k).ne.ex))
     &               then
                       f_hf=w4*(
     &                    (work(i+2,j,k))+work(i+2,j,k)
     &                -4*((work(i+1,j,k))+work(i+1,j,k))+6*work(i,j,k))
                     end if
                    else if (i.eq.Nx-1.and.
     &               phys_bdy_type(2).eq.even) then
                     if ((chr(i+1,j,k).ne.ex.and.
     &                    chr(i-2,j,k).ne.ex.and.chr(i-1,j,k).ne.ex))
     &               then
                       f_hf=w4*(
     &                     work(i-2,j,k)+(work(i,j,k))
     &                 -4*(work(i-1,j,k)+work(i+1,j,k))+6*work(i,j,k))
                     end if
                    else if (i.eq.Nx.and.phys_bdy_type(2).eq.even) then
                     if ((chr(i-2,j,k).ne.ex.and.chr(i-1,j,k).ne.ex))
     &               then
                       f_hf=w4*(
     &                    work(i-2,j,k)+(work(i-2,j,k))
     &                -4*(work(i-1,j,k)+(work(i-1,j,k)))+6*work(i,j,k))
                     end if
                    else if (i.gt.(2-bo1).and.i.lt.(nx-1)) then
                     if (chr(i-1,j,k).ne.ex.and.
     &                   chr(i+1,j,k).ne.ex.and.
     &                   chr(i+2,j,k).ne.ex)
     &               then
                       f_hf=w3*(
     &                    -work(i-1,j,k)+3*work(i,j,k)
     &                  -3*work(i+1,j,k)+work(i+2,j,k))
                     end if
                    else if (i.gt.2.and.i.lt.(nx-1+bo1)) then
                     if (chr(i-1,j,k).ne.ex.and.chr(i+1,j,k).ne.ex.and.
     &                   chr(i-2,j,k).ne.ex)
     &               then
                       f_hf=w3*(
     &                    -work(i+1,j,k)+3*work(i,j,k)
     &                  -3*work(i-1,j,k)+work(i-2,j,k))
                     end if
                    else if ((i.gt.(2-bo2).or.
     &                       (i.eq.2.and.chr(1,j,k).eq.ex))
     &                       .and.i.lt.(nx-1).and.
     &                  (chr(i+1,j,k).ne.ex.and.chr(i+2,j,k).ne.ex) )
     &              then
                       f_hf=w2*(
     &                     work(i,j,k)-2*work(i+1,j,k)+work(i+2,j,k))
                    else if (i.gt.2.and.(i.lt.(nx-1+bo2).or.
     &                       (i.eq.(nx-1).and.chr(nx,j,k).eq.ex)).and.
     &                  (chr(i-1,j,k).ne.ex.and.chr(i-2,j,k).ne.ex) )
     &              then
                       f_hf=w2*(
     &                     work(i,j,k)-2*work(i-1,j,k)+work(i-2,j,k))
                    end if

                    f(i,j,k)=f(i,j,k)-eps_eff*f_hf
                   end if

!LR commented this out
!                   if (.not.ind_sweeps.or.pass.eq.1) then
!                    f_hf=0
!                    if (i.gt.2.and.i.lt.(nx-1).and.
!     &                  ((chr(i-2,j,k).ne.ex.and.chr(i-1,j,k).ne.ex.and.
!     &                    chr(i+2,j,k).ne.ex.and.chr(i+1,j,k).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                     work(i-2,j,k)+work(i+2,j,k)
!     &                 -4*(work(i-1,j,k)+work(i+1,j,k))+6*work(i,j,k))
!                    else if (i.eq.2.and.phys_bdy_type(1).eq.odd.and.
!     &                  ((chr(i-1,j,k).ne.ex.and.
!     &                    chr(i+2,j,k).ne.ex.and.chr(i+1,j,k).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                    (-work(i,j,k))+work(i+2,j,k)
!     &                 -4*(work(i-1,j,k)+work(i+1,j,k))+6*work(i,j,k))
!                    else if (i.eq.1.and.phys_bdy_type(1).eq.odd.and.
!     &                 ((chr(i+2,j,k).ne.ex.and.chr(i+1,j,k).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                    (-work(i+2,j,k))+work(i+2,j,k)
!     &                -4*((-work(i+1,j,k))+work(i+1,j,k))+6*work(i,j,k))
!                    else if (i.eq.Nx-1.and.phys_bdy_type(2).eq.odd.and.
!     &                  ((chr(i+1,j,k).ne.ex.and.
!     &                    chr(i-2,j,k).ne.ex.and.chr(i-1,j,k).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                     work(i-2,j,k)+(-work(i,j,k))
!     &                 -4*(work(i-1,j,k)+work(i+1,j,k))+6*work(i,j,k))
!                    else if (i.eq.Nx.and.phys_bdy_type(2).eq.odd.and.
!     &                 ((chr(i-2,j,k).ne.ex.and.chr(i-1,j,k).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                    work(i-2,j,k)+(-work(i-2,j,k))
!     &                -4*(work(i-1,j,k)+(-work(i-1,j,k)))+6*work(i,j,k))
!                    else if (i.eq.2.and.phys_bdy_type(1).eq.even.and.
!     &                  ((chr(i-1,j,k).ne.ex.and.
!     &                    chr(i+2,j,k).ne.ex.and.chr(i+1,j,k).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                    (work(i,j,k))+work(i+2,j,k)
!     &                 -4*(work(i-1,j,k)+work(i+1,j,k))+6*work(i,j,k))
!                    else if (i.eq.1.and.phys_bdy_type(1).eq.even.and.
!     &                 ((chr(i+2,j,k).ne.ex.and.chr(i+1,j,k).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                    (work(i+2,j,k))+work(i+2,j,k)
!     &                -4*((work(i+1,j,k))+work(i+1,j,k))+6*work(i,j,k))
!                    else if (i.eq.Nx-1.and.phys_bdy_type(2).eq.even.and.
!     &                  ((chr(i+1,j,k).ne.ex.and.
!     &                    chr(i-2,j,k).ne.ex.and.chr(i-1,j,k).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                     work(i-2,j,k)+(work(i,j,k))
!     &                 -4*(work(i-1,j,k)+work(i+1,j,k))+6*work(i,j,k))
!                    else if (i.eq.Nx.and.phys_bdy_type(2).eq.even.and.
!     &                 ((chr(i-2,j,k).ne.ex.and.chr(i-1,j,k).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                    work(i-2,j,k)+(work(i-2,j,k))
!     &                -4*(work(i-1,j,k)+(work(i-1,j,k)))+6*work(i,j,k))
!                    else if (i.gt.(2-bo1).and.i.lt.(nx-1).and.
!     &                  (chr(i-1,j,k).ne.ex.and.chr(i+1,j,k).ne.ex.and.
!     &                   chr(i+2,j,k).ne.ex)) 
!     &              then
!                       f_hf=w3*(
!     &                    -work(i-1,j,k)+3*work(i,j,k)
!     &                  -3*work(i+1,j,k)+work(i+2,j,k))
!                    else if (i.gt.2.and.i.lt.(nx-1+bo1).and.
!     &                  (chr(i-1,j,k).ne.ex.and.chr(i+1,j,k).ne.ex.and.
!     &                   chr(i-2,j,k).ne.ex) ) 
!     &              then
!                       f_hf=w3*(
!     &                    -work(i+1,j,k)+3*work(i,j,k)
!     &                  -3*work(i-1,j,k)+work(i-2,j,k))
!                    else if ((i.gt.(2-bo2).or.
!     &                       (i.eq.2.and.chr(1,j,k).eq.ex))
!     &                       .and.i.lt.(nx-1).and.
!     &                  (chr(i+1,j,k).ne.ex.and.chr(i+2,j,k).ne.ex) )
!     &              then
!                       f_hf=w2*(
!     &                     work(i,j,k)-2*work(i+1,j,k)+work(i+2,j,k))
!                    else if (i.gt.2.and.(i.lt.(nx-1+bo2).or.
!     &                       (i.eq.(nx-1).and.chr(nx,j,k).eq.ex)).and.
!     &                  (chr(i-1,j,k).ne.ex.and.chr(i-2,j,k).ne.ex) )
!     &              then
!                       f_hf=w2*(
!     &                     work(i,j,k)-2*work(i-1,j,k)+work(i-2,j,k))
!                    end if
!
!                    f(i,j,k)=f(i,j,k)-eps_eff*f_hf
!                   end if

!LR: equivalent to the commented out piece below, but it does not go out of bounds
                   if (.not.ind_sweeps.or.pass.eq.2) then
                    f_hf=0
                    if (j.gt.2.and.j.lt.(ny-1)) then
                     if ((chr(i,j-2,k).ne.ex.and.chr(i,j-1,k).ne.ex.and.
     &                    chr(i,j+2,k).ne.ex.and.chr(i,j+1,k).ne.ex))
     &               then
                       f_hf=w4*(
     &                     work(i,j-2,k)+work(i,j+2,k)
     &                 -4*(work(i,j-1,k)+work(i,j+1,k))+6*work(i,j,k))
                     end if
                    else if (j.eq.2.and.phys_bdy_type(3).eq.odd) then
                     if ((chr(i,j-1,k).ne.ex.and.
     &                    chr(i,j+2,k).ne.ex.and.chr(i,j+1,k).ne.ex))
     &               then
                       f_hf=w4*(
     &                    (-work(i,j,k))+work(i,j+2,k)
     &                 -4*(work(i,j-1,k)+work(i,j+1,k))+6*work(i,j,k))
                     end if
                    else if (j.eq.1.and.phys_bdy_type(3).eq.odd) then
                     if ((chr(i,j+2,k).ne.ex.and.chr(i,j+1,k).ne.ex))
     &               then
                       f_hf=w4*(
     &                    (-work(i,j+2,k))+work(i,j+2,k)
     &                -4*((-work(i,j+1,k))+work(i,j+1,k))+6*work(i,j,k))
                     end if
                    else if (j.eq.Ny-1.and.phys_bdy_type(4).eq.odd) then
                     if ((chr(i,j+1,k).ne.ex.and.
     &                    chr(i,j-2,k).ne.ex.and.chr(i,j-1,k).ne.ex))
     &               then
                       f_hf=w4*(
     &                     work(i,j-2,k)+(-work(i,j,k))
     &                 -4*(work(i,j-1,k)+work(i,j+1,k))+6*work(i,j,k))
                     end if
                    else if (j.eq.Ny.and.phys_bdy_type(4).eq.odd) then
                     if ((chr(i,j-2,k).ne.ex.and.chr(i,j-1,k).ne.ex))
     &               then
                       f_hf=w4*(
     &                    work(i,j-2,k)+(-work(i,j-2,k))
     &                -4*(work(i,j-1,k)+(-work(i,j-1,k)))+6*work(i,j,k))
                     end if
                    else if (j.eq.2.and.phys_bdy_type(3).eq.even) then
                     if ((chr(i,j-1,k).ne.ex.and.
     &                    chr(i,j+2,k).ne.ex.and.chr(i,j+1,k).ne.ex))
     &               then
                       f_hf=w4*(
     &                    (work(i,j,k))+work(i,j+2,k)
     &                 -4*(work(i,j-1,k)+work(i,j+1,k))+6*work(i,j,k))
                     end if
                    else if (j.eq.1.and.phys_bdy_type(3).eq.even) then
                     if ((chr(i,j+2,k).ne.ex.and.chr(i,j+1,k).ne.ex))
     &               then
                       f_hf=w4*(
     &                    (work(i,j+2,k))+work(i,j+2,k)
     &                -4*((work(i,j+1,k))+work(i,j+1,k))+6*work(i,j,k))
                     end if
                    else if (j.eq.Ny-1.and.
     &                       phys_bdy_type(4).eq.even) then
                     if ((chr(i,j+1,k).ne.ex.and.
     &                    chr(i,j-2,k).ne.ex.and.chr(i,j-1,k).ne.ex))
     &               then
                       f_hf=w4*(
     &                     work(i,j-2,k)+(work(i,j,k))
     &                 -4*(work(i,j-1,k)+work(i,j+1,k))+6*work(i,j,k))
                     end if
                    else if (j.eq.Ny.and.phys_bdy_type(4).eq.even) then
                     if ((chr(i,j-2,k).ne.ex.and.chr(i,j-1,k).ne.ex))
     &               then
                       f_hf=w4*(
     &                    work(i,j-2,k)+(work(i,j-2,k))
     &                -4*(work(i,j-1,k)+(work(i,j-1,k)))+6*work(i,j,k))
                     end if
                    else if (j.gt.(2-bo1).and.j.lt.(ny-1)) then
                     if (chr(i,j-1,k).ne.ex.and.chr(i,j+1,k).ne.ex.and.
     &                   chr(i,j+2,k).ne.ex)
     &               then
                       f_hf=w3*(
     &                    -work(i,j-1,k)+3*work(i,j,k)
     &                  -3*work(i,j+1,k)+work(i,j+2,k))
                     end if
                    else if (j.gt.2.and.j.lt.(ny-1+bo1)) then
                     if (chr(i,j-1,k).ne.ex.and.chr(i,j+1,k).ne.ex.and.
     &                   chr(i,j-2,k).ne.ex)
     &              then
                       f_hf=w3*(
     &                    -work(i,j+1,k)+3*work(i,j,k)
     &                  -3*work(i,j-1,k)+work(i,j-2,k))
                     end if
                    else if ((j.gt.(2-bo2).or.
     &                       (j.eq.2.and.chr(i,1,k).eq.ex))
     &                       .and.j.lt.(ny-1).and.
     &                  (chr(i,j+1,k).ne.ex.and.chr(i,j+2,k).ne.ex) )
     &              then
                       f_hf=w2*(
     &                     work(i,j,k)-2*work(i,j+1,k)+work(i,j+2,k))
                    else if (j.gt.2.and.(j.lt.(ny-1+bo2).or.
     &                       (j.eq.(ny-1).and.chr(i,ny,k).eq.ex)).and.
     &                  (chr(i,j-1,k).ne.ex.and.chr(i,j-2,k).ne.ex) )
     &              then
                       f_hf=w2*(
     &                     work(i,j,k)-2*work(i,j-1,k)+work(i,j-2,k))
                    end if

                    f(i,j,k)=f(i,j,k)-eps_eff*f_hf
                   end if
!LR commented this out
!                   if (.not.ind_sweeps.or.pass.eq.2) then
!                    f_hf=0
!                    if (j.gt.2.and.j.lt.(ny-1).and.
!     &                  ((chr(i,j-2,k).ne.ex.and.chr(i,j-1,k).ne.ex.and.
!     &                    chr(i,j+2,k).ne.ex.and.chr(i,j+1,k).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                     work(i,j-2,k)+work(i,j+2,k)
!     &                 -4*(work(i,j-1,k)+work(i,j+1,k))+6*work(i,j,k))
!                    else if (j.eq.2.and.phys_bdy_type(3).eq.odd.and.
!     &                  ((chr(i,j-1,k).ne.ex.and.
!     &                    chr(i,j+2,k).ne.ex.and.chr(i,j+1,k).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                    (-work(i,j,k))+work(i,j+2,k)
!     &                 -4*(work(i,j-1,k)+work(i,j+1,k))+6*work(i,j,k))
!                    else if (j.eq.1.and.phys_bdy_type(3).eq.odd.and.
!     &                  ((chr(i,j+2,k).ne.ex.and.chr(i,j+1,k).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                    (-work(i,j+2,k))+work(i,j+2,k)
!     &                -4*((-work(i,j+1,k))+work(i,j+1,k))+6*work(i,j,k))
!                    else if (j.eq.Ny-1.and.phys_bdy_type(4).eq.odd.and.
!     &                  ((chr(i,j+1,k).ne.ex.and.
!     &                    chr(i,j-2,k).ne.ex.and.chr(i,j-1,k).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                     work(i,j-2,k)+(-work(i,j,k))
!     &                 -4*(work(i,j-1,k)+work(i,j+1,k))+6*work(i,j,k))
!                    else if (j.eq.Ny.and.phys_bdy_type(4).eq.odd.and.
!     &                  ((chr(i,j-2,k).ne.ex.and.chr(i,j-1,k).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                    work(i,j-2,k)+(-work(i,j-2,k))
!     &                -4*(work(i,j-1,k)+(-work(i,j-1,k)))+6*work(i,j,k))
!                    else if (j.eq.2.and.phys_bdy_type(3).eq.even.and.
!     &                  ((chr(i,j-1,k).ne.ex.and.
!     &                    chr(i,j+2,k).ne.ex.and.chr(i,j+1,k).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                    (work(i,j,k))+work(i,j+2,k)
!     &                 -4*(work(i,j-1,k)+work(i,j+1,k))+6*work(i,j,k))
!                    else if (j.eq.1.and.phys_bdy_type(3).eq.even.and.
!     &                  ((chr(i,j+2,k).ne.ex.and.chr(i,j+1,k).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                    (work(i,j+2,k))+work(i,j+2,k)
!     &                -4*((work(i,j+1,k))+work(i,j+1,k))+6*work(i,j,k))
!                    else if (j.eq.Ny-1.and.phys_bdy_type(4).eq.even.and.
!     &                  ((chr(i,j+1,k).ne.ex.and.
!     &                    chr(i,j-2,k).ne.ex.and.chr(i,j-1,k).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                     work(i,j-2,k)+(work(i,j,k))
!     &                 -4*(work(i,j-1,k)+work(i,j+1,k))+6*work(i,j,k))
!                    else if (j.eq.Ny.and.phys_bdy_type(4).eq.even.and.
!     &                  ((chr(i,j-2,k).ne.ex.and.chr(i,j-1,k).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                    work(i,j-2,k)+(work(i,j-2,k))
!     &                -4*(work(i,j-1,k)+(work(i,j-1,k)))+6*work(i,j,k))
!                    else if (j.gt.(2-bo1).and.j.lt.(ny-1).and.
!     &                  (chr(i,j-1,k).ne.ex.and.chr(i,j+1,k).ne.ex.and.
!     &                   chr(i,j+2,k).ne.ex)) 
!     &              then
!                       f_hf=w3*(
!     &                    -work(i,j-1,k)+3*work(i,j,k)
!     &                  -3*work(i,j+1,k)+work(i,j+2,k))
!                    else if (j.gt.2.and.j.lt.(ny-1+bo1).and.
!     &                  (chr(i,j-1,k).ne.ex.and.chr(i,j+1,k).ne.ex.and.
!     &                   chr(i,j-2,k).ne.ex) ) 
!     &              then
!                       f_hf=w3*(
!     &                    -work(i,j+1,k)+3*work(i,j,k)
!     &                  -3*work(i,j-1,k)+work(i,j-2,k))
!                    else if ((j.gt.(2-bo2).or.
!     &                       (j.eq.2.and.chr(i,1,k).eq.ex))
!     &                       .and.j.lt.(ny-1).and.
!     &                  (chr(i,j+1,k).ne.ex.and.chr(i,j+2,k).ne.ex) )
!     &              then
!                       f_hf=w2*(
!     &                     work(i,j,k)-2*work(i,j+1,k)+work(i,j+2,k))
!                    else if (j.gt.2.and.(j.lt.(ny-1+bo2).or.
!     &                       (j.eq.(ny-1).and.chr(i,ny,k).eq.ex)).and.
!     &                  (chr(i,j-1,k).ne.ex.and.chr(i,j-2,k).ne.ex) )
!     &              then
!                       f_hf=w2*(
!     &                     work(i,j,k)-2*work(i,j-1,k)+work(i,j-2,k))
!                    end if
!
!                    f(i,j,k)=f(i,j,k)-eps_eff*f_hf
!                   end if

!LR: equivalent to the commented out piece below, but it does not go out of bounds
                   if (.not.ind_sweeps.or.pass.eq.3) then
                    f_hf=0
                    if (k.gt.2.and.k.lt.(nz-1)) then
                     if ((chr(i,j,k-2).ne.ex.and.chr(i,j,k-1).ne.ex.and.
     &                    chr(i,j,k+2).ne.ex.and.chr(i,j,k+1).ne.ex))
     &               then
                       f_hf=w4*(
     &                     work(i,j,k-2)+work(i,j,k+2)
     &                 -4*(work(i,j,k-1)+work(i,j,k+1))+6*work(i,j,k))
                     end if
                    else if (k.eq.2.and.phys_bdy_type(5).eq.odd) then
                     if ((chr(i,j,k-1).ne.ex.and.
     &                    chr(i,j,k+2).ne.ex.and.chr(i,j,k+1).ne.ex))
     &               then
                       f_hf=w4*(
     &                    (-work(i,j,k))+work(i,j,k+2)
     &                 -4*(work(i,j,k-1)+work(i,j,k+1))+6*work(i,j,k))
                     end if
                    else if (k.eq.1.and.phys_bdy_type(5).eq.odd) then
                     if ((chr(i,j,k+2).ne.ex.and.chr(i,j,k+1).ne.ex))
     &               then
                       f_hf=w4*(
     &                    (-work(i,j,k+2))+work(i,j,k+2)
     &                -4*((-work(i,j,k+1))+work(i,j,k+1))+6*work(i,j,k))
                     end if
                    else if (k.eq.Nz-1.and.phys_bdy_type(6).eq.odd) then
                     if ((chr(i,j,k-2).ne.ex.and.chr(i,j,k-1).ne.ex.and.
     &                    chr(i,j,k+1).ne.ex))
     &               then
                       f_hf=w4*(
     &                     work(i,j,k-2)+(-work(i,j,k))
     &                 -4*(work(i,j,k-1)+work(i,j,k+1))+6*work(i,j,k))
                     end if
                    else if (k.eq.Nz.and.phys_bdy_type(6).eq.odd) then
                     if ((chr(i,j,k-2).ne.ex.and.chr(i,j,k-1).ne.ex))
     &               then
                       f_hf=w4*(
     &                    work(i,j,k-2)+(-work(i,j,k-2))
     &                -4*(work(i,j,k-1)+(-work(i,j,k-1)))+6*work(i,j,k))
                     end if
                    else if (k.eq.2.and.phys_bdy_type(5).eq.even) then
                     if ((chr(i,j,k-1).ne.ex.and.
     &                    chr(i,j,k+2).ne.ex.and.chr(i,j,k+1).ne.ex))
     &               then
                       f_hf=w4*(
     &                    (work(i,j,k))+work(i,j,k+2)
     &                 -4*(work(i,j,k-1)+work(i,j,k+1))+6*work(i,j,k))
                     end if
                    else if (k.eq.1.and.phys_bdy_type(5).eq.even) then
                     if ((chr(i,j,k+2).ne.ex.and.chr(i,j,k+1).ne.ex))
     &               then
                       f_hf=w4*(
     &                    (work(i,j,k+2))+work(i,j,k+2)
     &                -4*((work(i,j,k+1))+work(i,j,k+1))+6*work(i,j,k))
                     end if
                    else if (k.eq.Nz-1.and.
     &                       phys_bdy_type(6).eq.even) then
                     if ((chr(i,j,k-2).ne.ex.and.chr(i,j,k-1).ne.ex.and.
     &                    chr(i,j,k+1).ne.ex))
     &               then
                       f_hf=w4*(
     &                     work(i,j,k-2)+(work(i,j,k))
     &                 -4*(work(i,j,k-1)+work(i,j,k+1))+6*work(i,j,k))
                     end if
                    else if (k.eq.Nz.and.phys_bdy_type(6).eq.even) then
                     if ((chr(i,j,k-2).ne.ex.and.chr(i,j,k-1).ne.ex))
     &               then
                       f_hf=w4*(
     &                    work(i,j,k-2)+(work(i,j,k-2))
     &                -4*(work(i,j,k-1)+(work(i,j,k-1)))+6*work(i,j,k))
                     end if
                    else if (k.gt.(2-bo1).and.k.lt.(nz-1)) then
                     if (chr(i,j,k-1).ne.ex.and.chr(i,j,k+1).ne.ex.and.
     &                   chr(i,j,k+2).ne.ex)
     &               then
                       f_hf=w3*(
     &                    -work(i,j,k-1)+3*work(i,j,k)
     &                  -3*work(i,j,k+1)+work(i,j,k+2))
                     end if
                    else if (k.gt.2.and.k.lt.(nz-1+bo1)) then
                     if (chr(i,j,k-1).ne.ex.and.chr(i,j,k+1).ne.ex.and.
     &                   chr(i,j,k-2).ne.ex)
     &               then
                       f_hf=w3*(
     &                    -work(i,j,k+1)+3*work(i,j,k)
     &                  -3*work(i,j,k-1)+work(i,j,k-2))
                     end if
                    else if ((k.gt.(2-bo2).or.
     &                       (k.eq.2.and.chr(i,j,1).eq.ex))
     &                       .and.k.lt.(nz-1).and.
     &                  (chr(i,j,k+1).ne.ex.and.chr(i,j,k+2).ne.ex) )
     &              then
                       f_hf=w2*(
     &                     work(i,j,k)-2*work(i,j,k+1)+work(i,j,k+2))
                    else if (k.gt.2.and.(k.lt.(nz-1+bo2).or.
     &                       (k.eq.(nz-1).and.chr(i,j,nz).eq.ex)).and.
     &                  (chr(i,j,k-1).ne.ex.and.chr(i,j,k-2).ne.ex) )
     &              then
                       f_hf=w2*(
     &                     work(i,j,k)-2*work(i,j,k-1)+work(i,j,k-2))
                    end if

                    f(i,j,k)=f(i,j,k)-eps_eff*f_hf
                   end if

!LR commented this out
!                   if (.not.ind_sweeps.or.pass.eq.3) then
!                    f_hf=0
!                    if (k.gt.2.and.k.lt.(nz-1).and.
!     &                  ((chr(i,j,k-2).ne.ex.and.chr(i,j,k-1).ne.ex.and.
!     &                    chr(i,j,k+2).ne.ex.and.chr(i,j,k+1).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                     work(i,j,k-2)+work(i,j,k+2)
!     &                 -4*(work(i,j,k-1)+work(i,j,k+1))+6*work(i,j,k))
!                    else if (k.eq.2.and.phys_bdy_type(5).eq.odd.and.
!     &                  ((chr(i,j,k-1).ne.ex.and.
!     &                    chr(i,j,k+2).ne.ex.and.chr(i,j,k+1).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                    (-work(i,j,k))+work(i,j,k+2)
!     &                 -4*(work(i,j,k-1)+work(i,j,k+1))+6*work(i,j,k))
!                    else if (k.eq.1.and.phys_bdy_type(5).eq.odd.and.
!     &                  ((chr(i,j,k+2).ne.ex.and.chr(i,j,k+1).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                    (-work(i,j,k+2))+work(i,j,k+2)
!     &                -4*((-work(i,j,k+1))+work(i,j,k+1))+6*work(i,j,k))
!                    else if (k.eq.Nz-1.and.phys_bdy_type(6).eq.odd.and.
!     &                  ((chr(i,j,k-2).ne.ex.and.chr(i,j,k-1).ne.ex.and.
!     &                    chr(i,j,k+1).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                     work(i,j,k-2)+(-work(i,j,k))
!     &                 -4*(work(i,j,k-1)+work(i,j,k+1))+6*work(i,j,k))
!                    else if (k.eq.Nz.and.phys_bdy_type(6).eq.odd.and.
!     &                  ((chr(i,j,k-2).ne.ex.and.chr(i,j,k-1).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                    work(i,j,k-2)+(-work(i,j,k-2))
!     &                -4*(work(i,j,k-1)+(-work(i,j,k-1)))+6*work(i,j,k))
!                    else if (k.eq.2.and.phys_bdy_type(5).eq.even.and.
!     &                  ((chr(i,j,k-1).ne.ex.and.
!     &                    chr(i,j,k+2).ne.ex.and.chr(i,j,k+1).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                    (work(i,j,k))+work(i,j,k+2)
!     &                 -4*(work(i,j,k-1)+work(i,j,k+1))+6*work(i,j,k))
!                    else if (k.eq.1.and.phys_bdy_type(5).eq.even.and.
!     &                  ((chr(i,j,k+2).ne.ex.and.chr(i,j,k+1).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                    (work(i,j,k+2))+work(i,j,k+2)
!     &                -4*((work(i,j,k+1))+work(i,j,k+1))+6*work(i,j,k))
!                    else if (k.eq.Nz-1.and.phys_bdy_type(6).eq.even.and.
!     &                  ((chr(i,j,k-2).ne.ex.and.chr(i,j,k-1).ne.ex.and.
!     &                    chr(i,j,k+1).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                     work(i,j,k-2)+(work(i,j,k))
!     &                 -4*(work(i,j,k-1)+work(i,j,k+1))+6*work(i,j,k))
!                    else if (k.eq.Nz.and.phys_bdy_type(6).eq.even.and.
!     &                  ((chr(i,j,k-2).ne.ex.and.chr(i,j,k-1).ne.ex)) )
!     &              then
!                       f_hf=w4*(
!     &                    work(i,j,k-2)+(work(i,j,k-2))
!     &                -4*(work(i,j,k-1)+(work(i,j,k-1)))+6*work(i,j,k))
!                    else if (k.gt.(2-bo1).and.k.lt.(nz-1).and.
!     &                  (chr(i,j,k-1).ne.ex.and.chr(i,j,k+1).ne.ex.and.
!     &                   chr(i,j,k+2).ne.ex)) 
!     &              then
!                       f_hf=w3*(
!     &                    -work(i,j,k-1)+3*work(i,j,k)
!     &                  -3*work(i,j,k+1)+work(i,j,k+2))
!                    else if (k.gt.2.and.k.lt.(nz-1+bo1).and.
!     &                  (chr(i,j,k-1).ne.ex.and.chr(i,j,k+1).ne.ex.and.
!     &                   chr(i,j,k-2).ne.ex) ) 
!     &              then
!                       f_hf=w3*(
!     &                    -work(i,j,k+1)+3*work(i,j,k)
!     &                  -3*work(i,j,k-1)+work(i,j,k-2))
!                    else if ((k.gt.(2-bo2).or.
!     &                       (k.eq.2.and.chr(i,j,1).eq.ex))
!     &                       .and.k.lt.(nz-1).and.
!     &                  (chr(i,j,k+1).ne.ex.and.chr(i,j,k+2).ne.ex) )
!     &              then
!                       f_hf=w2*(
!     &                     work(i,j,k)-2*work(i,j,k+1)+work(i,j,k+2))
!                    else if (k.gt.2.and.(k.lt.(nz-1+bo2).or.
!     &                       (k.eq.(nz-1).and.chr(i,j,nz).eq.ex)).and.
!     &                  (chr(i,j,k-1).ne.ex.and.chr(i,j,k-2).ne.ex) )
!     &              then
!                       f_hf=w2*(
!     &                     work(i,j,k)-2*work(i,j,k-1)+work(i,j,k-2))
!                    end if
!
!                    f(i,j,k)=f(i,j,k)-eps_eff*f_hf
!                   end if

                 end if
              end do
           end do
         end do
        end do

        return
        end
c---------------------------------------------------------------------------
c applies a simple KO filter to f, where mask != mask_off, and modified
c near excision boundaries as follows:
c (or all boundaries if do_bdy=1, 
c  or neither if do_bdy=0 and internal flag no_uc_ex_bdy=.true.)
c
c with undivided operators Up(f) = f(i+1)-f(i)
c                          Um(f) = f(i)-f(i-1)
c
c Interior: eps*w4*(Up Um)^2
c left+1  : eps*w3*(Up Up - 2 Up Um)
c left    : eps*w2*(Up Um)
c right-1 : eps*w3*(Um Um - 2 Up Um)
c right   : eps*w2*(Up Um)
c
c NOTE: update DV version too
c
c TEST VERSION
c---------------------------------------------------------------------------
        subroutine dmdiss3d_ex_X(f,work,eps,do_bdy,phys_bdy_type,even,
     &                         odd,mask,mask_off,nx,ny,
     &                         nz,chr,ex,ind_sweeps)
        implicit none
        integer nx,ny,nz,do_bdy,phys_bdy_type(6),even,odd
        real*8 f(nx,ny,nz),work(nx,ny,nz),mask(nx,ny,nz),eps,mask_off
        real*8 chr(nx,ny,nz),ex
        logical ind_sweeps

        integer i,j,k,bo1,bo2
        real*8 eps_eff,f_hf,norm_f
        real*8 w4,w3,w2
c        parameter (w4=1.0d0/16.0d0,w3=1.0d0/16.0d0,w2=1.0d0/16.0d0)
c!        parameter (w4=1.0d0/16.0d0,w3=1.0d0/8.0d0,w2=1.0d0/4.0d0)
        parameter (w4=1.0d0/16.0d0,w3=1.0d0/16.0d0,w2=1.0d0/8.0d0)
        integer pass,npass
        logical no_uc_ex_diss
        parameter (no_uc_ex_diss=.true.)
        
        eps_eff=eps

        if (do_bdy.eq.0) then
           bo1=0
           bo2=0
           if (no_uc_ex_diss) then
              bo1=-max(nx,ny,nz)
              bo2=bo1
           end if
        else
           bo1=1
           bo2=2
        end if

        npass=1
        if (ny.gt.1.and.ind_sweeps) npass=2
        if (nz.gt.1.and.ind_sweeps) npass=3

        do pass=1,npass

         do i=1,nx
           do j=1,ny
              do k=1,nz
                 work(i,j,k)=f(i,j,k)
              end do
           end do
         end do

         do i=1,nx
           do j=1,ny
              do k=1,nz
                 if (mask(i,j,k).ne.mask_off.and.
     &               (chr(i,j,k).ne.ex)) then

                   if (.not.ind_sweeps.or.pass.eq.1) then
                    f_hf=0
                    if (i.gt.2.and.i.lt.(nx-1).and.
     &                  ((chr(i-2,j,k).ne.ex.and.chr(i-1,j,k).ne.ex.and.
     &                    chr(i+2,j,k).ne.ex.and.chr(i+1,j,k).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                     work(i-2,j,k)+work(i+2,j,k)
     &                 -4*(work(i-1,j,k)+work(i+1,j,k))+6*work(i,j,k))
                    else if (i.eq.2.and.phys_bdy_type(1).eq.odd.and.
     &                  ((chr(i-1,j,k).ne.ex.and.
     &                    chr(i+2,j,k).ne.ex.and.chr(i+1,j,k).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                    (-work(i,j,k))+work(i+2,j,k)
     &                 -4*(work(i-1,j,k)+work(i+1,j,k))+6*work(i,j,k))
                    else if (i.eq.1.and.phys_bdy_type(1).eq.odd.and.
     &                 ((chr(i+2,j,k).ne.ex.and.chr(i+1,j,k).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                    (-work(i+2,j,k))+work(i+2,j,k)
     &                -4*((-work(i+1,j,k))+work(i+1,j,k))+6*work(i,j,k))
                    else if (i.eq.Nx-1.and.phys_bdy_type(2).eq.odd.and.
     &                  ((chr(i+1,j,k).ne.ex.and.
     &                    chr(i-2,j,k).ne.ex.and.chr(i-1,j,k).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                     work(i-2,j,k)+(-work(i,j,k))
     &                 -4*(work(i-1,j,k)+work(i+1,j,k))+6*work(i,j,k))
                    else if (i.eq.Nx.and.phys_bdy_type(2).eq.odd.and.
     &                 ((chr(i-2,j,k).ne.ex.and.chr(i-1,j,k).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                    work(i-2,j,k)+(-work(i-2,j,k))
     &                -4*(work(i-1,j,k)+(-work(i-1,j,k)))+6*work(i,j,k))
                    else if (i.eq.2.and.phys_bdy_type(1).eq.even.and.
     &                  ((chr(i-1,j,k).ne.ex.and.
     &                    chr(i+2,j,k).ne.ex.and.chr(i+1,j,k).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                    (work(i,j,k))+work(i+2,j,k)
     &                 -4*(work(i-1,j,k)+work(i+1,j,k))+6*work(i,j,k))
                    else if (i.eq.1.and.phys_bdy_type(1).eq.even.and.
     &                 ((chr(i+2,j,k).ne.ex.and.chr(i+1,j,k).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                    (work(i+2,j,k))+work(i+2,j,k)
     &                -4*((work(i+1,j,k))+work(i+1,j,k))+6*work(i,j,k))
                    else if (i.eq.Nx-1.and.phys_bdy_type(2).eq.even.and.
     &                  ((chr(i+1,j,k).ne.ex.and.
     &                    chr(i-2,j,k).ne.ex.and.chr(i-1,j,k).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                     work(i-2,j,k)+(work(i,j,k))
     &                 -4*(work(i-1,j,k)+work(i+1,j,k))+6*work(i,j,k))
                    else if (i.eq.Nx.and.phys_bdy_type(2).eq.even.and.
     &                 ((chr(i-2,j,k).ne.ex.and.chr(i-1,j,k).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                    work(i-2,j,k)+(work(i-2,j,k))
     &                -4*(work(i-1,j,k)+(work(i-1,j,k)))+6*work(i,j,k))
                    else if (i.gt.(2-bo1).and.i.lt.(nx-1).and.
     &                  (chr(i-1,j,k).ne.ex.and.chr(i+1,j,k).ne.ex.and.
     &                   chr(i+2,j,k).ne.ex)) 
     &              then
                       f_hf=w3*(
     &                  -2*work(i-1,j,k)+5*work(i,j,k)
     &                  -4*work(i+1,j,k)+work(i+2,j,k))
                    else if (i.gt.2.and.i.lt.(nx-1+bo1).and.
     &                  (chr(i-1,j,k).ne.ex.and.chr(i+1,j,k).ne.ex.and.
     &                   chr(i-2,j,k).ne.ex) ) 
     &              then
                       f_hf=w3*(
     &                   -2*work(i+1,j,k)+5*work(i,j,k)
     &                   -4*work(i-1,j,k)+work(i-2,j,k))
                    else if ((i.gt.(2-bo2).or.
     &                       (i.eq.2.and.chr(1,j,k).eq.ex))
     &                       .and.i.lt.(nx-1).and.
     &                  (chr(i+1,j,k).ne.ex.and.chr(i+2,j,k).ne.ex) )
     &              then
                       f_hf=w2*(
     &                     work(i,j,k)-2*work(i+1,j,k)+work(i+2,j,k))
                    else if (i.gt.2.and.(i.lt.(nx-1+bo2).or.
     &                       (i.eq.(nx-1).and.chr(nx,j,k).eq.ex)).and.
     &                  (chr(i-1,j,k).ne.ex.and.chr(i-2,j,k).ne.ex) )
     &              then
                       f_hf=w2*(
     &                     work(i,j,k)-2*work(i-1,j,k)+work(i-2,j,k))
                    end if

                    f(i,j,k)=f(i,j,k)-eps_eff*f_hf
                   end if

                   if (.not.ind_sweeps.or.pass.eq.2) then
                    f_hf=0
                    if (j.gt.2.and.j.lt.(ny-1).and.
     &                  ((chr(i,j-2,k).ne.ex.and.chr(i,j-1,k).ne.ex.and.
     &                    chr(i,j+2,k).ne.ex.and.chr(i,j+1,k).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                     work(i,j-2,k)+work(i,j+2,k)
     &                 -4*(work(i,j-1,k)+work(i,j+1,k))+6*work(i,j,k))
                    else if (j.eq.2.and.phys_bdy_type(3).eq.odd.and.
     &                  ((chr(i,j-1,k).ne.ex.and.
     &                    chr(i,j+2,k).ne.ex.and.chr(i,j+1,k).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                    (-work(i,j,k))+work(i,j+2,k)
     &                 -4*(work(i,j-1,k)+work(i,j+1,k))+6*work(i,j,k))
                    else if (j.eq.1.and.phys_bdy_type(3).eq.odd.and.
     &                  ((chr(i,j+2,k).ne.ex.and.chr(i,j+1,k).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                    (-work(i,j+2,k))+work(i,j+2,k)
     &                -4*((-work(i,j+1,k))+work(i,j+1,k))+6*work(i,j,k))
                    else if (j.eq.Ny-1.and.phys_bdy_type(4).eq.odd.and.
     &                  ((chr(i,j+1,k).ne.ex.and.
     &                    chr(i,j-2,k).ne.ex.and.chr(i,j-1,k).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                     work(i,j-2,k)+(-work(i,j,k))
     &                 -4*(work(i,j-1,k)+work(i,j+1,k))+6*work(i,j,k))
                    else if (j.eq.Ny.and.phys_bdy_type(4).eq.odd.and.
     &                  ((chr(i,j-2,k).ne.ex.and.chr(i,j-1,k).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                    work(i,j-2,k)+(-work(i,j-2,k))
     &                -4*(work(i,j-1,k)+(-work(i,j-1,k)))+6*work(i,j,k))
                    else if (j.eq.2.and.phys_bdy_type(3).eq.even.and.
     &                  ((chr(i,j-1,k).ne.ex.and.
     &                    chr(i,j+2,k).ne.ex.and.chr(i,j+1,k).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                    (work(i,j,k))+work(i,j+2,k)
     &                 -4*(work(i,j-1,k)+work(i,j+1,k))+6*work(i,j,k))
                    else if (j.eq.1.and.phys_bdy_type(3).eq.even.and.
     &                  ((chr(i,j+2,k).ne.ex.and.chr(i,j+1,k).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                    (work(i,j+2,k))+work(i,j+2,k)
     &                -4*((work(i,j+1,k))+work(i,j+1,k))+6*work(i,j,k))
                    else if (j.eq.Ny-1.and.phys_bdy_type(4).eq.even.and.
     &                  ((chr(i,j+1,k).ne.ex.and.
     &                    chr(i,j-2,k).ne.ex.and.chr(i,j-1,k).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                     work(i,j-2,k)+(work(i,j,k))
     &                 -4*(work(i,j-1,k)+work(i,j+1,k))+6*work(i,j,k))
                    else if (j.eq.Ny.and.phys_bdy_type(4).eq.even.and.
     &                  ((chr(i,j-2,k).ne.ex.and.chr(i,j-1,k).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                    work(i,j-2,k)+(work(i,j-2,k))
     &                -4*(work(i,j-1,k)+(work(i,j-1,k)))+6*work(i,j,k))
                    else if (j.gt.(2-bo1).and.j.lt.(ny-1).and.
     &                  (chr(i,j-1,k).ne.ex.and.chr(i,j+1,k).ne.ex.and.
     &                   chr(i,j+2,k).ne.ex)) 
     &              then
                       f_hf=w3*(
     &                   -2*work(i,j-1,k)+5*work(i,j,k)
     &                   -4*work(i,j+1,k)+work(i,j+2,k))
                    else if (j.gt.2.and.j.lt.(ny-1+bo1).and.
     &                  (chr(i,j-1,k).ne.ex.and.chr(i,j+1,k).ne.ex.and.
     &                   chr(i,j-2,k).ne.ex) ) 
     &              then
                       f_hf=w3*(
     &                   -2*work(i,j+1,k)+5*work(i,j,k)
     &                   -4*work(i,j-1,k)+work(i,j-2,k))
                    else if ((j.gt.(2-bo2).or.
     &                       (j.eq.2.and.chr(i,1,k).eq.ex))
     &                       .and.j.lt.(ny-1).and.
     &                  (chr(i,j+1,k).ne.ex.and.chr(i,j+2,k).ne.ex) )
     &              then
                       f_hf=w2*(
     &                     work(i,j,k)-2*work(i,j+1,k)+work(i,j+2,k))
                    else if (j.gt.2.and.(j.lt.(ny-1+bo2).or.
     &                       (j.eq.(ny-1).and.chr(i,ny,k).eq.ex)).and.
     &                  (chr(i,j-1,k).ne.ex.and.chr(i,j-2,k).ne.ex) )
     &              then
                       f_hf=w2*(
     &                     work(i,j,k)-2*work(i,j-1,k)+work(i,j-2,k))
                    end if

                    f(i,j,k)=f(i,j,k)-eps_eff*f_hf
                   end if

                   if (.not.ind_sweeps.or.pass.eq.3) then
                    f_hf=0
                    if (k.gt.2.and.k.lt.(nz-1).and.
     &                  ((chr(i,j,k-2).ne.ex.and.chr(i,j,k-1).ne.ex.and.
     &                    chr(i,j,k+2).ne.ex.and.chr(i,j,k+1).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                     work(i,j,k-2)+work(i,j,k+2)
     &                 -4*(work(i,j,k-1)+work(i,j,k+1))+6*work(i,j,k))
                    else if (k.eq.2.and.phys_bdy_type(5).eq.odd.and.
     &                  ((chr(i,j,k-1).ne.ex.and.
     &                    chr(i,j,k+2).ne.ex.and.chr(i,j,k+1).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                    (-work(i,j,k))+work(i,j,k+2)
     &                 -4*(work(i,j,k-1)+work(i,j,k+1))+6*work(i,j,k))
                    else if (k.eq.1.and.phys_bdy_type(5).eq.odd.and.
     &                  ((chr(i,j,k+2).ne.ex.and.chr(i,j,k+1).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                    (-work(i,j,k+2))+work(i,j,k+2)
     &                -4*((-work(i,j,k+1))+work(i,j,k+1))+6*work(i,j,k))
                    else if (k.eq.Nz-1.and.phys_bdy_type(6).eq.odd.and.
     &                  ((chr(i,j,k-2).ne.ex.and.chr(i,j,k-1).ne.ex.and.
     &                    chr(i,j,k+1).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                     work(i,j,k-2)+(-work(i,j,k))
     &                 -4*(work(i,j,k-1)+work(i,j,k+1))+6*work(i,j,k))
                    else if (k.eq.Nz.and.phys_bdy_type(6).eq.odd.and.
     &                  ((chr(i,j,k-2).ne.ex.and.chr(i,j,k-1).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                    work(i,j,k-2)+(-work(i,j,k-2))
     &                -4*(work(i,j,k-1)+(-work(i,j,k-1)))+6*work(i,j,k))
                    else if (k.eq.2.and.phys_bdy_type(5).eq.even.and.
     &                  ((chr(i,j,k-1).ne.ex.and.
     &                    chr(i,j,k+2).ne.ex.and.chr(i,j,k+1).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                    (work(i,j,k))+work(i,j,k+2)
     &                 -4*(work(i,j,k-1)+work(i,j,k+1))+6*work(i,j,k))
                    else if (k.eq.1.and.phys_bdy_type(5).eq.even.and.
     &                  ((chr(i,j,k+2).ne.ex.and.chr(i,j,k+1).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                    (work(i,j,k+2))+work(i,j,k+2)
     &                -4*((work(i,j,k+1))+work(i,j,k+1))+6*work(i,j,k))
                    else if (k.eq.Nz-1.and.phys_bdy_type(6).eq.even.and.
     &                  ((chr(i,j,k-2).ne.ex.and.chr(i,j,k-1).ne.ex.and.
     &                    chr(i,j,k+1).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                     work(i,j,k-2)+(work(i,j,k))
     &                 -4*(work(i,j,k-1)+work(i,j,k+1))+6*work(i,j,k))
                    else if (k.eq.Nz.and.phys_bdy_type(6).eq.even.and.
     &                  ((chr(i,j,k-2).ne.ex.and.chr(i,j,k-1).ne.ex)) )
     &              then
                       f_hf=w4*(
     &                    work(i,j,k-2)+(work(i,j,k-2))
     &                -4*(work(i,j,k-1)+(work(i,j,k-1)))+6*work(i,j,k))
                    else if (k.gt.(2-bo1).and.k.lt.(nz-1).and.
     &                  (chr(i,j,k-1).ne.ex.and.chr(i,j,k+1).ne.ex.and.
     &                   chr(i,j,k+2).ne.ex)) 
     &              then
                       f_hf=w3*(
     &                  -2*work(i,j,k-1)+5*work(i,j,k)
     &                  -4*work(i,j,k+1)+work(i,j,k+2))
                    else if (k.gt.2.and.k.lt.(nz-1+bo1).and.
     &                  (chr(i,j,k-1).ne.ex.and.chr(i,j,k+1).ne.ex.and.
     &                   chr(i,j,k-2).ne.ex) ) 
     &              then
                       f_hf=w3*(
     &                   -2*work(i,j,k+1)+5*work(i,j,k)
     &                   -4*work(i,j,k-1)+work(i,j,k-2))
                    else if ((k.gt.(2-bo2).or.
     &                       (k.eq.2.and.chr(i,j,1).eq.ex))
     &                       .and.k.lt.(nz-1).and.
     &                  (chr(i,j,k+1).ne.ex.and.chr(i,j,k+2).ne.ex) )
     &              then
                       f_hf=w2*(
     &                     work(i,j,k)-2*work(i,j,k+1)+work(i,j,k+2))
                    else if (k.gt.2.and.(k.lt.(nz-1+bo2).or.
     &                       (k.eq.(nz-1).and.chr(i,j,nz).eq.ex)).and.
     &                  (chr(i,j,k-1).ne.ex.and.chr(i,j,k-2).ne.ex) )
     &              then
                       f_hf=w2*(
     &                     work(i,j,k)-2*work(i,j,k-1)+work(i,j,k-2))
                    end if

                    f(i,j,k)=f(i,j,k)-eps_eff*f_hf
                   end if

                 end if
              end do
           end do
         end do
        end do

        return
        end

c---------------------------------------------------------------------------
c applies a simple KO filter to f, where mask != mask_off.
c if do_bdy, then stencil modified as with excision below near boundaries
c
c if (do_ex) then chr is an excision mask ... and for now call
c dmdiss3d_ex
c 
c   NOTE: SPECIAL FLAG WITH EXCISION: 
c         if do_ex < 0, then eps is treated as an array, 
c         else if do_ex > 0, eps is treated as a scalar
c
c phys_bdy_type specifies wether a given boundary is even, odd (or something
c else, which means the flag is ignored)
c
c NOTE: update DV version too 
c---------------------------------------------------------------------------
        subroutine dmdiss3d(f,work,eps,do_bdy,phys_bdy_type,even,odd,
     &                      mask,mask_off,nx,ny,nz,chr,ex,do_ex)
        implicit none
        integer nx,ny,nz,do_ex,do_bdy,phys_bdy_type(6),even,odd
        real*8 f(nx,ny,nz),work(nx,ny,nz),mask(nx,ny,nz),eps(nx,ny,nz)
        real*8 mask_off
        real*8 chr(nx,ny,nz),ex

        integer i,j,k
        real*8 eps_eff,f_hf,norm_f,eps0
        real*8 w4,w3,w2
c        parameter (w4=1.0d0/16.0d0,w3=1.0d0/16.0d0,w2=1.0d0/16.0d0)
        parameter (w4=1.0d0/16.0d0,w3=1.0d0/8.0d0,w2=1.0d0/4.0d0)
        integer pass,npass,cdo_bdy(6)
        logical ind_sweeps
        parameter (ind_sweeps=.true.)

        if (do_ex.ne.0) then
           call dmdiss3d_ex(f,work,eps,do_bdy,phys_bdy_type,even,odd,
     &          mask,mask_off,nx,ny,nz,chr,ex,do_ex,ind_sweeps)
           return
        end if

        do i=1,6
           cdo_bdy(i)=do_bdy
           if (phys_bdy_type(i).eq.even.or.phys_bdy_type(i).eq.odd)
     &        cdo_bdy(i)=1
        end do

        eps_eff=eps(1,1,1)

        npass=1
        if (ny.gt.1.and.ind_sweeps) npass=2
        if (nz.gt.1.and.ind_sweeps) npass=3

        do pass=1,npass

         do i=1,nx
           do j=1,ny
              do k=1,nz
                 work(i,j,k)=f(i,j,k)
              end do
           end do
         end do

         do i=1,nx
           do j=1,ny
              do k=1,nz
                 if (mask(i,j,k).ne.mask_off) then
                   if (.not.ind_sweeps.or.pass.eq.1) then
                    f_hf=0
                    if (i.gt.2.and.i.lt.(nx-1)) then
                       f_hf=w4*(
     &                     work(i-2,j,k)+work(i+2,j,k)
     &                 -4*(work(i-1,j,k)+work(i+1,j,k))+6*work(i,j,k))
                    else 
                     if (cdo_bdy(1).ne.0.and.Nx.ge.5) then
                      if (i.eq.2.and.phys_bdy_type(1).eq.odd) then
                       f_hf=w4*(
     &                    (-work(i,j,k))+work(i+2,j,k)
     &                 -4*(work(i-1,j,k)+work(i+1,j,k))+6*work(i,j,k))
                      else if (i.eq.2.and.phys_bdy_type(1).eq.even) then
                       f_hf=w4*(
     &                    (work(i,j,k))+work(i+2,j,k)
     &                 -4*(work(i-1,j,k)+work(i+1,j,k))+6*work(i,j,k))
                      else if (i.eq.2) then
                       f_hf=w3*(
     &                    -work(i-1,j,k)+3*work(i,j,k)
     &                  -3*work(i+1,j,k)+work(i+2,j,k))
                      else if (i.eq.1.and.phys_bdy_type(1).eq.odd) then
                       f_hf=w4*(
     &                    (-work(i+2,j,k))+work(i+2,j,k)
     &                -4*((-work(i+1,j,k))+work(i+1,j,k))+6*work(i,j,k))
                      else if (i.eq.1.and.phys_bdy_type(1).eq.even) then
                       f_hf=w4*(
     &                     (work(i+2,j,k))+work(i+2,j,k)
     &                 -4*((work(i+1,j,k))+work(i+1,j,k))+6*work(i,j,k))
                      else if (i.eq.1) then
                       f_hf=w2*(
     &                     work(i,j,k)-2*work(i+1,j,k)+work(i+2,j,k))
                      end if
                     end if
                     if (cdo_bdy(2).ne.0.and.Nx.ge.5) then
                      if (i.eq.(nx-1).and.phys_bdy_type(2).eq.odd) then
                       f_hf=w4*(
     &                     work(i-2,j,k)+(-work(i,j,k))
     &                 -4*(work(i-1,j,k)+work(i+1,j,k))+6*work(i,j,k))
                      else if (i.eq.(nx-1).and.phys_bdy_type(2).eq.even) 
     &                then
                       f_hf=w4*(
     &                     work(i-2,j,k)+(work(i,j,k))
     &                 -4*(work(i-1,j,k)+work(i+1,j,k))+6*work(i,j,k))
                      else if (i.eq.(nx-1)) then
                       f_hf=w3*(
     &                    -work(i+1,j,k)+3*work(i,j,k)
     &                  -3*work(i-1,j,k)+work(i-2,j,k))
                      else if (i.eq.nx.and.phys_bdy_type(2).eq.odd) then
                       f_hf=w4*(
     &                     work(i-2,j,k)+(-work(i-2,j,k))
     &                -4*(work(i-1,j,k)+(-work(i-1,j,k)))+6*work(i,j,k))
                      else if (i.eq.nx.and.phys_bdy_type(2).eq.even) 
     &                then
                       f_hf=w4*(
     &                     work(i-2,j,k)+(work(i-2,j,k))
     &                 -4*(work(i-1,j,k)+(work(i-1,j,k)))+6*work(i,j,k))
                      else if (i.eq.nx) then
                       f_hf=w2*(
     &                     work(i,j,k)-2*work(i-1,j,k)+work(i-2,j,k))
                      end if
                     end if
                    end if

                    f(i,j,k)=f(i,j,k)-eps_eff*f_hf
                   end if

                   if (.not.ind_sweeps.or.pass.eq.2) then
                    f_hf=0
                    if (j.gt.2.and.j.lt.(ny-1)) then
                       f_hf=w4*(
     &                     work(i,j-2,k)+work(i,j+2,k)
     &                 -4*(work(i,j-1,k)+work(i,j+1,k))+6*work(i,j,k))
                    else 
                     if (cdo_bdy(3).ne.0.and.Ny.ge.5) then
                      if (j.eq.2.and.phys_bdy_type(3).eq.odd) then
                       f_hf=w4*(
     &                     (-work(i,j,k))+work(i,j+2,k)
     &                 -4*(work(i,j-1,k)+work(i,j+1,k))+6*work(i,j,k))
                      else if (j.eq.2.and.phys_bdy_type(3).eq.even) then
                       f_hf=w4*(
     &                    (work(i,j,k))+work(i,j+2,k)
     &                 -4*(work(i,j-1,k)+work(i,j+1,k))+6*work(i,j,k))
                      else if (j.eq.2) then
                       f_hf=w3*(
     &                    -work(i,j-1,k)+3*work(i,j,k)
     &                  -3*work(i,j+1,k)+work(i,j+2,k))
                      else if (j.eq.1.and.phys_bdy_type(3).eq.odd) then
                       f_hf=w4*(
     &                    (-work(i,j+2,k))+work(i,j+2,k)
     &                -4*((-work(i,j+1,k))+work(i,j+1,k))+6*work(i,j,k))
                      else if (j.eq.1.and.phys_bdy_type(3).eq.even) then
                       f_hf=w4*(
     &                     (work(i,j+2,k))+work(i,j+2,k)
     &                 -4*((work(i,j+1,k))+work(i,j+1,k))+6*work(i,j,k))
                      else if (j.eq.1) then
                       f_hf=w2*(
     &                     work(i,j,k)-2*work(i,j+1,k)+work(i,j+2,k))
                      end if
                     end if
                     if (cdo_bdy(4).ne.0.and.Ny.ge.5) then
                      if (j.eq.(ny-1).and.phys_bdy_type(4).eq.odd) then
                       f_hf=w4*(
     &                     work(i,j-2,k)+(-work(i,j,k))
     &                 -4*(work(i,j-1,k)+work(i,j+1,k))+6*work(i,j,k))
                      else if (j.eq.(ny-1).and.phys_bdy_type(4).eq.even)
     &                then
                       f_hf=w4*(
     &                     work(i,j-2,k)+(work(i,j,k))
     &                 -4*(work(i,j-1,k)+work(i,j+1,k))+6*work(i,j,k))
                      else if (j.eq.(ny-1)) then
                       f_hf=w3*(
     &                    -work(i,j+1,k)+3*work(i,j,k)
     &                  -3*work(i,j-1,k)+work(i,j-2,k))
                      else if (j.eq.ny.and.phys_bdy_type(4).eq.odd) then
                       f_hf=w4*(
     &                    work(i,j-2,k)+(-work(i,j-2,k))
     &                -4*(work(i,j-1,k)+(-work(i,j-1,k)))+6*work(i,j,k))
                      else if (j.eq.ny.and.phys_bdy_type(4).eq.even) 
     &                then
                       f_hf=w4*(
     &                     work(i,j-2,k)+(work(i,j-2,k))
     &                 -4*(work(i,j-1,k)+(work(i,j-1,k)))+6*work(i,j,k))
                      else if (j.eq.ny) then
                       f_hf=w2*(
     &                     work(i,j,k)-2*work(i,j-1,k)+work(i,j-2,k))
                      end if
                     end if
                    end if
         
                    f(i,j,k)=f(i,j,k)-eps_eff*f_hf
                   end if

                   if (.not.ind_sweeps.or.pass.eq.3) then
                    f_hf=0
                    if (k.gt.2.and.k.lt.(nz-1)) then
                       f_hf=w4*(
     &                     work(i,j,k-2)+work(i,j,k+2)
     &                 -4*(work(i,j,k-1)+work(i,j,k+1))+6*work(i,j,k))
                    else 
                     if (cdo_bdy(5).ne.0.and.Nz.ge.5) then
                      if (k.eq.2.and.phys_bdy_type(5).eq.odd) then
                       f_hf=w4*(
     &                    (-work(i,j,k))+work(i,j,k+2)
     &                 -4*(work(i,j,k-1)+work(i,j,k+1))+6*work(i,j,k))
                      else if (k.eq.2.and.phys_bdy_type(5).eq.even) then
                       f_hf=w4*(
     &                    (work(i,j,k))+work(i,j,k+2)
     &                 -4*(work(i,j,k-1)+work(i,j,k+1))+6*work(i,j,k))
                      else if (k.eq.2) then
                       f_hf=w3*(
     &                    -work(i,j,k-1)+3*work(i,j,k)
     &                  -3*work(i,j,k+1)+work(i,j,k+2))
                      else if (k.eq.1.and.phys_bdy_type(5).eq.odd) then
                       f_hf=w4*(
     &                    (-work(i,j,k+2))+work(i,j,k+2)
     &                -4*((-work(i,j,k+1))+work(i,j,k+1))+6*work(i,j,k))
                      else if (k.eq.1.and.phys_bdy_type(5).eq.even) then
                       f_hf=w4*(
     &                     (work(i,j,k+2))+work(i,j,k+2)
     &                 -4*((work(i,j,k+1))+work(i,j,k+1))+6*work(i,j,k))
                      else if (k.eq.1) then
                       f_hf=w2*(
     &                     work(i,j,k)-2*work(i,j,k+1)+work(i,j,k+2))
                      end if
                     end if
                     if (cdo_bdy(6).ne.0.and.Nz.ge.5) then
                      if (k.eq.(nz-1).and.phys_bdy_type(6).eq.odd) then
                       f_hf=w4*(
     &                     work(i,j,k-2)+(-work(i,j,k))
     &                 -4*(work(i,j,k-1)+work(i,j,k+1))+6*work(i,j,k))
                      else if (k.eq.(nz-1).and.phys_bdy_type(6).eq.even) 
     &                then
                       f_hf=w4*(
     &                     work(i,j,k-2)+(work(i,j,k))
     &                 -4*(work(i,j,k-1)+work(i,j,k+1))+6*work(i,j,k))
                      else if (k.eq.(nz-1)) then
                       f_hf=w3*(
     &                    -work(i,j,k+1)+3*work(i,j,k)
     &                  -3*work(i,j,k-1)+work(i,j,k-2))
                      else if (k.eq.nz.and.phys_bdy_type(6).eq.odd) then
                       f_hf=w4*(
     &                    work(i,j,k-2)+(-work(i,j,k-2))
     &                -4*(work(i,j,k-1)+(-work(i,j,k-1)))+6*work(i,j,k))
                      else if (k.eq.nz.and.phys_bdy_type(6).eq.even) 
     &                then
                       f_hf=w4*(
     &                     work(i,j,k-2)+(work(i,j,k-2))
     &                 -4*(work(i,j,k-1)+(work(i,j,k-1)))+6*work(i,j,k))
                      else if (k.eq.nz) then
                       f_hf=w2*(
     &                     work(i,j,k)-2*work(i,j,k-1)+work(i,j,k-2))
                      end if
                     end if
                    end if

                    f(i,j,k)=f(i,j,k)-eps_eff*f_hf
                   end if
                 end if
              end do
           end do
         end do
        end do

        return
        end

c-----------------------------------------------------------------------
c enforces linear approach to: which_bnd =
c
c 1 physical
c 2 AMR
c 3 both
c
c phys_bdy=1 .. physical
c         =0 .. AMR
c
c if (do_ex) then interprets ex as an excision mask
c-----------------------------------------------------------------------
        subroutine linbnd3d(f,phys_bdy,which_bnd,Nx,Ny,Nz,
     &                      chr,ex,do_ex)
        implicit none
        integer Nx,Ny,Nz,which_bnd,do_ex
        real*8 f(Nx,Ny,Nz),chr(Nx,Ny,Nz),ex
        integer phys_bdy(6)

        integer i,j,k,is,ie,js,je,ks,ke,pf
        logical all

        if (Nx.lt.3.or.Ny.eq.2.or.Nz.eq.2) then
           write(*,*) 'lin_bnd_3D error: invalid grid dimension'
           stop
        end if

        all=.false.
        if (which_bnd.eq.3) then
           all=.true.
        else if (which_bnd.eq.1) then
           pf=1
        else if (which_bnd.eq.2) then
           pf=0
        end if

        if (phys_bdy(1).eq.pf.or.all) then
           do j=min(Ny,2),max(1,Ny-1)
              do k=min(Nz,2),max(1,Nz-1)
                 if (do_ex.eq.0.or.
     &               (chr(3,j,k).ne.ex.and.chr(1,j,k).ne.ex))
     &              f(2,j,k)=0.5*(f(3,j,k)+f(1,j,k))
              end do
           end do
           if (Ny.gt.1.and.(phys_bdy(3).eq.pf.or.all)) then
              do k=min(Nz,2),max(1,Nz-1)
                 if (do_ex.eq.0.or.
     &               (chr(3,3,k).ne.ex.and.chr(3,1,k).ne.ex.and.
     &                chr(1,3,k).ne.ex.and.chr(1,1,k).ne.ex))
     &              f(2,2,k)=0.25*(f(3,3,k)+f(1,3,k)+f(3,1,k)+f(1,1,k))
              end do
           end if
           if (Ny.gt.1.and.(phys_bdy(4).eq.pf.or.all)) then
              do k=min(Nz,2),max(1,Nz-1)
                 if (do_ex.eq.0.or.
     &               (chr(3,Ny-2,k).ne.ex.and.chr(1,Ny-2,k).ne.ex.and.
     &                chr(3,Ny,k).ne.ex.and.chr(1,Ny,k).ne.ex))
     &              f(2,Ny-1,k)=
     &              0.25*(f(3,Ny-2,k)+f(1,Ny-2,k)+f(3,Ny,k)+f(1,Ny,k))
              end do
           end if
           if (Nz.gt.1.and.(phys_bdy(5).eq.pf.or.all)) then
              do j=2,Ny-1
                 if (do_ex.eq.0.or.
     &               (chr(3,j,3).ne.ex.and.chr(1,j,3).ne.ex.and.
     &                chr(3,j,1).ne.ex.and.chr(1,j,1).ne.ex))
     &              f(2,j,2)=0.25*(f(3,j,3)+f(1,j,3)+f(3,j,1)+f(1,j,1))
              end do
           end if
           if (Nz.gt.1.and.(phys_bdy(6).eq.pf.or.all)) then
              do j=2,Ny-1
                 if (do_ex.eq.0.or.
     &               (chr(3,j,Nz-2).ne.ex.and.chr(1,j,Nz-2).ne.ex.and.
     &                chr(3,j,Nz).ne.ex.and.chr(1,j,Nz).ne.ex))
     &              f(2,j,Nz-1)=0.25*(f(3,j,Nz-2)+f(1,j,Nz-2)+
     &                                f(3,j,Nz)+f(1,j,Nz))
              end do
           end if 

           if (Nz.gt.1.and.(
     &         (phys_bdy(3).eq.pf.and.phys_bdy(5).eq.pf).or.all.and.
     &         (do_ex.eq.0.or.
     &          (chr(3,3,3).ne.ex.and.chr(1,3,3).ne.ex.and.
     &           chr(3,1,3).ne.ex.and.chr(3,3,1).ne.ex.and.
     &           chr(1,1,3).ne.ex.and.chr(1,3,1).ne.ex.and.
     &           chr(3,1,1).ne.ex.and.chr(1,1,1).ne.ex))))
     &        f(2,2,2)=0.125*(f(3,3,3)+f(1,3,3)+f(3,1,3)+f(3,3,1)+
     &                        f(1,1,3)+f(1,3,1)+f(3,1,1)+f(1,1,1))
           if (Nz.gt.1.and.(
     &         (phys_bdy(4).eq.pf.and.phys_bdy(5).eq.pf).or.all.and.
     &         (do_ex.eq.0.or.
     &          (chr(3,Ny-2,3).ne.ex.and.chr(1,Ny-2,3).ne.ex.and.
     &           chr(3,Ny,3).ne.ex.and.chr(3,Ny-2,1).ne.ex.and.
     &           chr(1,Ny,3).ne.ex.and.chr(1,Ny-2,1).ne.ex.and.
     &           chr(3,Ny,1).ne.ex.and.chr(1,Ny,1).ne.ex))))
     &        f(2,Ny-1,2)=0.125*(f(3,Ny-2,3)+f(1,Ny-2,3)+f(3,Ny,3)+
     &                           f(3,Ny-2,1)+f(1,Ny,3)+f(1,Ny-2,1)+
     &                           f(3,Ny,1)+f(1,Ny,1))
           if (Nz.gt.1.and.(
     &         (phys_bdy(3).eq.pf.and.phys_bdy(6).eq.pf).or.all.and.
     &         (do_ex.eq.0.or.
     &          (chr(3,3,Nz-2).ne.ex.and.chr(1,3,Nz-2).ne.ex.and.
     &           chr(3,1,Nz-2).ne.ex.and.chr(3,3,Nz).ne.ex.and.
     &           chr(1,1,Nz-2).ne.ex.and.chr(1,3,Nz).ne.ex.and.
     &           chr(3,1,Nz).ne.ex.and.chr(1,1,Nz).ne.ex))))
     &        f(2,2,Nz-1)=0.125*(f(3,3,Nz-2)+f(1,3,Nz-2)+f(3,1,Nz-2)+
     &                           f(3,3,Nz)+f(1,1,Nz-2)+f(1,3,Nz)+
     &                           f(3,1,Nz)+f(1,1,Nz))
           if (Nz.gt.1.and.(
     &         (phys_bdy(4).eq.pf.and.phys_bdy(6).eq.pf).or.all.and.
     &         (do_ex.eq.0.or.
     &          (chr(3,Ny-2,Nz-2).ne.ex.and.chr(1,Ny-2,Nz-2).ne.ex.and.
     &           chr(3,Ny,Nz-2).ne.ex.and.chr(3,Ny-2,Nz).ne.ex.and.
     &           chr(1,Ny,Nz-2).ne.ex.and.chr(1,Ny-2,Nz).ne.ex.and.
     &           chr(3,Ny,Nz).ne.ex.and.chr(1,Ny,Nz).ne.ex))))
     &        f(2,Ny-1,Nz-1)=0.125*(f(3,Ny-2,Nz-2)+f(1,Ny-2,Nz-2)+
     &                              f(3,Ny,Nz-2)+f(3,Ny-2,Nz)+
     &                              f(1,Ny,Nz-2)+f(1,Ny-2,Nz)+
     &                              f(3,Ny,Nz)+f(1,Ny,Nz))
        end if

        if (phys_bdy(2).eq.pf.or.all) then
           do j=min(Ny,2),max(1,Ny-1)
              do k=min(Nz,2),max(1,Nz-1)
                 if (do_ex.eq.0.or.
     &               (chr(Nx-2,j,k).ne.ex.and.chr(Nx,j,k).ne.ex))
     &           f(Nx-1,j,k)=0.5*(f(Nx-2,j,k)+f(Nx,j,k))
              end do
           end do
           if (Ny.gt.1.and.(phys_bdy(3).eq.pf.or.all)) then
              do k=min(Nz,2),max(1,Nz-1)
                 if (do_ex.eq.0.or.
     &               (chr(Nx-2,3,k).ne.ex.and.chr(Nx,3,k).ne.ex.and.
     &                chr(Nx-2,1,k).ne.ex.and.chr(Nx,1,k).ne.ex))
     &           f(Nx-1,2,k)=0.25*(f(Nx-2,3,k)+f(Nx,3,k)+
     &                             f(Nx-2,1,k)+f(Nx,1,k))
              end do
           end if
           if (Ny.gt.1.and.(phys_bdy(4).eq.pf.or.all)) then
              do k=min(Nz,2),max(1,Nz-1)
                 if (do_ex.eq.0.or.
     &            (chr(Nx-2,Ny-2,k).ne.ex.and.chr(Nx,Ny-2,k).ne.ex.and.
     &             chr(Nx-2,Ny,k).ne.ex.and.chr(Nx,Ny,k).ne.ex))
     &           f(Nx-1,Ny-1,k)=0.25*(f(Nx-2,Ny-2,k)+f(Nx,Ny-2,k)+
     &                                f(Nx-2,Ny,k)+f(Nx,Ny,k))
              end do
           end if
           if (Nz.gt.1.and.(phys_bdy(5).eq.pf.or.all)) then
              do j=2,Ny-1
                 if (do_ex.eq.0.or.
     &            (chr(Nx-2,j,3).ne.ex.and.chr(Nx,j,3).ne.ex.and.
     &             chr(Nx-2,j,1).ne.ex.and.chr(Nx,j,1).ne.ex))
     &           f(Nx-1,j,2)=0.25*(f(Nx-2,j,3)+f(Nx,j,3)+
     &                             f(Nx-2,j,1)+f(Nx,j,1))
              end do
           end if
           if (Nz.gt.1.and.(phys_bdy(6).eq.pf.or.all)) then
              do j=2,Ny-1
                 if (do_ex.eq.0.or.
     &            (chr(Nx-2,j,Nz-2).ne.ex.and.chr(Nx,j,Nz-2).ne.ex.and.
     &             chr(Nx-2,j,Nz).ne.ex.and.chr(Nx,j,Nz).ne.ex))
     &           f(Nx-1,j,Nz-1)=0.25*(f(Nx-2,j,Nz-2)+f(Nx,j,Nz-2)+
     &                                f(Nx-2,j,Nz)+f(Nx,j,Nz))
              end do
           end if
           if (Nz.gt.1.and.(
     &         (phys_bdy(3).eq.pf.and.phys_bdy(5).eq.pf).or.all.and.
     &         (do_ex.eq.0.or.
     &          (chr(Nx-2,3,3).ne.ex.and.chr(Nx,3,3).ne.ex.and.
     &           chr(Nx-2,1,3).ne.ex.and.chr(Nx-2,3,1).ne.ex.and.
     &           chr(Nx,1,3).ne.ex.and.chr(Nx,3,1).ne.ex.and.
     &           chr(Nx-2,1,1).ne.ex.and.chr(Nx,1,1).ne.ex))))
     &        f(Nx-1,2,2)=0.125*(f(Nx-2,3,3)+f(Nx,3,3)+
     &                           f(Nx-2,1,3)+f(Nx-2,3,1)+
     &                           f(Nx,1,3)+f(Nx,3,1)+
     &                           f(Nx-2,1,1)+f(Nx,1,1))
           if (Nz.gt.1.and.(
     &         (phys_bdy(4).eq.pf.and.phys_bdy(5).eq.pf).or.all.and.
     &         (do_ex.eq.0.or.
     &          (chr(Nx-2,Ny-2,3).ne.ex.and.chr(Nx,Ny-2,3).ne.ex.and.
     &           chr(Nx-2,Ny,3).ne.ex.and.chr(Nx-2,Ny-2,1).ne.ex.and.
     &           chr(Nx,Ny,3).ne.ex.and.chr(Nx,Ny-2,1).ne.ex.and.
     &           chr(Nx-2,Ny,1).ne.ex.and.chr(Nx,Ny,1).ne.ex))))
     &        f(Nx-1,Ny-1,2)=0.125*(f(Nx-2,Ny-2,3)+f(Nx,Ny-2,3)+
     &                              f(Nx-2,Ny,3)+f(Nx-2,Ny-2,1)+
     &                              f(Nx,Ny,3)+f(Nx,Ny-2,1)+
     &                              f(Nx-2,Ny,1)+f(Nx,Ny,1))
           if (Nz.gt.1.and.(
     &         (phys_bdy(3).eq.pf.and.phys_bdy(6).eq.pf).or.all.and.
     &         (do_ex.eq.0.or.
     &          (chr(Nx-2,3,Nz-2).ne.ex.and.chr(Nx,3,Nz-2).ne.ex.and.
     &           chr(Nx-2,1,Nz-2).ne.ex.and.chr(Nx-2,3,Nz).ne.ex.and.
     &           chr(Nx,1,Nz-2).ne.ex.and.chr(Nx,3,Nz).ne.ex.and.
     &           chr(Nx-2,1,Nz).ne.ex.and.chr(Nx,1,Nz).ne.ex))))
     &        f(Nx-1,2,Nz-1)=0.125*(f(Nx-2,3,Nz-2)+f(Nx,3,Nz-2)+
     &                              f(Nx-2,1,Nz-2)+f(Nx-2,3,Nz)+
     &                              f(Nx,1,Nz-2)+f(Nx,3,Nz)+
     &                              f(Nx-2,1,Nz)+f(Nx,1,Nz))
           if (Nz.gt.1.and.(
     &         (phys_bdy(4).eq.pf.and.phys_bdy(6).eq.pf).or.all.and.
     &         (do_ex.eq.0.or.
     &     (chr(Nx-2,Ny-2,Nz-2).ne.ex.and.chr(Nx,Ny,Nz).ne.ex.and.
     &      chr(Nx,Ny-2,Nz-2).ne.ex.and.chr(Nx-2,Ny,Nz-2).ne.ex.and.
     &      chr(Nx-2,Ny-2,Nz).ne.ex.and.chr(Nx,Ny,Nz-2).ne.ex.and.
     &      chr(Nx,Ny-2,Nz).ne.ex.and.chr(Nx-2,Ny,Nz).ne.ex))))
     &        f(Nx-1,Ny-1,Nz-1)=0.125*(f(Nx-2,Ny-2,Nz-2)+f(Nx,Ny,Nz)+
     &                          f(Nx,Ny-2,Nz-2)+f(Nx-2,Ny,Nz-2)+
     &                          f(Nx-2,Ny-2,Nz)+f(Nx,Ny,Nz-2)+
     &                          f(Nx,Ny-2,Nz)+f(Nx-2,Ny,Nz))
        end if

        if (Ny.gt.1.and.(phys_bdy(3).eq.pf.or.all)) then
           do i=3,Nx-2
              do k=min(Nz,2),max(1,Nz-1)
                 if (do_ex.eq.0.or.
     &               (chr(i,3,k).ne.ex.and.chr(i,1,k).ne.ex))
     &           f(i,2,k)=0.5*(f(i,3,k)+f(i,1,k))
              end do
           end do
           if (Nz.gt.1.and.(phys_bdy(5).eq.pf.or.all)) then
              do i=3,Nx-2
                 if (do_ex.eq.0.or.
     &            (chr(i,3,3).ne.ex.and.chr(i,1,3).ne.ex.and.
     &             chr(i,3,1).ne.ex.and.chr(i,1,1).ne.ex))
     &           f(i,2,2)=0.25*(f(i,3,3)+f(i,1,3)+f(i,3,1)+f(i,1,1))
              end do
           end if
           if (Nz.gt.1.and.(phys_bdy(6).eq.pf.or.all)) then
              do i=3,Nx-2
                 if (do_ex.eq.0.or.
     &            (chr(i,3,Nz-2).ne.ex.and.chr(i,1,Nz-2).ne.ex.and.
     &             chr(i,3,Nz).ne.ex.and.chr(i,1,Nz).ne.ex))
     &           f(i,2,Nz-1)=0.25*(f(i,3,Nz-2)+f(i,1,Nz-2)+
     &                             f(i,3,Nz)+f(i,1,Nz))
              end do
           end if
        end if

        if (Ny.gt.1.and.(phys_bdy(4).eq.pf.or.all)) then
           do i=3,Nx-2
              do k=min(Nz,2),max(1,Nz-1)
                 if (do_ex.eq.0.or.
     &               (chr(i,Ny-2,k).ne.ex.and.chr(i,Ny,k).ne.ex))
     &           f(i,Ny-1,k)=0.5*(f(i,Ny-2,k)+f(i,Ny,k))
              end do
           end do
           if (Nz.gt.1.and.(phys_bdy(5).eq.pf.or.all)) then
              do i=3,Nx-2
                 if (do_ex.eq.0.or.
     &            (chr(i,Ny-2,3).ne.ex.and.chr(i,Ny,3).ne.ex.and.
     &             chr(i,Ny-2,1).ne.ex.and.chr(i,Ny,1).ne.ex))
     &           f(i,Ny-1,2)=0.25*(f(i,Ny-2,3)+f(i,Ny,3)+
     &                             f(i,Ny-2,1)+f(i,Ny,1))
              end do
           end if
           if (Nz.gt.1.and.(phys_bdy(6).eq.pf.or.all)) then
              do i=3,Nx-2
                 if (do_ex.eq.0.or.
     &            (chr(i,Ny-2,Nz-2).ne.ex.and.chr(i,Ny,Nz-2).ne.ex.and.
     &             chr(i,Ny-2,Nz).ne.ex.and.chr(i,Ny,Nz).ne.ex))
     &           f(i,Ny-1,Nz-1)=0.25*(f(i,Ny-2,Nz-2)+f(i,Ny,Nz-2)+
     &                                f(i,Ny-2,Nz)+f(i,Ny,Nz))
              end do
           end if
        end if

        if (Nz.gt.1.and.(phys_bdy(5).eq.pf.or.all)) then
           do i=3,Nx-2
              do j=3,Ny-2
                 if (do_ex.eq.0.or.
     &               (chr(i,j,3).ne.ex.and.chr(i,j,1).ne.ex))
     &           f(i,j,2)=0.5*(f(i,j,3)+f(i,j,1))
              end do
           end do
        end if

        if (Nz.gt.1.and.(phys_bdy(6).eq.pf.or.all)) then
           do i=3,Nx-2
              do j=3,Ny-2
                 if (do_ex.eq.0.or.
     &               (chr(i,j,Nz-2).ne.ex.and.chr(i,j,Nz).ne.ex))
     &           f(i,j,Nz-1)=0.5*(f(i,j,Nz-2)+f(i,j,Nz))
              end do
           end do
        end if

        return
        end

c---------------------------------------------------------------------------
c repopulates all grid points of f *one* point in from an excision surface,
c via at most (linear[2]/quadratic[3]/cubic[4]) extrapolation from adjacent 
c cells 
c
c   * --- o --- o --- o --- o 
c
c   f_ex  f1    f2    f3    f4
c
c 1st order: f_ex = f1
c
c 2nd order: f_ex = 2*f1 - f2
c
c 3rd order: f_ex = 3*f1 - 3*f2 + f3
c
c 4th order: f_ex = 4*f1 - 6*f2 + 4*f3 - f4
c
c io = interpolation order
c---------------------------------------------------------------------------
        subroutine dmrepop3d1(f,chr,ex,io,nx,ny,nz)
        implicit none
        integer nx,ny,nz,io
        real*8 f(nx,ny,nz)
        real*8 chr(nx,ny,nz),ex

        integer i,j,k,wx,wy,wz
        real*8 fx,fy,fz

        if (io.lt.1.or.io.gt.4) then
           write(*,*) 'dmrepop3d1: only straight-copy(1), '
           write(*,*) 'linear(2), quadratic(3)'
           write(*,*) 'and cubic(4) supported at this time'
           stop
        end if

        do i=1,Nx
           do j=1,Ny
              do k=1,Nz
                 if (chr(i,j,k).eq.ex) then
                    wx=0
                    wy=0
                    wz=0
                    fx=0
                    fy=0
                    fz=0

!LR: equivalent to the commented out piece below, but it does not go out of bounds
                    if ((i.eq.2.or.(i.ge.2.and.io.eq.1))) then
                     if (chr(i-1,j,k).ne.ex) then
                       fx=f(i-1,j,k)
                       wx=1
                     end if
                    else if ((i.eq.3.or.(i.ge.3.and.io.eq.2))) then
                     if (chr(i-1,j,k).ne.ex.and.chr(i-2,j,k).ne.ex) then
                       fx=2*f(i-1,j,k)-f(i-2,j,k)
                       wx=2
                     end if
                    else if ((i.eq.4.or.(i.ge.4.and.io.eq.3))) then
                     if (chr(i-1,j,k).ne.ex.and.chr(i-2,j,k).ne.ex.and.
     &                  chr(i-3,j,k).ne.ex) then
                       fx=3*(f(i-1,j,k)-f(i-2,j,k))+f(i-3,j,k)
                       wx=3
                     end if
                    else if ((i.ge.5.and.io.eq.4)) then
                     if (chr(i-1,j,k).ne.ex.and.chr(i-2,j,k).ne.ex.and.
     &                  chr(i-3,j,k).ne.ex.and.chr(i-4,j,k).ne.ex) then
                       fx= 4*(f(i-1,j,k)+f(i-3,j,k))
     &                    -6*f(i-2,j,k)-f(i-4,j,k)
                       wx=4
                     end if
                    else if ((i.eq.Nx-1.or.
     &                       (i.le.Nx-1.and.io.eq.1))) then
                     if (chr(i+1,j,k).ne.ex) then
                       fx=f(i+1,j,k)
                       wx=1
                     end if
                    else if ((i.eq.Nx-2.or.
     &                       (i.le.Nx-2.and.io.eq.2))) then
                     if (chr(i+1,j,k).ne.ex.and.chr(i+2,j,k).ne.ex) then
                       fx=2*f(i+1,j,k)-f(i+2,j,k)
                       wx=2
                     end if
                    else if ((i.eq.Nx-3.or.
     &                       (i.le.Nx-3.and.io.eq.3))) then
                     if (chr(i+1,j,k).ne.ex.and.chr(i+2,j,k).ne.ex.and.
     &                  chr(i+3,j,k).ne.ex) then
                       fx=3*(f(i+1,j,k)-f(i+2,j,k))+f(i+3,j,k)
                       wx=3
                     end if
                    else if ((i.le.Nx-4.and.io.eq.4)) then
                     if (chr(i+1,j,k).ne.ex.and.chr(i+2,j,k).ne.ex.and.
     &                  chr(i+3,j,k).ne.ex.and.chr(i+4,j,k).ne.ex) then
                       fx= 4*(f(i+1,j,k)+f(i+3,j,k))
     &                    -6*f(i+2,j,k)-f(i+4,j,k)
                       wx=4
                     end if
                    end if
!LR commented this out
!                    if ((i.eq.2.or.(i.ge.2.and.io.eq.1)).and.
!     &                  chr(i-1,j,k).ne.ex) then
!                       fx=f(i-1,j,k)
!                       wx=1
!                    else if ((i.eq.3.or.(i.ge.3.and.io.eq.2)).and.
!     &                  chr(i-1,j,k).ne.ex.and.chr(i-2,j,k).ne.ex) then
!                       fx=2*f(i-1,j,k)-f(i-2,j,k)
!                       wx=2
!                    else if ((i.eq.4.or.(i.ge.4.and.io.eq.3)).and.
!     &                  chr(i-1,j,k).ne.ex.and.chr(i-2,j,k).ne.ex.and.
!     &                  chr(i-3,j,k).ne.ex) then
!                       fx=3*(f(i-1,j,k)-f(i-2,j,k))+f(i-3,j,k)
!                       wx=3
!                    else if ((i.ge.5.and.io.eq.4).and.
!     &                  chr(i-1,j,k).ne.ex.and.chr(i-2,j,k).ne.ex.and.
!     &                  chr(i-3,j,k).ne.ex.and.chr(i-4,j,k).ne.ex) then
!                       fx= 4*(f(i-1,j,k)+f(i-3,j,k))
!     &                    -6*f(i-2,j,k)-f(i-4,j,k)
!                       wx=4
!                    else if ((i.eq.Nx-1.or.(i.le.Nx-1.and.io.eq.1)).and.
!     &                  chr(i+1,j,k).ne.ex) then
!                       fx=f(i+1,j,k)
!                       wx=1
!                    else if ((i.eq.Nx-2.or.(i.le.Nx-2.and.io.eq.2)).and.
!     &                  chr(i+1,j,k).ne.ex.and.chr(i+2,j,k).ne.ex) then
!                       fx=2*f(i+1,j,k)-f(i+2,j,k)
!                       wx=2
!                    else if ((i.eq.Nx-3.or.(i.le.Nx-3.and.io.eq.3)).and.
!     &                  chr(i+1,j,k).ne.ex.and.chr(i+2,j,k).ne.ex.and.
!     &                  chr(i+3,j,k).ne.ex) then
!                       fx=3*(f(i+1,j,k)-f(i+2,j,k))+f(i+3,j,k)
!                       wx=3
!                    else if ((i.le.Nx-4.and.io.eq.4).and.
!     &                  chr(i+1,j,k).ne.ex.and.chr(i+2,j,k).ne.ex.and.
!     &                  chr(i+3,j,k).ne.ex.and.chr(i+4,j,k).ne.ex) then
!                       fx= 4*(f(i+1,j,k)+f(i+3,j,k))
!     &                    -6*f(i+2,j,k)-f(i+4,j,k)
!                       wx=4
!                    end if

!LR: equivalent to the commented out piece below, but it does not go out of bounds
                    if ((j.eq.2.or.(j.ge.2.and.io.eq.1))) then
                     if (chr(i,j-1,k).ne.ex) then
                       fy=f(i,j-1,k)
                       wy=1
                     end if
                    else if ((j.eq.3.or.(j.ge.3.and.io.eq.2))) then
                     if (chr(i,j-1,k).ne.ex.and.chr(i,j-2,k).ne.ex) then
                       fy=2*f(i,j-1,k)-f(i,j-2,k)
                       wy=2
                     end if
                    else if ((j.eq.4.or.(j.ge.4.and.io.eq.3))) then
                     if (chr(i,j-1,k).ne.ex.and.chr(i,j-2,k).ne.ex.and.
     &                  chr(i,j-3,k).ne.ex) then
                       fy=3*(f(i,j-1,k)-f(i,j-2,k))+f(i,j-3,k)
                       wy=3
                     end if
                    else if ((j.ge.5.and.io.eq.4)) then
                     if (chr(i,j-1,k).ne.ex.and.chr(i,j-2,k).ne.ex.and.
     &                  chr(i,j-3,k).ne.ex.and.chr(i,j-4,k).ne.ex) then
                       fy= 4*(f(i,j-1,k)+f(i,j-3,k))
     &                    -6*f(i,j-2,k)-f(i,j-4,k)
                       wy=4
                     end if
                    else if ((j.eq.Ny-1.or.
     &                       (j.le.Ny-1.and.io.eq.1))) then
                     if (chr(i,j+1,k).ne.ex) then
                       fy=f(i,j+1,k)
                       wy=1
                     end if
                    else if ((j.eq.Ny-2.or.
     &                       (j.le.Ny-2.and.io.eq.2))) then
                     if (chr(i,j+1,k).ne.ex.and.chr(i,j+2,k).ne.ex) then
                       fy=2*f(i,j+1,k)-f(i,j+2,k)
                       wy=2
                     end if
                    else if ((j.eq.Ny-3.or.
     &                       (j.le.Ny-3.and.io.eq.3))) then
                     if (chr(i,j+1,k).ne.ex.and.chr(i,j+2,k).ne.ex.and.
     &                  chr(i,j+3,k).ne.ex) then
                       fy=3*(f(i,j+1,k)-f(i,j+2,k))+f(i,j+3,k)
                       wy=3
                     end if
                    else if ((j.le.Ny-4.and.io.eq.4)) then
                     if (chr(i,j+1,k).ne.ex.and.chr(i,j+2,k).ne.ex.and.
     &                  chr(i,j+3,k).ne.ex.and.chr(i,j+4,k).ne.ex) then
                       fy= 4*(f(i,j+1,k)+f(i,j+3,k))
     &                    -6*f(i,j+2,k)-f(i,j+4,k)
                       wy=4
                     end if
                    end if
!LR commented this out
!                    if ((j.eq.2.or.(j.ge.2.and.io.eq.1)).and.
!     &                  chr(i,j-1,k).ne.ex) then
!                       fy=f(i,j-1,k)
!                       wy=1
!                    else if ((j.eq.3.or.(j.ge.3.and.io.eq.2)).and.
!     &                  chr(i,j-1,k).ne.ex.and.chr(i,j-2,k).ne.ex) then
!                       fy=2*f(i,j-1,k)-f(i,j-2,k)
!                       wy=2
!                    else if ((j.eq.4.or.(j.ge.4.and.io.eq.3)).and.
!     &                  chr(i,j-1,k).ne.ex.and.chr(i,j-2,k).ne.ex.and.
!     &                  chr(i,j-3,k).ne.ex) then
!                       fy=3*(f(i,j-1,k)-f(i,j-2,k))+f(i,j-3,k)
!                       wy=3
!                    else if ((j.ge.5.and.io.eq.4).and.
!     &                  chr(i,j-1,k).ne.ex.and.chr(i,j-2,k).ne.ex.and.
!     &                  chr(i,j-3,k).ne.ex.and.chr(i,j-4,k).ne.ex) then
!                       fy= 4*(f(i,j-1,k)+f(i,j-3,k))
!     &                    -6*f(i,j-2,k)-f(i,j-4,k)
!                       wy=4
!                    else if ((j.eq.Ny-1.or.(j.le.Ny-1.and.io.eq.1)).and.
!     &                  chr(i,j+1,k).ne.ex) then
!                       fy=f(i,j+1,k)
!                       wy=1
!                    else if ((j.eq.Ny-2.or.(j.le.Ny-2.and.io.eq.2)).and.
!     &                  chr(i,j+1,k).ne.ex.and.chr(i,j+2,k).ne.ex) then
!                       fy=2*f(i,j+1,k)-f(i,j+2,k)
!                       wy=2
!                    else if ((j.eq.Ny-3.or.(j.le.Ny-3.and.io.eq.3)).and.
!     &                  chr(i,j+1,k).ne.ex.and.chr(i,j+2,k).ne.ex.and.
!     &                  chr(i,j+3,k).ne.ex) then
!                       fy=3*(f(i,j+1,k)-f(i,j+2,k))+f(i,j+3,k)
!                       wy=3
!                    else if ((j.le.Ny-4.and.io.eq.4).and.
!     &                  chr(i,j+1,k).ne.ex.and.chr(i,j+2,k).ne.ex.and.
!     &                  chr(i,j+3,k).ne.ex.and.chr(i,j+4,k).ne.ex) then
!                       fy= 4*(f(i,j+1,k)+f(i,j+3,k))
!     &                    -6*f(i,j+2,k)-f(i,j+4,k)
!                       wy=4
!                    end if

!LR: equivalent to the commented out piece below, but it does not go out of bounds
                    if ((k.eq.2.or.(k.ge.2.and.io.eq.1))) then
                     if (chr(i,j,k-1).ne.ex) then
                       fz=f(i,j,k-1)
                       wz=1
                     end if
                    else if ((k.eq.3.or.(k.ge.3.and.io.eq.2))) then
                     if (chr(i,j,k-1).ne.ex.and.chr(i,j,k-2).ne.ex) then
                       fz=2*f(i,j,k-1)-f(i,j,k-2)
                       wz=2
                     end if
                    else if ((k.eq.4.or.(k.ge.4.and.io.eq.3))) then
                     if (chr(i,j,k-1).ne.ex.and.chr(i,j,k-2).ne.ex.and.
     &                  chr(i,j,k-3).ne.ex) then
                       fz=3*(f(i,j,k-1)-f(i,j,k-2))+f(i,j,k-3)
                       wz=3
                     end if
                    else if ((k.ge.5.and.io.eq.4)) then
                     if (chr(i,j,k-1).ne.ex.and.chr(i,j,k-2).ne.ex.and.
     &                  chr(i,j,k-3).ne.ex.and.chr(i,j,k-4).ne.ex) then
                       fz= 4*(f(i,j,k-1)+f(i,j,k-3))
     &                    -6*f(i,j,k-2)-f(i,j,k-4)
                       wz=4
                     end if
                    else if ((k.eq.Nz-1.or.
     &                       (k.le.Nz-1.and.io.eq.1))) then
                     if (chr(i,j,k+1).ne.ex) then
                       fz=f(i,j,k+1)
                       wz=1
                     end if
                    else if ((k.eq.Nz-2.or.
     &                       (k.le.Nz-2.and.io.eq.2))) then
                     if (chr(i,j,k+1).ne.ex.and.chr(i,j,k+2).ne.ex) then
                       fz=2*f(i,j,k+1)-f(i,j,k+2)
                       wz=2
                     end if
                    else if ((k.eq.Nz-3.or.
     &                       (k.le.Nz-3.and.io.eq.3))) then
                     if (chr(i,j,k+1).ne.ex.and.chr(i,j,k+2).ne.ex.and.
     &                  chr(i,j,k+3).ne.ex) then
                       fz=3*(f(i,j,k+1)-f(i,j,k+2))+f(i,j,k+3)
                       wz=3
                     end if
                    else if ((k.le.Nz-4.and.io.eq.4)) then
                     if (chr(i,j,k+1).ne.ex.and.chr(i,j,k+2).ne.ex.and.
     &                  chr(i,j,k+3).ne.ex.and.chr(i,j,k+4).ne.ex) then
                       fz= 4*(f(i,j,k+1)+f(i,j,k+3))
     &                    -6*f(i,j,k+2)-f(i,j,k+4)
                       wz=4
                     end if
                    end if

!LR commented this out                    
!                    if ((k.eq.2.or.(k.ge.2.and.io.eq.1)).and.
!     &                  chr(i,j,k-1).ne.ex) then
!                       fz=f(i,j,k-1)
!                       wz=1
!                    else if ((k.eq.3.or.(k.ge.3.and.io.eq.2)).and.
!     &                  chr(i,j,k-1).ne.ex.and.chr(i,j,k-2).ne.ex) then
!                       fz=2*f(i,j,k-1)-f(i,j,k-2)
!                       wz=2
!                    else if ((k.eq.4.or.(k.ge.4.and.io.eq.3)).and.
!     &                  chr(i,j,k-1).ne.ex.and.chr(i,j,k-2).ne.ex.and.
!     &                  chr(i,j,k-3).ne.ex) then
!                       fz=3*(f(i,j,k-1)-f(i,j,k-2))+f(i,j,k-3)
!                       wz=3
!                    else if ((k.ge.5.and.io.eq.4).and.
!     &                  chr(i,j,k-1).ne.ex.and.chr(i,j,k-2).ne.ex.and.
!     &                  chr(i,j,k-3).ne.ex.and.chr(i,j,k-4).ne.ex) then
!                       fz= 4*(f(i,j,k-1)+f(i,j,k-3))
!     &                    -6*f(i,j,k-2)-f(i,j,k-4)
!                       wz=4
!                    else if ((k.eq.Nz-1.or.(k.le.Nz-1.and.io.eq.1)).and.
!     &                  chr(i,j,k+1).ne.ex) then
!                       fz=f(i,j,k+1)
!                       wz=1
!                    else if ((k.eq.Nz-2.or.(k.le.Nz-2.and.io.eq.2)).and.
!     &                  chr(i,j,k+1).ne.ex.and.chr(i,j,k+2).ne.ex) then
!                       fz=2*f(i,j,k+1)-f(i,j,k+2)
!                       wz=2
!                    else if ((k.eq.Nz-3.or.(k.le.Nz-3.and.io.eq.3)).and.
!     &                  chr(i,j,k+1).ne.ex.and.chr(i,j,k+2).ne.ex.and.
!     &                  chr(i,j,k+3).ne.ex) then
!                       fz=3*(f(i,j,k+1)-f(i,j,k+2))+f(i,j,k+3)
!                       wz=3
!                    else if ((k.le.Nz-4.and.io.eq.4).and.
!     &                  chr(i,j,k+1).ne.ex.and.chr(i,j,k+2).ne.ex.and.
!     &                  chr(i,j,k+3).ne.ex.and.chr(i,j,k+4).ne.ex) then
!                       fz= 4*(f(i,j,k+1)+f(i,j,k+3))
!     &                    -6*f(i,j,k+2)-f(i,j,k+4)
!                       wz=4
!                    end if

                    if (wx+wy+wz.ne.0) then
                       !-----------------------------------------------
                       ! only blend extrapolations of the same order
                       !-----------------------------------------------
                       if (wx.gt.wy) wy=0
                       if (wx.gt.wz) wz=0
                       if (wy.gt.wx) wx=0
                       if (wy.gt.wz) wz=0
                       if (wz.gt.wx) wx=0
                       if (wz.gt.wy) wy=0
                       f(i,j,k)=(wx*fx+wy*fy+wz*fz)/(wx+wy+wz)
                    end if
                 end if
              end do
           end do
        end do
                 
        return
        end

c---------------------------------------------------------------------------
c repopulates the characteristic array *one* point in from an excision 
c surface ... NOTE : changes values of unexcised points
c---------------------------------------------------------------------------
        subroutine dmrepop3dc1(chr,ex,nx,ny,nz)
        implicit none
        integer nx,ny,nz
        real*8 chr(nx,ny,nz),ex

        integer i,j,k
        real*8 un_ex0,un_ex1

        !-------------------------------------------------------------------
        ! we want a definite value for 'original' unexcised pieces ... 
        ! set to ex-1
        !-------------------------------------------------------------------
        un_ex0=ex-1
        un_ex1=ex-2
        do i=1,Nx
           do j=1,Ny
              do k=1,Nz
                 if (chr(i,j,k).ne.ex) chr(i,j,k)=un_ex0
              end do
           end do
        end do

        do i=1,Nx
           do j=1,Ny
              do k=1,Nz

!LR: equivalent to the commented out piece below, but it does not go out of bounds
                 if (chr(i,j,k).eq.ex) then
                  if (i.ge.2) then
                   if (chr(i-1,j,k).eq.un_ex0) then
                     chr(i,j,k)=un_ex1
                   end if
                  else if (i.le.Nx-1) then
                   if (chr(i+1,j,k).eq.un_ex0) then
                     chr(i,j,k)=un_ex1
                   end if
                  else if (j.ge.2) then
                   if (chr(i,j-1,k).eq.un_ex0) then
                     chr(i,j,k)=un_ex1
                   end if
                  else if (j.le.Ny-1) then
                   if (chr(i,j+1,k).eq.un_ex0) then
                     chr(i,j,k)=un_ex1
                   end if
                  else if (k.ge.2) then
                   if (chr(i,j,k-1).eq.un_ex0) then
                     chr(i,j,k)=un_ex1
                   end if
                  else if (k.le.Nz-1) then
                   if (chr(i,j,k+1).eq.un_ex0) then
                     chr(i,j,k)=un_ex1
                   end if
                  end if
                 end if
!LR commented this out
!                 if (chr(i,j,k).eq.ex.and.
!     &               ((i.ge.2.and.chr(i-1,j,k).eq.un_ex0).or.
!     &                (i.le.Nx-1.and.chr(i+1,j,k).eq.un_ex0).or.
!     &                (j.ge.2.and.chr(i,j-1,k).eq.un_ex0).or.
!     &                (j.le.Ny-1.and.chr(i,j+1,k).eq.un_ex0).or.
!     &                (k.ge.2.and.chr(i,j,k-1).eq.un_ex0).or.
!     &                (k.le.Nz-1.and.chr(i,j,k+1).eq.un_ex0))) then
!                     chr(i,j,k)=un_ex1
!                 end if

              end do
           end do
        end do

        do i=1,Nx
           do j=1,Ny
              do k=1,Nz
                 if (chr(i,j,k).ne.ex) chr(i,j,k)=un_ex0
              end do
           end do
        end do
                 
        return
        end

c---------------------------------------------------------------------------
c The equivalent of m_util_r8.f for cell centered grids
c---------------------------------------------------------------------------


c---------------------------------------------------------------------------
c dminterp3d_c interpolates a segment of grid u1 into grid u2
c
c offset in u1 = [i1,j1,k1] with stride [si1,sj1,sk1] (nominally [1,1,1])
c offset in u2 = [i2,j2,k2] with rho    [irho,jrho,krho] 
c size of segment=ni x nj x nk in u1 of the coarse
c si1, sj1, sk1 is the jump (strive)
c
c ord is type of interpolation as enumerated in ops.inc
c
c if do_ex is non-zero, then chr1() specifies an excision mask
c on the coarse level : i.e. if chr1(i,j,k)=ex then that
c cell is excised, and should not be used in the interpolation
c
c For MC interpolation, copy boundary conditions will be applied
c to the stencil for excised zones.  This has the effect of clipping
c to zero slope.
c
c KNOWN ISSUE with excision:  Introducing a fine grid atop a grid
c containing the excision zone can result in loss of rest-mass, etc.
c This is because the mass is doled out to all child cells equally.
c If one or more of those child cells is excised, then that mass
c will be lost and unaccounted for.  This is unlikely to be a very
c large effect.  It will only affect those grid functions which 
c use FIRST_ORDER_EXTENSIVE interpolation.
c
c dim is the dimension of the grids. A combination of qrho, nq1 and nq2 
c could be used to figure out whether q is not present (i.e. the
c dimension is less than q), q is present but not refined,
c or a slice of thickness 1 in q on the coarse cell is to be refined.
c However, passing the dimension is cleaner (these ambiguities
c don't arise in the vertex centered routine)
c
c This routine assumes that the destination grid is aligned with 
c the source grid.  This should always be possible since u2 here
c is a temporary data array.  
c---------------------------------------------------------------------------
      subroutine dminterp3d_c(u1,u2,nx1,ny1,nz1,nx2,ny2,nz2,
     &     i1,j1,k1,i2,j2,k2,ni,nj,nk,
     &     si1,sj1,sk1,irho,jrho,krho,ord,
     &     chr1,ex,do_ex,dim)
      implicit none
      integer nx1,ny1,nz1,nx2,ny2,nz2,dim
      integer i1,j1,k1,i2,j2,k2,ni,nj,nk
      integer si1,sj1,sk1,irho,jrho,krho,ord,do_ex
      real*8 u1(nx1,ny1,nz1),u2(nx2,ny2,nz2)
      real*8 chr1(nx1,ny1,nz1),ex
      real*8 mc_slope
      
      logical       ltrace
      parameter   ( ltrace = .false. )
      
      character*12  cdnm
      parameter   ( cdnm = 'dminterp3d_c' )
      
      
c-----------------------------------------------------------------------
c     Local variables
c-----------------------------------------------------------------------
      
      integer i,j,k,l,m,n
      integer ic,  jc,  kc 
      integer nni, nnj, nnk
      real*8 dux, duy, duz
      real*8 v1,v2,v3
      include 'ops.inc'
      
c=======================================================================
c=======================================================================
c=======================================================================
c=======================================================================

      if ( irho .ne. 2 .or. jrho .ne. 2 .or. krho .ne. 2 ) then
         write(0,*) cdnm,': irho = ', irho
         write(0,*) cdnm,': jrho = ', jrho
         write(0,*) cdnm,': krho = ', krho
         write(0,*) cdnm,': Only 2:1 refinement allowed' 
      endif
      
      if ( si1 .ne. 1 .or. sj1 .ne. 1 .or. sk1 .ne. 1 ) then
         write(0,*) cdnm,': si1 = ', si1
         write(0,*) cdnm,': sj1 = ', sj1
         write(0,*) cdnm,': sk1 = ', sk1
      endif
      if ( ltrace ) then 
         write(0,*) cdnm,': i1  = ', i1
         write(0,*) cdnm,': j1  = ', j1
         write(0,*) cdnm,': k1  = ', k1
         write(0,*) cdnm,': i2  = ', i2
         write(0,*) cdnm,': j2  = ', j2
         write(0,*) cdnm,': k2  = ', k2
         write(0,*) cdnm,': ni  = ', ni
         write(0,*) cdnm,': nj  = ', nj
         write(0,*) cdnm,': nk  = ', nk
         write(0,*) cdnm,': dim = ', dim
      endif

      if ( 2*ni+i2-1 .le. 2*nx2-1 ) then
         nni = 2*ni+i2-1
      else
         nni = 2*nx2-1
      endif 
      if ( 2*nj+j2-1 .le. 2*ny2-1 ) then
         nnj = 2*nj+j2-1
      else
         nnj = 2*ny2-1
      endif 
      if ( 2*nj+j2-1 .le. 2*ny2-1 ) then
         nnj = 2*nj+j2-1
      else
         nnj = 2*ny2-1
      endif 
      if ( 2*nk+k2-1 .le. 2*nz2-1 ) then
         nnk = 2*nk+k2-1
      else
         nnk = 2*nz2-1
      endif 
      if ( ord .eq. PAMR_FIRST_ORDER_CONS ) then 
         if ( dim .eq. 1 ) then
            do i = i2, nni, 2
               ic = (i-i2)/2 + i1 
               if (do_ex.eq.0 .or. chr1(ic,1,1).ne.ex) then
                  u2(i,   1,   1  ) = u1(ic, 1, 1)           
                  u2(i+1, 1,   1  ) = u1(ic, 1, 1)           
               else
                  u2(i,   1,   1  ) = 0.d0
                  u2(i+1, 1,   1  ) = 0.d0
               end if
            enddo
         elseif ( dim .eq. 2 ) then
            do i = i2, nni, 2
               do j = j2, nnj, 2
                  ic = (i-i2)/2 + i1
                  jc = (j-j2)/2 + j1
                  if (do_ex.eq.0 .or. chr1(ic,jc,1).ne.ex) then
                     u2(i,   j,   1  ) = u1(ic, jc, 1)           
                     u2(i+1, j,   1  ) = u1(ic, jc, 1)           
                     u2(i,   j+1, 1  ) = u1(ic, jc, 1)           
                     u2(i+1, j+1, 1  ) = u1(ic, jc, 1)           
                  else
                     u2(i,   j,   1  ) = 0.d0
                     u2(i+1, j,   1  ) = 0.d0
                     u2(i,   j+1, 1  ) = 0.d0
                     u2(i+1, j+1, 1  ) = 0.d0
                  end if
               enddo
            enddo
         elseif ( dim .eq. 3 ) then
c-----------------------------------------------------------------------
c     Note that the assumption of 2:1 refinment is heavily used 
c-----------------------------------------------------------------------
            do i = i2, nni, 2
               do j = j2, nnj, 2
                  do k = k2, nnk, 2
c-----------------------------------------------------------------------
c     Calculation of the coarse cell index
c-----------------------------------------------------------------------
                     ic = (i-i2)/2 + i1
                     jc = (j-j2)/2 + j1
                     kc = (k-k2)/2 + k1
c-----------------------------------------------------------------------
c     The eight new cells are interpolated with the coarse cell value.
c-----------------------------------------------------------------------
                     if (do_ex.eq.0 .or. chr1(ic,jc,kc).ne.ex) then
                        u2(i,   j,   k  ) = u1(ic, jc, kc)           
                        u2(i+1, j,   k  ) = u1(ic, jc, kc)           
                        u2(i,   j+1, k  ) = u1(ic, jc, kc)           
                        u2(i+1, j+1, k  ) = u1(ic, jc, kc)           
                        u2(i,   j,   k+1) = u1(ic, jc, kc)           
                        u2(i+1, j,   k+1) = u1(ic, jc, kc)           
                        u2(i,   j+1, k+1) = u1(ic, jc, kc)           
                        u2(i+1, j+1, k+1) = u1(ic, jc, kc)           
                     else
                        u2(i,   j,   k  ) = 0.d0
                        u2(i+1, j,   k  ) = 0.d0
                        u2(i,   j+1, k  ) = 0.d0
                        u2(i+1, j+1, k  ) = 0.d0
                        u2(i,   j,   k+1) = 0.d0
                        u2(i+1, j,   k+1) = 0.d0
                        u2(i,   j+1, k+1) = 0.d0
                        u2(i+1, j+1, k+1) = 0.d0
                     end if
                  enddo
               enddo
            enddo
         endif
         
      elseif ( ord .eq. PAMR_FIRST_ORDER_EXTENSIVE ) then
         if ( dim .eq. 1 ) then
            do i = i2, nni, 2
               ic = (i-i2)/2 + i1 
               if (do_ex.eq.0 .or. chr1(ic,1,1).ne.ex) then
                  u2(i,   1,   1  ) = u1(ic, 1, 1)/2.d0           
                  u2(i+1, 1,   1  ) = u1(ic, 1, 1)/2.d0          
               else            
                  u2(i,   1,   1  ) = 0.d0
                  u2(i+1, 1,   1  ) = 0.d0
               end if
            enddo
         elseif ( dim .eq. 2 ) then
            do i = i2, nni, 2
               do j = j2, nnj, 2
                  ic = (i-i2)/2 + i1
                  jc = (j-j2)/2 + j1
                  if (do_ex.eq.0 .or. chr1(ic,jc,1).ne.ex) then
                     u2(i,   j,   1  ) = u1(ic, jc, 1)/4.d0           
                     u2(i+1, j,   1  ) = u1(ic, jc, 1)/4.d0                      
                     u2(i,   j+1, 1  ) = u1(ic, jc, 1)/4.d0                      
                     u2(i+1, j+1, 1  ) = u1(ic, jc, 1)/4.d0                
                 else
                     u2(i,   j,   1  ) = 0.d0
                     u2(i+1, j,   1  ) = 0.d0
                     u2(i,   j+1, 1  ) = 0.d0
                     u2(i+1, j+1, 1  ) = 0.d0
                  end if      
               enddo
            enddo
         elseif ( dim .eq. 3 ) then
c-----------------------------------------------------------------------
c     Note that the assumption of 2:1 refinment is heavily used 
c-----------------------------------------------------------------------
            do i = i2, nni, 2
               do j = j2, nnj, 2
                  do k = k2, nnk, 2
c-----------------------------------------------------------------------
c     Calculation of the coarse cell index
c-----------------------------------------------------------------------
                     ic = (i-i2)/2 + i1
                     jc = (j-j2)/2 + j1
                     kc = (k-k2)/2 + k1
c-----------------------------------------------------------------------
c     The eight new cells are interpolated with the coarse cell value over 8
c-----------------------------------------------------------------------
                     if (do_ex.eq.0 .or. chr1(ic,jc,kc).ne.ex) then
                        u2(i,   j,   k  ) = u1(ic, jc, kc)/8.d0           
                        u2(i+1, j,   k  ) = u1(ic, jc, kc)/8.d0           
                        u2(i,   j+1, k  ) = u1(ic, jc, kc)/8.d0                      
                        u2(i+1, j+1, k  ) = u1(ic, jc, kc)/8.d0                      
                        u2(i,   j,   k+1) = u1(ic, jc, kc)/8.d0                      
                        u2(i+1, j,   k+1) = u1(ic, jc, kc)/8.d0                      
                        u2(i,   j+1, k+1) = u1(ic, jc, kc)/8.d0
                        u2(i+1, j+1, k+1) = u1(ic, jc, kc)/8.d0
                     else
                        u2(i,   j,   k  ) = 0.d0
                        u2(i+1, j,   k  ) = 0.d0
                        u2(i,   j+1, k  ) = 0.d0
                        u2(i+1, j+1, k  ) = 0.d0
                        u2(i,   j,   k+1) = 0.d0
                        u2(i+1, j,   k+1) = 0.d0
                        u2(i,   j+1, k+1) = 0.d0
                        u2(i+1, j+1, k+1) = 0.d0
                     end if
                  enddo
               enddo
            enddo
         endif
         
      elseif ( ord .eq. PAMR_MC ) then
         
         if ( dim .eq. 1 ) then

            do i = i2, nni, 2

               ic = (i-i2)/2 + i1 
               ! calculate MC derivative
               ! use one-sided derivative at edges.
               if (i.eq.i2) then
                  if (do_ex.eq.0 .or. chr1(ic+1,1,1).ne.ex) then
                     dux = (u1(ic+1,1,1) - u1(ic,1,1))/2.d0
                  else
                     dux = 0.d0
                  end if
               else if (i.eq.nni .or. i.eq.nni-1) then
                  if (do_ex.eq.0 .or. chr1(ic-1,1,1).ne.ex) then
                     dux = (u1(ic,1,1) - u1(ic-1,1,1))/2.d0
                  else
                     dux = 0.d0
                  end if
               else
                  if (do_ex.eq.0 .or. (chr1(ic-1,1,1).ne.ex .and. 
     &                 chr1(ic+1,1,1).ne.ex) ) then
                     dux = 
     &                    mc_slope(u1(ic-1,1,1),u1(ic,1,1),u1(ic+1,1,1))
                  else
                     dux = 0.d0
                  end if
               end if

               do l=0,1
                  if (do_ex.eq.0 .or. chr1(ic,1,1).ne.ex) then
                     u2(i+l,1,1) = u1(ic,1,1) + (0.5d0*l - 0.25d0) * dux
                  else
                     u2(i+l,1,1) = 0.d0
                  end if
               end do

            enddo
          
         elseif ( dim .eq. 2 ) then

            do j = j2, nnj, 2
               do i = i2, nni, 2

                  ic = (i-i2)/2 + i1 
                  jc = (j-j2)/2 + j1 
                  ! calculate MC derivative
                  ! use one-sided derivative at edges.
                  if (i.eq.i2) then
                     if (do_ex.eq.0 .or. chr1(ic+1,jc,1).ne.ex) then
                        dux = (u1(ic+1,jc,1) - u1(ic,jc,1))/2.d0
                     else
                        dux = 0.d0
                     end if
                  else if (i.eq.nni .or. i.eq.nni-1) then
                     if (do_ex.eq.0 .or. chr1(ic-1,jc,1).ne.ex) then
                        dux = (u1(ic,jc,1) - u1(ic-1,jc,1))/2.d0
                     else
                        dux = 0.d0
                     end if
                  else
                     if (do_ex.eq.0 .or. (chr1(ic-1,jc,1).ne.ex .and.
     &                    chr1(ic+1,jc,1).ne.ex) ) then
                        dux = mc_slope(u1(ic-1,jc,1),u1(ic,jc,1),
     &                       u1(ic+1,jc,1))
                     else
                        dux = 0.d0
                     end if
                  end if

                  if (j.eq.j2) then
                     if (do_ex.eq.0 .or. chr1(ic,jc+1,1).ne.ex) then
                        duy = (u1(ic,jc+1,1) - u1(ic,jc,1))/2.d0
                     else
                        duy = 0.d0
                     end if
                  else if (j.eq.nnj .or. j.eq.nnj-1) then
                     if (do_ex.eq.0 .or. chr1(ic,jc-1,1).ne.ex) then
                        duy = (u1(ic,jc,1) - u1(ic,jc-1,1))/2.d0
                     else
                        duy = 0.d0
                     end if
                  else
                     if (do_ex.eq.0 .or. (chr1(ic,jc-1,1).ne.ex .and.
     &                    chr1(ic,jc+1,1).ne.ex) ) then
                        duy = mc_slope(u1(ic,jc-1,1),u1(ic,jc,1),
     &                       u1(ic,jc+1,1))
                     else
                        duy = 0.d0
                     end if
                  end if
                  
                  do m=0,1
                     do l=0,1
                        if (do_ex.eq.0 .or. chr1(ic,jc,1).ne.ex) then
                           u2(i+l,j+m,1) = u1(ic,jc,1) + 
     &                          (0.5d0*l - 0.25d0) * dux +
     &                          (0.5d0*m - 0.25d0) * duy
                        else
                           u2(i+l,j+m,1) = 0.d0
                        end if
                     end do
                  end do
                  
               end do
            end do

         elseif ( dim .eq. 3 ) then

            do k = k2, nnk, 2
               do j = j2, nnj, 2
                  do i = i2, nni, 2

                     ic = (i-i2)/2 + i1 
                     jc = (j-j2)/2 + j1 
                     kc = (k-k2)/2 + k1 
                     ! calculate MC derivative
                     ! use one-sided derivative at edges.
                     if (i.eq.i2) then
                        if (do_ex.eq.0 .or. chr1(ic+1,jc,kc).ne.ex) then
                           dux = (u1(ic+1,jc,kc) - u1(ic,jc,kc))/2.d0
                        else
                           dux = 0.d0
                        end if
                     else if (i.eq.nni .or. i.eq.nni-1) then
                        if (do_ex.eq.0 .or. chr1(ic-1,jc,kc).ne.ex) then
                           dux = (u1(ic,jc,kc) - u1(ic-1,jc,kc))/2.d0
                        else
                           dux = 0.d0
                        end if
                     else
                        if (do_ex.eq.0 .or. (chr1(ic-1,jc,kc).ne.ex 
     &                       .and. chr1(ic+1,jc,kc).ne.ex) ) then
                           dux = mc_slope(u1(ic-1,jc,kc),u1(ic,jc,kc),
     &                          u1(ic+1,jc,kc))
                        else
                           dux = 0.d0
                        end if
                     end if
                     
                     if (j.eq.j2) then
                        if (do_ex.eq.0 .or. chr1(ic,jc+1,kc).ne.ex) then
                           duy = (u1(ic,jc+1,kc) - u1(ic,jc,kc))/2.d0
                        else
                           duy = 0.d0
                        end if
                     else if (j.eq.nnj .or. j.eq.nnj-1) then
                        if (do_ex.eq.0 .or. chr1(ic,jc-1,kc).ne.ex) then
                           duy = (u1(ic,jc,kc) - u1(ic,jc-1,kc))/2.d0
                        else
                           duy = 0.d0
                        end if
                     else
                        if (do_ex.eq.0 .or. (chr1(ic,jc-1,kc).ne.ex 
     &                       .and. chr1(ic,jc+1,kc).ne.ex) ) then
                           duy = mc_slope(u1(ic,jc-1,kc),u1(ic,jc,kc),
     &                          u1(ic,jc+1,kc))
                        else
                           duy = 0.d0
                        end if
                     end if

                     if (k.eq.k2) then
                        if (do_ex.eq.0 .or. chr1(ic,jc,kc+1).ne.ex) then
                           duz = (u1(ic,jc,kc+1) - u1(ic,jc,kc))/2.d0
                        else
                           duz = 0.d0
                        end if
                     else if (k.eq.nnk .or. k.eq.nnk-1) then
                        if (do_ex.eq.0 .or. chr1(ic,jc,kc-1).ne.ex) then
                           duz = (u1(ic,jc,kc) - u1(ic,jc,kc-1))/2.d0
                        else
                           duz = 0.d0
                        end if
                     else
                        if (do_ex.eq.0 .or. (chr1(ic,jc,kc-1).ne.ex 
     &                       .and. chr1(ic,jc,kc+1).ne.ex) ) then
                           duz = mc_slope(u1(ic,jc,kc-1),u1(ic,jc,kc),
     &                          u1(ic,jc,kc+1))
                        else
                           duz = 0.d0
                        end if
                     end if
                     
                     do n=0,1
                        do m=0,1
                           do l=0,1
                              if (do_ex.eq.0 .or. chr1(ic,jc,kc).ne.ex) 
     &                             then
                                 u2(i+l,j+m,k+n) = u1(ic,jc,kc) + 
     &                                (0.5d0*l - 0.25d0) * dux +
     &                                (0.5d0*m - 0.25d0) * duy +
     &                                (0.5d0*n - 0.25d0) * duz 
                              else
                                 u2(i+l,j+m,k+n) = 0.d0
                              end if
                           end do
                        end do
                     end do
                  
                  end do
               end do
            end do

         end if
         
      elseif ( ord .eq. PAMR_SECOND_ORDER ) then
         write(0,*) cdnm,': NOT implemented yet.'
      elseif ( ord. eq. PAMR_NO_INTERP ) then
         return
      else
         write(0,*) cdnm,': NOT implemented yet.'
      endif
      if ( ltrace ) then
         write(0,*) cdnm,': Leaving Routine'
      endif
      return
      end

c---------------------------------------------------------------------------
c repopulates all cells of f *one* point in from an excision surface,
c via *at most* io'th order extrapolation from adjacent 
c cells (lower order will be used if adjacent cells aren't available)
c---------------------------------------------------------------------------
        subroutine dmrepop3d1_c(f,chr,ex,io,nx,ny,nz)
        implicit none
        integer nx,ny,nz,io
        real*8 f(nx,ny,nz)
        real*8 chr(nx,ny,nz),ex
        integer i,j,k,wx,wy,wz
        real*8 fx,fy,fz

        if (io.lt.1.or.io.gt.2) then
           write(*,*) 'dmrepop3d_c1: only straight-copy(1) '
           write(*,*) 'supported at this time'
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

                    if (i.ge.2 .and. chr(i-1,j,k).ne.ex) then
                       fx=f(i-1,j,k)
                       wx=1
                    else if (i.le.Nx-1 .and. chr(i+1,j,k).ne.ex) then
                       fx=f(i+1,j,k)
                       wx=1
                    end if

                    if (j.ge.2 .and. chr(i,j-1,k).ne.ex) then
                       fy=f(i,j-1,k)
                       wy=1
                    else if (j.le.Ny-1 .and. chr(i,j+1,k).ne.ex) then
                       fy=f(i,j+1,k)
                       wy=1
                    end if
                    
                    if (k.ge.2 .and. chr(i,j,k-1).ne.ex) then
                       fz=f(i,j,k-1)
                       wz=1
                    else if (k.le.Nz-1 .and. chr(i,j,k+1).ne.ex) then
                       fz=f(i,j,k+1)
                       wz=1
                    end if
                    if (wx+wy+wz.ne.0) then
                       !-----------------------------------------------
                       ! For hydro variables, it's best not to try to 
                       ! average the copied values because of the primitive
                       ! variable inversion.  Just choose one.
                       ! I'll give preference to x, then y, then z.  This
                       ! will, of course, introduce an asymmetry near
                       ! the excision zone.  But it is unlikely to make
                       ! much of a difference.  
                       !-----------------------------------------------
                       if (wx .ne. 0) then
                          wy = 0
                          wz = 0
                       else if (wy .ne. 0) then
                          wz = 0
                       end if
                       f(i,j,k)=(wx*fx+wy*fy+wz*fz)/(wx+wy+wz)
                    end if
                 end if
              end do
           end do
        end do

        return
        end

c---------------------------------------------------------------------------
c transfers cell-centered grid f_c(nxc,nyc,...) to vertex centered one 
c f_v(nxv,nyv,...)
c
c NOTE: nxc=nxv-1, etc., *except* for lower dimensional cases, where
c the dummy higher dimensional array sizes should be set to 1.
c
c order is type of interpolation as enumerated in ops.inc 
c So far, only second order is available.  If higher orders
c are ever implemented, we will have to drop to lower order if 
c adjacent to excised points/cells.
c
c If do_ex is not zero, then chr_c & chr_v specify excision masks
c which define cells/points that are excised if the corresponding
c array value is equal to ex.
c---------------------------------------------------------------------------
        subroutine dm_c_to_v(f_c,f_v,chr_c,chr_v,
     &                       nxc,nyc,nzc,nxv,nyv,nzv,
     &                       ord,ex,do_ex,dim)
        implicit none
        integer nxc,nyc,nzc,nxv,nyv,nzv,do_ex,ord,dim
        real*8 f_c(nxc,nyc,nzc),f_v(nxv,nyv,nzv)
        real*8 chr_c(nxc,nyc,nzc),chr_v(nxv,nyv,nzv)
        real*8 ex

        logical     ltrace 
        parameter ( ltrace  = .false. )

        character*9  cdnm
        parameter  ( cdnm = 'dm_c_to_v' )
        include 'ops.inc'

c=======================================================================
c=======================================================================
c=======================================================================
c=======================================================================

        if ( ltrace ) then
           write(0,*) cdnm,': nxc= ', nxc
           write(0,*) cdnm,': nyc= ', nyc
           write(0,*) cdnm,': nzc= ', nzc
           write(0,*) cdnm,': nxv= ', nxv
           write(0,*) cdnm,': nyv= ', nyv
           write(0,*) cdnm,': nzv= ', nzv
           write(0,*) cdnm,': dim= ', dim
        endif
        if (nxc .lt. 3 .or. (dim .gt. 1 .and. nyc .lt. 3) .or. 
     &       (dim .gt. 2 .and. nzc .lt. 3)) then
           write(*,*) 'grid too small for c_to_v!'
           write(*,*) 'nxc = ', nxc, ' nyc = ', nyc, ' nzc = ', nzc
           stop
        end if
        
        if (ord .eq. PAMR_C_TO_V_NO_TRANSFER) then
           return
        else if (.not.(ord.eq.PAMR_C_TO_V_SECOND_ORDER)) then
           write(*,*) 'dm_c_to_v: ord = ', ord, ' not yet implemented'
           stop
        end if

        if (dim .eq. 1) then
           call dm_c_to_v_1d(f_c,f_v,chr_v,nxc,nxv,ex,do_ex)
        else if (dim .eq. 2) then
           call dm_c_to_v_2d(f_c,f_v,chr_v,nxc,nyc,nxv,nyv,ex,do_ex)
        else if (dim .eq. 3) then
           call dm_c_to_v_3d(f_c,f_v,chr_v,nxc,nyc,nzc,nxv,nyv,nzv,
     &          ex,do_ex)
        else
           write(*,*) 'dm_c_to_v: dim = ', dim, ' not supported'
           stop
        end if

        return
        end

c---------------------------------------------------------------------------
c For 3D arrays, performs a 2nd-order interpolation (averaging/extrap.) of 
c cell-centered function (f_c) into a vertex-centered function (f_v). 
c
c Note that, for an un-excised vertex, there will always be eight surrouding
c unexcised cells, since excised cells are defined as being interior to 
c the boundary formed by the first layer of excised vertices.  This is why
c we only check to see if the vertex in question is excised.  
c
c This routine assumes that the excision zone will not come within two cell
c widths of an AMR boundary.  Otherwise, the technique of using diagonals
c to extrapolate to the surface will fail.  
c 
c  -- does not check to see if  ny*,nz* > 1
c---------------------------------------------------------------------------
      subroutine dm_c_to_v_3d(f_c,f_v,chr_v,nxc,nyc,nzc,nxv,nyv,nzv,
     &     ex,do_ex)
        implicit none
        integer nxc,nyc,nzc,nxv,nyv,nzv,ord
        integer do_ex
        real*8 ex
        real*8 f_c(nxc,nyc,nzc),f_v(nxv,nyv,nzv)
        real*8 chr_v(nxv,nyv,nzv)
        integer i,j,k, ntotc, ntotv, startc, startv
        integer is, js, ks, ie, je, ke
        integer i_c1, j_c1, k_c1, i_c2, j_c2, k_c2
        integer i_offset1, j_offset1, k_offset1
        integer i_offset2, j_offset2, k_offset2
        integer i_dir, j_dir, k_dir
        integer ibeg(2), jbeg(2), kbeg(2), iend(2),jend(2),kend(2)
        integer ind_diff

        real*8 half,fourth,eighth
        parameter  ( eighth       =  0.1250 0000 0000 0000 d0,
     &               fourth       =  0.2500 0000 0000 0000 d0,
     &               half         =  0.5000 0000 0000 0000 d0 )

        
        !-------------------------------------------------------------------
        !  Order of interpolation and extrapolation
        !-------------------------------------------------------------------
        ord = 2
        !
        !-------------------------------------------------------------------
        !  Interior Points:Interpolation: (averaging nearest-neighbors)
        !-------------------------------------------------------------------

        do k=2,nzv-1
           do j=2,nyv-1
              do i=2,nxv-1
                 if (do_ex.eq.1 .and. chr_v(i,j,k).eq.ex) then
                    f_v(i,j,k) = 0.d0
                 else
                    f_v(i,j,k) = eighth * ( f_c(i-1, j-1, k-1) +
     &                                      f_c(i  , j-1, k-1) +
     &                                      f_c(i-1, j  , k-1) +
     &                                      f_c(i  , j  , k-1) +
     &                                      f_c(i-1, j-1, k  ) +
     &                                      f_c(i  , j-1, k  ) +
     &                                      f_c(i-1, j  , k  ) +
     &                                      f_c(i  , j  , k  )  )
                 end if
              end do
           end do
        end do

        !-------------------------------------------------------------------
        !  Face points: (average the extrapolations of 4 diagonals)
        !-------------------------------------------------------------------

        ntotv = nxv*nyv*nzv
        ntotc = nxc*nyc*nzc

        !- zero-out faces and corners so that we can just add 
        !- extrapolations to them without worrying too much about 
        !- whether we've visited a specific point or not
        do k = 1 , nzv
           do j = 1 , nyv
              do i = 1, nxv, (nxv-1)
                 f_v(i,j,k) = 0.0d0 
              end do
           end do
        end do

        do k = 1 , nzv
           do j = 1 , nyv, (nyv-1)
              do i = 1, nxv
                 f_v(i,j,k) = 0.0d0 
              end do
           end do
        end do

        do k = 1 , nzv, (nzv-1)
           do j = 1 , nyv
              do i = 1, nxv
                 f_v(i,j,k) = 0.0d0 
              end do
           end do
        end do

        !-------------------------------------------------------------------
        !- 8 kinds of diagonals to use for extrapolation:
        !-   +/-x,+/-y,+/-z 
        !-   each diagonal has it's set of points for which
        !-   it can be used to extrapolate to.  We'll 
        !-   do this diagonal-specifically, so we don't have to 
        !-   keep track (too closely) where in the cube we are
        !-------------------------------------------------------------------

        ibeg(1) = ord+1
        jbeg(1) = ord+1
        kbeg(1) = ord+1
        ibeg(2) = 1
        jbeg(2) = 1
        kbeg(2) = 1

        iend(1) = nxv
        jend(1) = nyv
        kend(1) = nzv
        iend(2) = nxv - ord
        jend(2) = nyv - ord
        kend(2) = nzv - ord

        !-------------------------------------------------------------------
        ! Looping over the vectors:  
        !         [i,j,k]_dir = 1  means " + " [x,y,z]
        !         [i,j,k]_dir = 2  means " - " [x,y,z]
        !-------------------------------------------------------------------
        
        do k_dir = 1, 2
           do j_dir = 1, 2
              do i_dir = 1, 2 

                 is = ibeg(i_dir)
                 ie = iend(i_dir)
                 js = jbeg(j_dir)
                 je = jend(j_dir)
                 ks = kbeg(k_dir)
                 ke = kend(k_dir)

                 if (i_dir .eq. 1) then
                    i_offset1 = -1
                    i_offset2 = -2
                 else
                    i_offset1 = 0
                    i_offset2 = 1
                 end if

                 if (j_dir .eq. 1) then
                    j_offset1 = -1
                    j_offset2 = -2
                 else
                    j_offset1 = 0
                    j_offset2 = 1
                 end if

                 if (k_dir .eq. 1) then
                    k_offset1 = -1
                    k_offset2 = -2
                 else
                    k_offset1 = 0
                    k_offset2 = 1
                 end if

                 ! For each diagonal, you will need to integrate over
                 ! three faces, a constant i face (east or west), a
                 ! constant j face (north or south), and a constant k
                 ! face (up or down)  

                 ! do the k-constant face (top or bottom of box)

                 if (k_dir .eq. 1) then
                    k = nzv
                 else
                    k = 1
                 end if

                 do j = js, je
                   i = is
                   i_c1 = i + i_offset1
                   j_c1 = j + j_offset1
                   k_c1 = k + k_offset1
                   i_c2 = i + i_offset2
                   j_c2 = j + j_offset2
                   k_c2 = k + k_offset2
                   startv = i + nxv * (j-1) + nxv*nyv * (k-1)
                   startc = i_c1 + nxc * (j_c1-1) + nxc*nyc * (k_c1-1)
                   ind_diff = i_c2 + nxc * (j_c2-1) + nxc*nyc * (k_c2-1) 
     &                  - startc
                   
                   call dvextrap2q1d(f_c,f_v,half,startc,startv,
     &                  ntotc,ntotv,
     &                  (nxv-ord),1,1,0,ind_diff,1)
                end do

                ! do the j-constant face (top or bottom of box)

                if (j_dir .eq. 1) then
                   j = nyv
                else
                   j = 1
                end if

                 do k = ks, ke
                   i = is
                   i_c1 = i + i_offset1
                   j_c1 = j + j_offset1
                   k_c1 = k + k_offset1
                   i_c2 = i + i_offset2
                   j_c2 = j + j_offset2
                   k_c2 = k + k_offset2
                   startv = i + nxv * (j-1) + nxv*nyv * (k-1)
                   startc = i_c1 + nxc * (j_c1-1) + nxc*nyc * (k_c1-1)
                   ind_diff = i_c2 + nxc * (j_c2-1) + nxc*nyc * (k_c2-1) 
     &                  - startc

                   call dvextrap2q1d(f_c,f_v,half,startc,startv,
     &                  ntotc,ntotv,
     &                  (nxv-ord),1,1,0,ind_diff,1)
                end do
                
                ! do the i-constant face (top or bottom of box)
                ! Since the loop interior to dvextrap2q1d is in 
                ! j, we need to set the strides to the appropriate
                ! values

                if (i_dir .eq. 1) then
                   i = nxv
                else
                   i = 1
                end if

                do k = ks, ke
                   j = js
                   i_c1 = i + i_offset1
                   j_c1 = j + j_offset1
                   k_c1 = k + k_offset1
                   i_c2 = i + i_offset2
                   j_c2 = j + j_offset2
                   k_c2 = k + k_offset2
                   startv = i + nxv * (j-1) + nxv*nyv * (k-1)
                   startc = i_c1 + nxc * (j_c1-1) + nxc*nyc * (k_c1-1)
                   ind_diff = i_c2 + nxc * (j_c2-1) + nxc*nyc * (k_c2-1) 
     &                  - startc
                   
                   call dvextrap2q1d(f_c,f_v,half,startc,startv,
     &                  ntotc,ntotv,
     &                  (nyv-ord),nxc,nxv,0,ind_diff,1)
                end do
       
              end do
           end do
        end do

        !- Now normalize the appropriate points by the number 
        !- of extrapolations performed for that point. 
        !- On any given face, the points in the middle are 
        !- extrapolated by four diagonals.  There are four 
        !- regions where 2 diagonals contribute, and four 
        !- where only one contributes.  The pattern looks like
        !- this (roughly):  

        !-           --------------------------
        !-           |  1  |     2      |  1  |
        !-           --------------------------
        !-           |     |            |     |
        !-           |     |            |     |
        !-           |  2  |     4      |  2  |
        !-           |     |            |     |
        !-           |     |            |     |
        !-           --------------------------
        !-           |  1  |     2      |  1  |
        !-           --------------------------

        !- In handling the faces below, we will avoid all edges
        !- and corners.

        ! x-faces:
        do k = 1+ord, nzv-ord
           do j = 1+ord, nyv-ord
              do i = 1, nxv, (nxv-1)
                 f_v(i,j,k) = fourth * f_v(i,j,k)
              end do
           end do
        end do

        do k = 2, ord
           do j = 1+ord, nyv-ord
              do i = 1, nxv, (nxv-1)
                 f_v(i,j,k) = half * f_v(i,j,k)
              end do
           end do
        end do

        do k = nzv-ord+1, nzv-1
           do j = 1+ord, nyv-ord
              do i = 1, nxv, (nxv-1)
                 f_v(i,j,k) = half * f_v(i,j,k)
              end do
           end do
        end do

        do k = 1+ord, nzv-ord
           do j = 2, ord
              do i = 1, nxv, (nxv-1)
                 f_v(i,j,k) = half * f_v(i,j,k)
              end do
           end do

           do j = nyv-ord+1, nyv-1
              do i = 1, nxv, (nxv-1)
                 f_v(i,j,k) = half * f_v(i,j,k)
              end do
           end do
        end do

        ! y-faces:
        do k = 1+ord, nzv-ord
           do j = 1, nyv, (nyv-1)
              do i = 1+ord, nxv-ord
                 f_v(i,j,k) = fourth * f_v(i,j,k)
              end do
           end do
        end do

        do k = 2, ord
           do j = 1, nyv, (nyv-1)
              do i = 1+ord, nxv-ord
                 f_v(i,j,k) = half * f_v(i,j,k)
              end do
           end do
        end do

        do k = nzv-ord+1, nzv-1
           do j = 1, nyv, (nyv-1)
              do i = 1+ord, nxv-ord
                 f_v(i,j,k) = half * f_v(i,j,k)
              end do
           end do
        end do

        do k = 1+ord, nzv-ord
           do j = 1, nyv, (nyv-1)
              do i = 2, ord
                 f_v(i,j,k) = half * f_v(i,j,k)
              end do

              do i = nxv-ord+1, nxv-1
                 f_v(i,j,k) = half * f_v(i,j,k)
              end do
           end do
        end do

        ! z-faces:

        do k = 1, nzv, (nzv-1)
           do j = 1+ord, nyv-ord
              do i = 1+ord, nxv-ord
                 f_v(i,j,k) = fourth * f_v(i,j,k)
              end do
           end do

           do j = 2, ord
              do i = 1+ord, nxv-ord
                 f_v(i,j,k) = half * f_v(i,j,k)
              end do
           end do

           do j = nyv-ord+1, nyv-1
              do i = 1+ord, nxv-ord
                 f_v(i,j,k) = half * f_v(i,j,k)
              end do
           end do

           do j = 1+ord, nyv-ord
              do i = 2, ord
                 f_v(i,j,k) = half * f_v(i,j,k)
              end do
           end do

           do j = 1+ord, nyv-ord
              do i = nxv-ord+1, nxv-1
                 f_v(i,j,k) = half * f_v(i,j,k)
              end do
           end do

        end do

        !- In addition, the edges of the box have been 
        !- double counted, and the vertices of the box
        !- have been triple counted.

        do k = 2, ord
           do j = 1, nyv, (nyv-1)
              do i = 1, nxv, (nxv-1)
                 f_v(i,j,k) = f_v(i,j,k)/2.d0
              end do
           end do
        end do

        do k = ord+1, nzv-ord
           do j = 1, nyv, (nyv-1)
              do i = 1, nxv, (nxv-1)
                 ! here one factor comes from 
                 ! the double counting, and the other
                 ! comes from the fact that two
                 ! diagonals contributed.
                 f_v(i,j,k) = f_v(i,j,k)/2.d0
                 f_v(i,j,k) = f_v(i,j,k)/2.d0
              end do
           end do
        end do

        do k = nzv-ord+1, nzv-1
           do j = 1, nyv, (nyv-1)
              do i = 1, nxv, (nxv-1)
                 f_v(i,j,k) = f_v(i,j,k)/2.d0
              end do
           end do
        end do

        do k = 1, nzv, (nzv-1)
           do j = 2, ord
              do i = 1, nxv, (nxv-1)
                 f_v(i,j,k) = f_v(i,j,k)/2.d0
              end do
           end do

           do j = ord+1, nyv-ord
              do i = 1, nxv, (nxv-1)
                 f_v(i,j,k) = f_v(i,j,k)/2.d0
                 f_v(i,j,k) = f_v(i,j,k)/2.d0
              end do
           end do
           
           do j = nyv-ord+1, nyv-1
              do i = 1, nxv, (nxv-1)
                 f_v(i,j,k) = f_v(i,j,k)/2.d0
              end do
           end do
        end do

        do k = 1, nzv, (nzv-1)
           do j = 1, nyv, (nyv-1)
              do i = 2, ord
                 f_v(i,j,k) = f_v(i,j,k)/2.d0
              end do

              do i = ord+1,nxv-ord
                 f_v(i,j,k) = f_v(i,j,k)/2.d0
                 f_v(i,j,k) = f_v(i,j,k)/2.d0
              end do

              do i = nxv-ord+1,nxv-1
                 f_v(i,j,k) = f_v(i,j,k)/2.d0
              end do
           end do
        end do

        do k = 1, nzv, (nzv-1)
           do j = 1, nyv, (nyv-1)
              do i = 1, nxv, (nxv-1)
                 f_v(i,j,k) = f_v(i,j,k)/3.d0
              end do
           end do
        end do
        
        return
        end

c---------------------------------------------------------------------------
c For 2D arrays, performs a 2nd-order interpolation (averaging/extrap.) of 
c cell-centered function (f_c) into a vertex-centered function (f_v). 
c 
c This routine assumes that the excision zone will not come within two cell
c widths of an AMR boundary.  Otherwise, the technique of using diagonals
c to extrapolate to the surface will fail.  
c
c  -- does not check to see if  ny* > 1
c---------------------------------------------------------------------------
        subroutine dm_c_to_v_2d(f_c,f_v,chr_v,nxc,nyc,nxv,nyv,ex,do_ex)
        implicit none
        integer nxc,nyc,nxv,nyv,ord,do_ex
        real*8 f_c(nxc,nyc),f_v(nxv,nyv)
        real*8 chr_v(nxv,nyv)
        real*8 ex
        integer i,j, ntotc, ntotv, startc, startv
        integer i_dir, j_dir, i_c1, j_c1, i_c2, j_c2
        integer ibeg(2), jbeg(2), iend(2),jend(2)
        integer i_offset1, j_offset1, i_offset2, j_offset2
        integer ind_diff, is, ie, js, je

        real*8 half,fourth
        parameter  ( fourth       =  0.2500 0000 0000 0000 d0,
     &               half         =  0.5000 0000 0000 0000 d0 )
        
        !-------------------------------------------------------------------
        !  Order of interpolation and extrapolation
        !-------------------------------------------------------------------
        ord = 2

        !-------------------------------------------------------------------
        !  Interior Points: (averaging nearest-neighbors)
        !-------------------------------------------------------------------

        do j=2,nyv-1
           do i=2,nxv-1
              if (do_ex.eq.1 .and. chr_v(i,j).eq.ex) then
                 f_v(i,j) = 0.d0
              else
                 f_v(i,j) = fourth * ( f_c(i-1, j-1) +
     &                                 f_c(i  , j-1) +
     &                                 f_c(i-1, j  ) +
     &                                 f_c(i  , j  )  )
              end if
           end do
        end do

        !-------------------------------------------------------------------
        !  Edge points: (average the extrapolations of 2 diagonals)
        !-------------------------------------------------------------------

        ntotv = nxv*nyv
        ntotc = nxc*nyc

        !- zero-out edges and corners so that we can just add 
        !- extrapolations to them without worrying too much about 
        !- whether we've visited a specific point or not
        do j = 1 , nyv
           do i = 1, nxv, (nxv-1)
              f_v(i,j) = 0.0d0 
           end do
        end do

        do j = 1 , nyv, (nyv-1)
           do i = 1, nxv
              f_v(i,j) = 0.0d0 
           end do
        end do

        !-------------------------------------------------------------------
        !- 4 kinds of diagonals to use for extrapolation:
        !-   NE,SE,NW, and SW
        !- (This assumes j increases to the north and i 
        !- increases to the east.)
        !-------------------------------------------------------------------
        
        ibeg(1) = ord+1
        jbeg(1) = ord+1
        ibeg(2) = 1
        jbeg(2) = 1

        iend(1) = nxv
        jend(1) = nyv
        iend(2) = nxv - ord
        jend(2) = nyv - ord

        !-------------------------------------------------------------------
        ! Looping over the vectors:  
        !         [i,j]_dir = 1  means " + " [x,y]
        !         [i,j]_dir = 2  means " - " [x,y]
        !-------------------------------------------------------------------
        do j_dir = 1, 2
           do i_dir = 1, 2 

              is = ibeg(i_dir)
              ie = iend(i_dir)
              js = jbeg(j_dir)
              je = jend(j_dir)

              if (i_dir .eq. 1) then
                 i_offset1 = -1
                 i_offset2 = -2
              else
                 i_offset1 = 0
                 i_offset2 = 1
              end if
              
              if (j_dir .eq. 1) then
                 j_offset1 = -1
                 j_offset2 = -2
              else
                 j_offset1 = 0
                 j_offset2 = 1
              end if

              ! For each diagonal, you will need to integrate over
              ! two faces, a constant i face (east or west), a
              ! constant j face (north or south)
              ! do the j-constant face (top or bottom of box)

              if (j_dir .eq. 1) then
                 j = nyv
              else
                 j = 1
              end if

              i = is
              i_c1 = i + i_offset1
              j_c1 = j + j_offset1
              i_c2 = i + i_offset2
              j_c2 = j + j_offset2
              startv = i + nxv * (j-1) 
              startc = i_c1 + nxc * (j_c1-1) 
              ind_diff = i_c2 + nxc * (j_c2-1) - startc
              call dvextrap2q1d(f_c,f_v,half,startc,startv,
     &             ntotc,ntotv,
     &             (nxv-ord),1,1,0,ind_diff,1)
                
              ! do the i-constant face (east or west sides of box)
              ! Since the loop interior to dvextrap2q1d is in 
              ! j, we need to set the strides to the appropriate
              ! values

              if (i_dir .eq. 1) then
                 i = nxv
              else
                 i = 1
              end if
              
              j = js
              i_c1 = i + i_offset1
              j_c1 = j + j_offset1
              i_c2 = i + i_offset2
              j_c2 = j + j_offset2
              startv = i + nxv * (j-1) 
              startc = i_c1 + nxc * (j_c1-1) 
              ind_diff = i_c2 + nxc * (j_c2-1) - startc
              call dvextrap2q1d(f_c,f_v,half,startc,startv,
     &             ntotc,ntotv,
     &             (nyv-ord),nxc,nxv,0,ind_diff,1)
       
           end do
        end do

        !- Now normalize the appropriate points by the number 
        !- of extrapolations performed for that point. 

        do j = 1 , nyv, (nyv-1)
           do i = 1+ord, nxv-ord
              f_v(i,j) = half * f_v(i,j) 
           end do
        end do
        do j = 1+ord, nyv-ord
           do i = 1, nxv, (nxv-1)
              f_v(i,j) = half * f_v(i,j) 
           end do
        end do

        !- We ended up visiting each of the corners twice, so 
        !- correct for this as well.
        do j = 1, nyv, (nyv-1)
           do i = 1, nxv, (nxv-1)
              f_v(i,j) = half * f_v(i,j)
           end do
        end do

        return
        end

c---------------------------------------------------------------------------
c For 1D arrays, performs a 2nd-order interpolation (averaging/extrap.) of 
c cell-centered function (f_c) into a vertex-centered function (f_v). 
c 
c This routine assumes that the excision zone will not come within two cell
c widths of an AMR boundary.  Otherwise, the technique of using diagonals
c to extrapolate to the surface will fail.  
c
c  -- does not check to see if  ny* > 1
c---------------------------------------------------------------------------
        subroutine dm_c_to_v_1d(f_c,f_v,chr_v,nxc,nxv,ex,do_ex)
        implicit none
        integer nxc,nxv,ord,do_ex
        real*8 f_c(nxc),f_v(nxv),chr_v(nxv)
        real*8 ex
        integer i, ntotc, ntotv, startc, startv
        integer i_dir
        integer ibeg(2), iend(2)
        integer ind_diff1, ind_diff2

        real*8 half
        parameter  ( half         =  0.5000 0000 0000 0000 d0 )

        
        !-------------------------------------------------------------------
        !  Order of interpolation and extrapolation
        !-------------------------------------------------------------------
        ord = 2

        !-------------------------------------------------------------------
        !  Interior Points: (averaging nearest-neighbors)
        !-------------------------------------------------------------------

        do i=2,nxv-1
           if (do_ex.eq.1 .and. chr_v(i).eq.ex) then
              f_v(i) = 0.d0
           else
              f_v(i) = half * ( f_c(i-1) + f_c(i) )
           end if
        end do

        !-------------------------------------------------------------------
        !  Edge points: (extrapolations)
        !-------------------------------------------------------------------
        
        !-- left boundary:
        call dvextrap2q1d(f_c,f_v,half,1,1,nxc,nxv,1,1,1,0,1,0)
        !-- right boundary:
        call dvextrap2q1d(f_c,f_v,half,nxc,nxv,nxc,nxv,1,1,1,0,-1,0)

        return
        end

c-----------------------------------------------------------------------
c Performs 2nd-order extrapolation to point v2(i) from v1(i+s11) 
c  and v1(i+s12).  It is assumed that |s11|<|s12|, so that v1(i+s12) is 
c  furthest away from the point at which we are extrapolating. 
c  dr = (x - x(i+s11)/dx) where  x() some monotonically-increasing 
c  uniformly-discretized coordinate w/ index i  and with spacing dx, 
c  i.e. dr is the dimensionless spacing between the extrapolatee
c  point and the nearest extrapolater point.  
c
c  saveold = 0  (overwrite v2)
c          = otherwise (add on to v2)
c-----------------------------------------------------------------------
        subroutine dvextrap2q1d(v1,v2,dr,i1,i2,n1,n2,n,s1,s2,di1,di2,
     &                          saveold)
        implicit none
        integer i1,i2,n1,n2,n,s1,s2,di1,di2,saveold
        real*8  v1(n1),v2(n2),dr

        integer i,j1,j2
        real*8  dr2

        
        if( saveold .eq. 0 ) then
           do i = 1, n
              j1 = (i-1)*s1 + i1
              j2 = (i-1)*s2 + i2
              v2(j2) =  v1(j1+di1)
     &                       + dr * ( v1(j1+di1) - v1(j1+di2) )
           end do
        else
           do i = 1, n 
              j1 = (i-1)*s1 + i1
              j2 = (i-1)*s2 + i2
              v2(j2) =  v2(j2)  + v1(j1+di1)
     &                       + dr * ( v1(j1+di1) - v1(j1+di2) )
           end do
        end if

        return 
        end 

c---------------------------------------------------------------------------
c the `inverse' of the above, where the copy is from vertex to cell 
c centered
c 
c  -- drops to first order near excision zone.
c
c order is as enumerated in 'ops.inc'
c  -- only ord=2 (averaging) implemented
c---------------------------------------------------------------------------
        subroutine dm_v_to_c(f_v,f_c,chr_v,chr_c,
     &                       nxv,nyv,nzv,nxc,nyc,nzc,
     &                       ord,ex,do_ex,dim)
        implicit none
        integer nxc,nyc,nzc,nxv,nyv,nzv,do_ex,ord,dim
        real*8 f_c(nxc,nyc,nzc),f_v(nxv,nyv,nzv)
        real*8 chr_c(nxc,nyc,nzc),chr_v(nxv,nyv,nzv)
        real*8 ex

        include 'ops.inc'

        if (ord .eq. PAMR_V_TO_C_NO_TRANSFER) then
           return
        else if (.not.(ord.eq.PAMR_V_TO_C_SECOND_ORDER)) then
           write(*,*) 'dm_v_to_c: ord = ', ord, ' not yet implemented'
           stop
        end if
           
        if (dim.eq.1) then
           call dm_v_to_c_1d(f_v,f_c,chr_v,chr_c,nxv,nxc,ex,do_ex)
        else if (dim .eq. 2) then
           call dm_v_to_c_2d(f_v,f_c,chr_v,chr_c,nxv,nyv,nxc,nyc,
     &          ex,do_ex)
        else if (dim .eq. 3) then
           call dm_v_to_c_3d(f_v,f_c,chr_v,chr_c,nxv,nyv,nzv,
     &          nxc,nyc,nzc,ex,do_ex)
        else 
           write(*,*) 'dm_v_to_c: dim = ', dim, ' not supported'
           stop
        end if
        
        return
        end
c---------------------------------------------------------------------------
c For 3D arrays, performs a 2nd-order interpolation (averaging) of 
c vertex-centered function (f_v) into a cell-centered function (f_c). 
c 
c  -- does not check to see if  ny*,nz* > 1
c---------------------------------------------------------------------------
        subroutine dm_v_to_c_3d(f_v,f_c,chr_v,chr_c,nxv,nyv,nzv,
     &     nxc,nyc,nzc,ex,do_ex)
        implicit none
        integer nxc,nyc,nzc,nxv,nyv,nzv,do_ex
        real*8 f_c(nxc,nyc,nzc),f_v(nxv,nyv,nzv)
        real*8 chr_c(nxc,nyc,nzc),chr_v(nxv,nyv,nzv)
        real*8 ex,n_ave
        integer i,j,k,li,lj,lk

        do k=1,nzc
           do j=1,nyc
              do i=1,nxc
                 if (do_ex .eq. 0) then
                    f_c(i,j,k) = (f_v(i,j,k) + f_v(i+1,j,k) + 
     &                   f_v(i,j+1,k) + f_v(i+1,j+1,k) + 
     &                   f_v(i,j,k+1) + f_v(i+1,j,k+1) + 
     &                   f_v(i,j+1,k+1) + f_v(i+1,j+1,k+1))
                    f_c(i,j,k) = f_c(i,j,k)/8.d0
                 else
                    if (chr_c(i,j,k).eq.ex) then
                       f_c(i,j,k) = 0.d0
                    else
                       f_c(i,j,k) = 0.d0
                       n_ave      = 0.d0
                       
                       do lk=k,k+1
                          do lj=j,j+1
                             do li=i,i+1
                                if (chr_v(li,lj,lk).ne.ex) then
                                   f_c(i,j,k) = f_c(i,j,k)+f_v(li,lj,lk)
                                   n_ave      = n_ave     +1.d0
                                end if
                             end do
                          end do
                       end do
                       
                       if (n_ave.eq.0.d0) then
                          f_c(i,j,k) = 0.d0
                       else
                          f_c(i,j,k) = f_c(i,j,k)/n_ave
                       end if
                    end if
                 end if
              end do
           end do
        end do

        return
        end

c---------------------------------------------------------------------------
c Like dm_v_to_c_3d() but for rank-2 arrays (ignoring dummy ranks), 
c interpolation (averaging) of vertex-centered function (f_v) into a 
c cell-centered function (f_c). 
c
c  -- does not check to see if  ny* > 1
c---------------------------------------------------------------------------
        subroutine dm_v_to_c_2d(f_v,f_c,chr_v,chr_c,nxv,nyv,nxc,nyc,
     &     ex,do_ex)
        implicit none
        integer nxc,nyc,nxv,nyv,do_ex
        real*8 f_c(nxc,nyc),f_v(nxv,nyv),chr_c(nxc,nyc),chr_v(nxv,nyv)
        real*8 ex,n_ave
        integer i,j,li,lj

        do j=1,nyc
           do i=1,nxc
              if (do_ex.ne.0 .and. chr_c(i,j).eq.ex) then
                 f_c(i,j) = 0.d0
              else
                 f_c(i,j) = 0.d0
                 n_ave    = 0.d0

                 do li=i,i+1
                    do lj=j,j+1
                       if (do_ex.eq.0 .or. chr_v(li,lj).ne.ex) then
                          f_c(i,j) = f_c(i,j) + f_v(li,lj)
                          n_ave    = n_ave    + 1.d0
                       end if
                    end do
                 end do

                 if (n_ave.eq.0.d0) then
                    f_c(i,j) = 0.d0
                 else
                    f_c(i,j) = f_c(i,j)/n_ave
                 end if
              end if
           end do
        end do

        return
        end
c---------------------------------------------------------------------------
c Like dm_v_to_c_3d() but for rank-1 arrays (ignoring dummy ranks), 
c performs a 2nd-order interpolation (averaging) of vertex-centered 
c function (f_v) into a cell-centered function (f_c). 
c---------------------------------------------------------------------------
        subroutine dm_v_to_c_1d(f_v,f_c,chr_v,chr_c,nxv,nxc,ex,do_ex)
        implicit none
        integer nxc,nxv,do_ex
        real*8 f_c(nxc),f_v(nxv),chr_c(nxc),chr_v(nxv)
        real*8 ex, n_ave
        integer i,li

        do i=1,nxc
           if (do_ex.ne.0 .and. chr_c(i).eq.ex) then
              f_c(i) = 0.d0
           else
              f_c(i) = 0.d0
              n_ave  = 0.d0

              do li=i,i+1
                 if (do_ex.eq.0 .or. chr_v(li).ne.ex) then
                    f_c(i) = f_c(i) + f_v(li)
                    n_ave  = n_ave  + 1.d0
                 end if
              end do

              if (n_ave.eq.0.d0) then
                 f_c(i) = 0.d0
              else
                 f_c(i) = f_c(i)/n_ave
              end if
           end if
        end do

        return
        end

c-----------------------------------------------------------------------
c coarsens via simple averaging from fine grid uf of size [nx_f,ny_f,nz_f]
c to coarse grid uc of size [nx_c,ny_c,nz_c], with spatial refinement 
c ratio rho. 
c
c We average over unexcised fine cells only.  If all fine cells are excised, 
c then we insert a zero.  But this coarse cell will, of course, also be
c excised.  
c
c offset into fine (coarse) grid is [i_f,j_f,k_f] ([i_c,j_c,k_c]), and size
c of region on coarse grid is [ni_c,nj_c,nk_c]
c dim is the dimensionality of the system        
c-----------------------------------------------------------------------
      subroutine dmavg3d_c(u_f,u_c,chr,nx_f,ny_f,nz_f,
     &                     nx_c,ny_c,nz_c,
     &                     i_f,j_f,k_f,i_c,j_c,k_c,ni_c,nj_c,nk_c,
     &                     rho,dim,ex,do_ex)
      implicit none
      integer nx_f,ny_f,nz_f,nx_c,ny_c,nz_c,do_ex
      integer i_f,j_f,k_f,i_c,j_c,k_c,ni_c,nj_c,nk_c,rho
      real*8 u_c(nx_c,ny_c,nz_c),u_f(nx_f,ny_f,nz_f)
      real*8 chr(nx_f,ny_f,nz_f)
      real*8 ex
      logical       ltrace
      parameter   ( ltrace = .false. )
      integer dim
      integer i,j,k,ii,jj,kk,li,lj,lk
      integer maxi,maxj,maxk
      real*8 n_ave
      include 'ops.inc'
      
      if (rho .ne. 2) then
         write(*,*) 'dmavg3d_c: only a 2:1 ratio supported'
         stop
      end if
      
      if( ni_c .gt. nx_c .or. nj_c .gt. ny_c .or. nk_c .gt. nz_c) then
         write(*,*) 'Offset grid larger than maximum grid size'
         write(*,*) 'ni_c = ', ni_c, 'nx_c = ',nx_c
         write(*,*) 'nj_c = ', nj_c, 'ny_c = ',ny_c
         write(*,*) 'nk_c = ', nk_c, 'nz_c = ',nz_c
         stop
      end if
      
      if (ltrace) then
         write(*,*) "shift c = ",i_c, j_c, k_c
         write(*,*) "shift f = ",i_f, j_f, k_f
         write(*,*) "cw = ",ni_c, nj_c, nk_c
         write(*,*) "cwmax= ",nx_c, ny_c, nz_c
      end if

      if(dim .eq. 1) then
         do i = 1,ni_c
            
            ii = i_f + 2*(i) - 1
            
            n_ave=0.d0
            u_c(i-1+i_c,nj_c,nk_c) = 0.d0
            do li=ii-1,ii
               if (do_ex.eq.0 .or. chr(li,1,1).ne.ex) then
                  n_ave = n_ave + 1.d0
                  u_c(i-1+i_c,nj_c,nk_c) = u_c(i-1+i_c,nj_c,nk_c) +
     &                                     u_f(li,1,1)
               end if
            end do

            if (n_ave.eq.0.d0) then
               u_c(i-1+i_c,nj_c,nk_c) = 0.d0
            else
               u_c(i-1+i_c,nj_c,nk_c) = u_c(i-1+i_c,nj_c,nk_c)/n_ave
            end if
         enddo
         
      else if(dim .eq. 2) then
         
         do i = 1,ni_c
            do j = 1, nj_c
               
               ii = i_f + 2*(i) - 1 
               jj = j_f + 2*(j) - 1 

               n_ave=0.d0
               u_c(i-1+i_c,j-1+j_c,nk_c) = 0.d0
               do li=ii-1,ii
                  do lj=jj-1,jj
                     if (do_ex.eq.0 .or. chr(li,lj,1).ne.ex) then
                        n_ave = n_ave + 1.d0
                        u_c(i-1+i_c,j-1+j_c,nk_c) = 
     &                       u_c(i-1+i_c,j-1+j_c,nk_c) +
     &                       u_f(li,lj,1)
                     end if
                  end do
               end do

               if (n_ave.eq.0.d0) then
                  u_c(i-1+i_c,j-1+j_c,nk_c) = 0.d0
               else
                  u_c(i-1+i_c,j-1+j_c,nk_c) = 
     &                 u_c(i-1+i_c,j-1+j_c,nk_c)/n_ave
               end if
            enddo
         enddo
         
      else if(dim .eq. 3) then
         
         do i = 1,ni_c
            do j = 1, nj_c
               do k = 1, nk_c
                  
                  ii = i_f + 2*(i) - 1
                  jj = j_f + 2*(j) - 1
                  kk = k_f + 2*(k) - 1

                  n_ave=0.d0
                  u_c(i-1+i_c,j-1+j_c,k-1+k_c) = 0.d0

                  do li=ii-1,ii
                     do lj=jj-1,jj
                        do lk=kk-1,kk 
                           if (do_ex.eq.0 .or. chr(li,lj,lk).ne.ex) then
                              n_ave=n_ave+1.d0
                              u_c(i-1+i_c,j-1+j_c,k-1+k_c) = 
     &                             u_c(i-1+i_c,j-1+j_c,k-1+k_c) +
     &                             u_f(li,lj,lk)
                           end if
                        end do
                     end do
                  end do

                  if (n_ave.eq.0.d0) then
                     u_c(i-1+i_c,j-1+j_c,k-1+k_c) = 0.d0
                  else
                     u_c(i-1+i_c,j-1+j_c,k-1+k_c) = 
     &                    u_c(i-1+i_c,j-1+j_c,k-1+k_c)/n_ave
                  end if                  
               enddo
            enddo
         enddo
         
      else
         write(*,*) 'The dimension of the system must be 3 or less'
         stop
      end if
      
      return
      end

c-----------------------------------------------------------------------
c coarsens via simple addition from fine grid uf of size [nx_f,ny_f,nz_f]
c to coarse grid uc of size [nx_c,ny_c,nz_c], with spatial refinement 
c ratio rho. 
c
c This is to be used for extensive quantities, such as the restmass 
c inside a given cell.  Note that this subroutine ignores excision.
c This is because excised cells will have been zeroed out and will make
c no contribution.
c-----------------------------------------------------------------------
      subroutine dmadd3d_c(u_f,u_c,nx_f,ny_f,nz_f,
     &                     nx_c,ny_c,nz_c,
     &                     i_f,j_f,k_f,i_c,j_c,k_c,ni_c,nj_c,nk_c,
     &                     rho,dim)
      implicit none
      integer nx_f,ny_f,nz_f,nx_c,ny_c,nz_c
      integer i_f,j_f,k_f,i_c,j_c,k_c,ni_c,nj_c,nk_c,rho
      real*8 u_c(nx_c,ny_c,nz_c),u_f(nx_f,ny_f,nz_f)
      
      logical       ltrace
      parameter   ( ltrace = .false. )

      integer dim
      integer i,j,k,ii,jj,kk
      integer maxi,maxj,maxk
      include 'ops.inc'

      if (rho .ne. 2) then
         write(*,*) 'dmadd3d_c: only a 2:1 ratio supported'
         stop
      end if
      
      if( ni_c .gt. nx_c .or. nj_c .gt. ny_c .or. nk_c .gt. nz_c) then
         write(*,*) 'Offset grid larger than maximum grid size'
         write(*,*) 'ni_c = ', ni_c, 'nx_c = ',nx_c
         write(*,*) 'nj_c = ', nj_c, 'ny_c = ',ny_c
         write(*,*) 'nk_c = ', nk_c, 'nz_c = ',nz_c
         stop
      end if
      
      if (ltrace) then
         write(*,*) "shift c = ",i_c, j_c, k_c
         write(*,*) "shift f = ",i_f, j_f, k_f
         write(*,*) "cw = ",ni_c, nj_c, nk_c
         write(*,*) "cwmax= ",nx_c, ny_c, nz_c
      end if
      
      if(dim .eq. 1) then
         do i = 1,ni_c
            
            ii = i_f + 2*(i) - 1
            
            u_c(i-1+i_c,nj_c,nk_c) = ( u_f(ii,1,1)
     &           + u_f(ii-1,1,1))
            
         enddo
         
      else if(dim .eq. 2) then
         
         do i = 1,ni_c
            do j = 1, nj_c
               
               ii = i_f + 2*(i) - 1 
               jj = j_f + 2*(j) - 1 
               
               u_c(i-1+i_c,j-1+j_c,nk_c) = ( u_f(ii,jj,1) 
     &              + u_f(ii-1,jj,1)
     &              + u_f(ii,jj-1,1) + u_f(ii-1,jj-1,1) )
               
            enddo
         enddo
         
      else if(dim .eq. 3) then
         
         do i = 1,ni_c
            do j = 1, nj_c
               do k = 1, nk_c
                  
                  ii = i_f + 2*(i) - 1
                  jj = j_f + 2*(j) - 1
                  kk = k_f + 2*(k) - 1
c                  write(*,*) i,j,k, "doubles ", ii,jj,kk
                  
                  u_c(i-1+i_c,j-1+j_c,k-1+k_c) = ( u_f(ii,jj,kk)
     &                 + u_f(ii-1,jj,kk)
     &                 + u_f(ii,jj-1,kk) + u_f(ii,jj,kk-1)
     &                 + u_f(ii-1,jj-1,kk) + u_f(ii-1,jj,kk-1)
     &                 + u_f(ii,jj-1,kk-1) + u_f(ii-1,jj-1,kk-1))
                  
               enddo
            enddo
         enddo
      else
         write(*,*) 'The dimension of the system must be 3 or less'
         stop
      end if
      
      return
      end

      ! A function for computing the mc slope given three
      ! adjacent points.

      real*8 function mc_slope(fm1,f,fp1)
      real*8 fm1,f,fp1
      real*8 dfp, dfm, dfc, temp
      dfp = fp1 - f
      dfm = f - fm1
      dfc = 0.5d0 * (dfp + dfm)
      
      if (dfp*dfm .lt. 0.d0) then 
         mc_slope = 0.d0
      else
         if (abs(2.d0*dfp) .lt. abs(2.d0*dfm)) then
            temp = 2.d0*dfp
         else
            temp = 2.d0*dfm
         end if
         
         if (abs(temp) .lt. abs(dfc)) then
            mc_slope = temp
         else
            mc_slope = dfc
         end if
      end if

      return

      end


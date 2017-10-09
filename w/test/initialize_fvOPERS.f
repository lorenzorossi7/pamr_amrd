        subroutine initialize_fv(fv, amp, 
     &                               x0,   y0,   z0,
     &                               sigx, sigy, sigz,  
     &                               x,    y,    z,
     &                               nx,   ny,   nz,
     &                               type)

        implicit      none

        logical       ltrace
        parameter   ( ltrace = .true. )
 
        character*13  cdnm
        parameter   ( cdnm = 'initialize_fv' )

        integer       Nx,   Ny,  Nz
        real*8        fv  (Nx,  Ny,  Nz) 
        real*8        x( Nx )
        real*8        y( Ny )
        real*8        z( Nz )
        real*8        x0,   y0,   z0 
        real*8        sigx, sigy, sigz
        integer       i,    j,    k

        real*8        amp
        integer       type

c=======================================================================
c=======================================================================
c=======================================================================
c=======================================================================
c=======================================================================

        if ( ltrace ) then
            write(0,*) cdnm,':   Nx = ', Nx 
            write(0,*) cdnm,':   Ny = ', Ny
            write(0,*) cdnm,':   Nz = ', Nz 
            write(0,*) cdnm,':   x0 = ', x0
            write(0,*) cdnm,':   y0 = ', y0
            write(0,*) cdnm,':   z0 = ', z0
            write(0,*) cdnm,': sigx = ', sigx
            write(0,*) cdnm,': sigy = ', sigy
            write(0,*) cdnm,': sigz = ', sigz
            write(0,*) cdnm,':  amp = ', amp
            write(0,*) cdnm,': type = ', type
        endif
        if ( type .eq. 0 ) then
          do i = 1, Nx
             do j = 1, Ny
                do k = 1, Nz
                   fv(i,j,k) = amp * exp(-(x(i)-x0)**2/sigx**2) * 
     &                               exp(-(y(j)-y0)**2/sigy**2) * 
     &                               exp(-(z(k)-z0)**2/sigz**2) 
                enddo
             enddo
          enddo
        elseif ( type .eq. 1 ) then
          do i = 1, Nx 
             do j = 1, Ny 
                do k = 1, Nz 
                   fv(i,j,k) = x(i) + 2.0d0*y(j) + 3.0d0*z(k)
                enddo
             enddo
          enddo
        else
          write(0,*) cdnm,': type = ',type, ' NOT IMPLEMENTED'
          stop
        endif
        if ( ltrace ) then
          write(0,*) cdnm,': Leaving Routine ...'
        endif
        return
        end


        subroutine output_f(f, fname,  level,
     &                      x,    y,    z,
     &                      nx,   ny,   nz,
     &                      rank)

        implicit      none

        logical       ltrace
        parameter   ( ltrace = .true. )
 
        character*13  cdnm
        parameter   ( cdnm = 'output_f' )

        integer           Nx, Ny, Nz
        real*8        f ( Nx, Ny, Nz ) 

        integer       rank

        real*8        x ( Nx )
        real*8        y ( Ny )
        real*8        z ( Nz )

        real*8        bbox(6)
        integer       shape(3)    

        integer       level
 
        integer       ret,   gft_out_brief

        character*5   fname
        character*7   filename

c=======================================================================
c=======================================================================
c=======================================================================
c=======================================================================
c=======================================================================

        if ( ltrace ) then
            write(0,*) cdnm,':   fname = ', fname
            write(0,*) cdnm,':   Nx = ', Nx 
            write(0,*) cdnm,':   Ny = ', Ny
            write(0,*) cdnm,':   Nz = ', Nz 
        endif
        shape(1) = Nx
        shape(2) = Ny
        shape(3) = Nz

        bbox(1)  = x(1)
        bbox(2)  = x(Nx)
        bbox(3)  = y(1)
        bbox(4)  = y(Ny)
        bbox(5)  = z(1)
        bbox(6)  = z(Nz)
        call gft_out_set_bbox(bbox,rank)
901     format(a,'_',i1)
        write(filename,901) fname,level
        ret = gft_out_brief(filename,0.0d0,shape,rank,f)

        return
        end


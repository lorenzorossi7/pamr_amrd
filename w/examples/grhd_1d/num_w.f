c----------------------------------------------------------------------
c numerical routines for the General Relativistic Hydro Code in 1D
c----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     Do one hydro step
c     On the metric variables being passed in
c     Subscript ft denotes "flux time"--the physical time from
c     which the fluxes are being calculated.  The unsubscripted
c     metric variables correspond to the final time for this 
c     evolution substep.
c-----------------------------------------------------------------------
      subroutine hydro_1step(d,sr,tau,d_p,sr_p,tau_p,
     &     rho,v,u,P,T,cs,jr,rhops,ham_source,x,x_c,gamma,limiter,Ye,dt,
     &     rho0,Nx,psi,psi_c,alpha,alpha_c,beta,beta_c,
     &     fc_mask,d_fcs,sr_fcs,tau_fcs,iter,eos_flag,phys_bdy,
     &     deltat1,deltat2,rho_atm,p_atm,u_atm,t_atm,lorentz_max)
      implicit none
      integer Nx,Ny,Nz
c     Conserved variables
      real*8 d(Nx),d_p(Nx)
      real*8 sr(Nx),sr_p(Nx)
      real*8 tau(Nx),tau_p(Nx)
c     Primitive variables
      real*8 v(Nx),u(Nx),rho(Nx)
      real*8 P(Nx), T(Nx), cs(Nx)
      real*8 jr(Nx), rhops(Nx), ham_source(Nx)
      real*8 psi(Nx+1),alpha(Nx+1),beta(Nx+1)
      real*8 psi_c(Nx),alpha_c(Nx),beta_c(Nx)
      
      real*8 x(Nx+1),x_c(Nx),gamma,Ye, rho0
      real*8 deltat1,deltat2, dt
      integer eos_flag, iter
      integer phys_bdy(6)
      real*8 rho_atm, p_atm, u_atm, t_atm
      
c     Flux correction stuff.
      integer fc_mask(Nx)
      integer lfc_mask, save_fcs
      real*8 d_fcs(Nx), tau_fcs(Nx), sr_fcs(Nx)
c     the local sign with which fluxes enter the flux correction.  
      real*8 lsign  
                                
      include 'fc_mask.inc'

c     Local variables.
      real*8 rhol(Nx),rhor(Nx)
      real*8 vl(Nx),vr(Nx)
      real*8 Pl(Nx),Pr(Nx)
      real*8 fd(Nx),fsr(Nx)
      real*8 ftau(Nx)
      real*8 dalpha, dbeta, sg, sr_sav, tau_sav
      real*8 dlogpsi, rhohloc, ploc, utloc
      real*8 rholoc, vloc, uloc, tloc, csloc, wloc
      real*8 lorentz_max, factor, vmp, vmm
      real*8 alphaloc, psiloc, betaloc, psi4loc
      real*8 sr_metric_source, tau_metric_source
      real*8 dloc, srloc, tauloc, xloc, jrloc, rhopsloc, ham_sourceloc
      integer i, limiter
      real*8 dx
      integer il

      dx=(x(2)-x(1))

      save_fcs = 0
      if (iter==2) save_fcs = 1
      
c     Use EOS to calculate pressure
      if (eos_flag == 0) then
         do i=1,Nx
            P(i) = (gamma-1.d0)*u(i)
         end do
      end if

      do i =1,Nx
c     zero fluxes
         fd(i)   = 0.d0
         fsr(i)  = 0.d0
         ftau(i) = 0.d0
      end do

c     Do primitive variable inversion.
      if (eos_flag == 0) then
         call find_primitives(d,sr,tau,rho,v,u,jr,rhops,
     &        ham_source,gamma,alpha_c,beta_c,psi_c,x_c,Nx,rho_atm,
     &        lorentz_max) 
      else if (eos_flag == 1) then
         call find_primitives_table(d,sr,tau,rho,v,u,jr,rhops,
     &        ham_source,P,T,cs,Ye,rho0,alpha_c,beta_c,psi_c,x_c,Nx,
     &        deltat1,deltat2,rho_atm,p_atm,u_atm,t_atm,lorentz_max)
      end if         
      
c     Do some outer boundary conditions.
      if (phys_bdy(2) .ne. 0) then
         rholoc   = rho(Nx-1)
         vloc     = max(v(Nx-1),0.d0)
         uloc     = u(Nx-1)
         ploc     = P(Nx-1)
         tloc     = T(Nx-1)
         csloc    = cs(Nx-1)
         xloc     = x_c(Nx)
         alphaloc = alpha_c(Nx)
         betaloc  = beta_c(Nx)
         psiloc   = psi_c(Nx)
c     recompute the conservative variables there.
         call compute_conservatives(dloc,tauloc,srloc,rholoc,ploc,
     &        uloc,vloc,alphaloc,betaloc,psiloc,xloc,gamma,eos_flag,
     &        lorentz_max)
         sg = xloc**2*psiloc**6
         wloc = dloc / (sg*rholoc)
         jrloc = srloc/sg
         rhopsloc = (rholoc+ploc+uloc)*(-1.d0+2.d0*wloc**2)+2.d0*ploc
         ham_sourceloc = (rholoc+ploc+uloc)*wloc**2-ploc
         
c     assign outgoing values
         rho(Nx)        = rholoc
         v(Nx)          = vloc
         u(Nx)          = uloc
         P(Nx)          = ploc
         T(Nx)          = tloc
         cs(Nx)         = csloc
         d(Nx)          = dloc
         tau(Nx)        = tauloc
         sr(Nx)         = srloc
         jr(Nx)         = jrloc
         rhops(Nx)      = rhopsloc
         ham_source(Nx) = ham_sourceloc
      end if

      if (iter < 3) then
c     call subroutine to advect in x-direction
         call hydro_advection(limiter,rho,v,P,u,cs, 
     &        fd,fsr,ftau,gamma,dx,Nx,x,x_c,eos_flag,
     &        alpha,beta,psi,phys_bdy,lorentz_max)
         
c     Finite differencing!
c     Assume zero flux through outer boundaries
         do i=1,Nx
            if (i==Nx) then
               d_p(i)   = 1.d0/dx * fd(i)
               sr_p(i)  = 1.d0/dx * fsr(i)
               tau_p(i) = 1.d0/dx * ftau(i)
            else
               d_p(i)   = - 1.d0/dx * (fd(i+1)-fd(i)) 
               sr_p(i)  = - 1.d0/dx * (fsr(i+1)-fsr(i)) 
               tau_p(i) = - 1.d0/dx * (ftau(i+1)-ftau(i)) 
            end if
c     Save flux corrections if necessary.
c     We first cast the real value of the mask to an integer.
            lfc_mask = fc_mask(i)
            lsign = -1.0
c     if (iand(lfc_mask,SIGN_FLAG) .ne. 0) lsign = 1.0
c     Note that lsign = 1 corresponds to the
c     natural sign, i.e., the sign with which
c     the fluxes contribute for evolution.
c     This natural sign should be used for type 
c     B cells.  
            if ((iand(lfc_mask,PXFLAG).ne.0) .and. 
     &           save_fcs==1) then
c     saving corrections at + x
               d_fcs(i)   = d_fcs(i)  - 
     &              lsign * dt/dx * fd(i+1)
               sr_fcs(i)  = sr_fcs(i) - 
     &              lsign * dt/dx * fsr(i+1)
               tau_fcs(i) = tau_fcs(i)  - 
     &              lsign * dt/dx * ftau(i+1)
            end if
            if ((iand(lfc_mask,MXFLAG).ne.0) .and. 
     &           save_fcs==1) then
c     saving corrections at - x
               d_fcs(i)   = d_fcs(i)  +
     &              lsign * dt/dx * fd(i)
               sr_fcs(i)  = sr_fcs(i) + 
     &              lsign * dt/dx * fsr(i)
               tau_fcs(i) = tau_fcs(i)  + 
     &              lsign * dt/dx * ftau(i)
            end if
         
c     Add metric source terms to the hydro evolution equations.
c     Note that these derivatives are cell-centered.
            dalpha  = (alpha(i+1)    - alpha(i))   /dx
            dbeta   = (beta(i+1)     - beta(i))    /dx
            
            dlogpsi = (log(psi(i+1)) - log(psi(i)))/dx
            sg = x_c(i)**2*alpha_c(i)*psi_c(i)**6
            psiloc = psi_c(i)
            
            ploc = P(i)
            rhohloc = rho(i) + u(i) + P(i)
            alphaloc = alpha_c(i)
            psi4loc = psiloc**4
            betaloc = beta_c(i)
            
            utloc  = -alphaloc**2 + psi4loc*(betaloc**2 + 
     &           2.d0*betaloc*v(i) + v(i)**2)
            if (utloc .ge. 0.d0) then
               write(*,*) 'superluminality inside hydro_1step at i = ', 
     &              i
               factor = (1.d0-1.d0/lorentz_max**2)
               factor = sqrt(factor)
c     Find the maximum allowed ingoing and outgoing
c     velocities--for our given speed limit.
               vmp    = alphaloc/psiloc**2*factor - betaloc
               vmm    = -alphaloc/psiloc**2*factor - betaloc
c     Choose the one closest to our guess.
               factor = abs(vmp - v(i))/abs(vmm - v(i))
               if (factor > 1) then
                  v(i) = vmm
               else
                  v(i) = vmp
               end if
               utloc = lorentz_max/alphaloc
            else
               utloc = 1.d0/sqrt(-utloc)
            end if
            
            sr_metric_source = sr(i)*dbeta - tau(i)*dalpha +
     &           2.d0*sg*( (rhohloc*(-1.d0+(alphaloc*utloc)**2) 
     &           + 3.d0*ploc)*dlogpsi + ploc/x_c(i) )
            
            tau_metric_source = 2.d0/3.d0*x_c(i)**2*psiloc**6*rhohloc * 
     &           (-1.d0+(alphaloc*utloc)**2)*(dbeta - betaloc/x_c(i))
     &           - sr(i)*dalpha/psiloc**4
            
            sr_p(i)  = sr_p(i)  + sr_metric_source
            tau_p(i) = tau_p(i) + tau_metric_source
         end do
      end if

      return
      end

c-----------------------------------------------------------------------
c     Do one psi step
c-----------------------------------------------------------------------
      subroutine psi_1step(psi,psi_p,beta,x,Nx,phys_bdy,iter)
      implicit none
      integer Nx
      real*8 psi(Nx+1),psi_p(Nx+1),beta(Nx+1),x(Nx+1)
      integer i, iter
      real*8 twodx,psip,psil,dx, betap
      integer phys_bdy(6)
      integer il

      dx = x(2)-x(1)
      twodx=2.d0*dx

c     Neumann condition at zero
      if (phys_bdy(1) .ne. 0) psi(1) = (4.d0*psi(2)-
     &     psi(3))/3.d0
      
c     Falloff condtion at outer boundary
c     Warning!  This is not right.  You should
c     Interpolate along a characteristic.
      if (phys_bdy(2) .ne. 0) psi(Nx+1) = 3.d0*psi(Nx) -
     &     3.d0*psi(Nx-1)+psi(Nx-2)

      if (iter < 3) then
         do i = 2, Nx
            psip = (psi(i+1)-psi(i-1))/twodx
            psil = psi(i)
            betap = (beta(i+1)-beta(i-1))/twodx
            psi_p(i) = beta(i)*(psil/(3*x(i))+psip)+psil*betap/6.d0
         end do
      end if
      
      return
      end

c-------------------------------------------------------------------------------------
c     This is the subroutine which calculates fluxes.
c-------------------------------------------------------------------------------------
      subroutine hydro_advection(limiter,rho,v,P,u,cs, 
     &     fd,fsr,ftau,gamma,dx,Nx,x,x_c,eos_flag, 
     &     alpha,beta,psi,phys_bdy,lorentz_max)
      implicit none
      integer Nx
      real*8 rho(Nx), P(Nx), u(Nx)
      real*8 v(Nx), cs(Nx)
      real*8 fd(Nx),fsr(Nx), ftau(Nx)
      real*8 x(Nx+1), x_c(Nx)
      real*8 alpha(Nx+1),beta(Nx+1),psi(Nx+1)
      real*8 gamma, dx
      integer eos_flag
      integer phys_bdy(6)
      
c     Local variables.
      real*8 rhol(Nx),rhor(Nx)
      real*8 vl(Nx),vr(Nx)
      real*8 Pl(Nx),Pr(Nx)
      real*8 ul(Nx),ur(Nx)
      real*8 csl(Nx),csr(Nx)
      
      real*8 rhol_l, rhor_l, Pl_l, Pr_l
      real*8 vl_l, vr_l
      real*8 fdl, fdr, fsrl, fsrr
      real*8 ftaul, ftaur, vil_l, vir_l
      real*8 cmax_l, cmin_l
      real*8 u_tl, u_rl
      real*8 u_tr, u_rr
      real*8 utl, url
      real*8 utr, urr
      real*8 hl_l,hr_l, ul_l, ur_l
      real*8 cpr,cpl,cmr,cml,csl_l,csr_l
      real*8 dl_l, taul_l, srl_l
      real*8 dr_l, taur_l, srr_l
      real*8 alpha_l, beta_l
      real*8 psi_l, psi_r, psi4_l, psi4_r
      
      real*8 sg
      integer i
      integer dir_flag, limiter
      real*8 lorentz_max, factor, vmp, vmm
      real*8 SYM, ANTI
      parameter (SYM=1.d0, ANTI=-1.d0)

c     Calculate left and right states
      call find_ur_ul(limiter,rho,rhol,rhor,dx,Nx,SYM,phys_bdy)
      call find_ur_ul(limiter,P,Pl,Pr,dx,Nx,SYM,phys_bdy)
      call find_ur_ul(limiter,v,vl,vr,dx,Nx,ANTI,phys_bdy)
      if (eos_flag == 1) then
         call find_ur_ul(limiter,u,ul,ur,dx,Nx,SYM,phys_bdy)
         call find_ur_ul(limiter,cs,csl,csr,dx,Nx,SYM,phys_bdy)
      end if

c     calculate and return fluxes
c     For now, simply apply the HLL formula.
      do i=2,Nx
         rhol_l  = rhol(i)
         rhor_l  = rhor(i)
         vl_l    = vl(i)
         vr_l    = vr(i)
         Pl_l    = Pl(i)
         Pr_l    = Pr(i)
         psi_l   = psi(i)
         alpha_l = alpha(i)
         beta_l  = beta(i)
         psi4_l  = psi_l**4
         sg      = x(i)**2*psi4_l*psi_l**2*alpha_l

         if (eos_flag == 0) then
            hl_l = 1.d0 + gamma * Pl_l / (gamma-1.d0) / rhol_l                 
            hr_l = 1.d0 + gamma * Pr_l / (gamma-1.d0) / rhor_l
         else if (eos_flag == 1) then 
            ul_l = ul(i)
            ur_l = ur(i)
            hl_l = 1.d0 + (Pl_l + ul_l) / rhol_l                 
            hr_l = 1.d0 + (Pr_l + ur_l) / rhor_l
         end if

c
c     Calculate u_\mu,L locally 
c
         utl = -alpha_l**2 + psi4_l*(beta_l**2 + 
     &        2.d0*beta_l*vl_l + vl_l**2)
         if (utl .ge. 0.d0) then
            write(*,*) 'superluminality on the left! i = ', i
            factor = (1.d0-1.d0/lorentz_max**2)
            factor = sqrt(factor)
c     Find the maximum allowed ingoing and outgoing
c     velocities--for our given speed limit.
            vmp    = alpha_l/psi_l**2*factor - beta_l
            vmm    = -alpha_l/psi_l**2*factor - beta_l
c     Choose the one closest to our guess.
            factor = abs(vmp - vl_l)/abs(vmm - vl_l)
            if (factor > 1) then
               vl_l = vmm
            else
               vl_l = vmp
            end if
            utl = lorentz_max/alpha_l
         else               
            utl = 1.d0/sqrt(-utl)
         end if
         url = utl*vl_l
         u_rl = utl*psi4_l*(vl_l+beta_l)
c     conservatives on left side.
         dl_l   = sg*rhol_l*utl
         taul_l = sg*alpha_l*(rhol_l*hl_l*utl**2-Pl_l/alpha_l**2)
         srl_l  = sg*rhol_l*hl_l*utl*u_rl

c         
c     Calculate u_\mu,R locally 
c
         utr = -alpha_l**2 + psi4_l*(beta_l**2 + 
     &        2.d0*beta_l*vr_l + vr_l**2)
         if (utr .ge. 0.d0) then
            write(*,*) 'superluminality on the right! i = ', i
            factor = (1.d0-1.d0/lorentz_max**2)
            factor = sqrt(factor)
            vmp    = alpha_l/psi_l**2*factor - beta_l
            vmm    = -alpha_l/psi_l**2*factor - beta_l
            factor = abs(vmp - vr_l)/abs(vmm - vr_l)
            if (factor > 1) then
               vr_l = vmm
            else
               vr_l = vmp
            end if
            utr = lorentz_max/alpha_l
         else
            utr = 1.d0/sqrt(-utr)
         end if
         urr = utr*vr_l
         u_rr = utr*psi4_l*(vr_l+beta_l)
c     conservatives on right side.
         dr_l   = sg*rhor_l*utr
         taur_l = sg*alpha_l*(rhor_l*hr_l*utr**2-Pr_l/alpha_l**2)
         srr_l  = sg*rhor_l*hr_l*utr*u_rr

c     Calculate fluxes on the right and left sides of the interface.
         fdl   = sg*rhol_l*utl*vl_l
         fdr   = sg*rhor_l*utr*vr_l
         ftaul = fdl*hl_l*utl*alpha_l + sg*beta_l*Pl_l/alpha_l
         ftaur = fdr*hr_l*utr*alpha_l + sg*beta_l*Pr_l/alpha_l
         fsrl  = fdl*hl_l*u_rl + sg*Pl_l
         fsrr  = fdr*hr_l*u_rr + sg*Pr_l
               
c     Calculate local comoving sound speed.
         if (eos_flag == 0) then
            csl_l = sqrt(gamma*Pl_l/(rhol_l*hl_l))
            csr_l = sqrt(gamma*Pr_l/(rhor_l*hr_l))
         else if (eos_flag == 1) then
            csl_l = csl(i)
            csr_l = csr(i)
         end if
               
c     Calculate c_+L, c_+R, c_-L, and c_-R
         call find_cp_cm(cpl, cml, csl_l, alpha_l, beta_l, psi_l,
     &        utl, url)
         call find_cp_cm(cpr, cmr, csr_l, alpha_l, beta_l, psi_l,
     &        utr, urr)
         cmax_l = max(0.d0,cpr,cpl)
         cmin_l = -min(0.d0,cmr,cml)
         
c     Apply HLL formula
         fd(i) = cmin_l * fdr + cmax_l * fdl - 
     &        cmax_l*cmin_l*(dr_l - dl_l)
         fd(i) = fd(i) / (cmax_l + cmin_l)

         fsr(i) = cmin_l * fsrr + cmax_l * fsrl - 
     &        cmax_l*cmin_l*(srr_l - srl_l)
         fsr(i) = fsr(i) / (cmax_l + cmin_l)
         
         ftau(i) =  cmin_l * ftaur + cmax_l * ftaul - 
     &        cmax_l*cmin_l*(taur_l - taul_l)
         ftau(i) = ftau(i) / (cmax_l + cmin_l)
         
      end do
      return
      end

c-------------------------------------------------------------
c      Reconstruction step.  Note:
c        limiter = 0 for zero slope
c        limiter = 1 for MC limiter
c-------------------------------------------------------------
      subroutine find_ur_ul(limiter,f,fl,fr,dx,Nx,sym,phys_bdy)
      implicit none
      integer Nx,i
      real*8 f(Nx), fr(Nx)
      real*8 fl(Nx), f_slopes(Nx)
      real*8 dx, dxo2
      real*8 sym
      integer phys_bdy(6)
      integer limiter            
      real*8 xa(4), ya(4), interpolant, xloc, dy

c     The mc_slopes routine doesn't
c     divide by dx, so you shouldn't multiply by it here.
      dxo2 = 1.d0 / 2.d0
      
      if (limiter == 0) then
         do i = 1, Nx
            f_slopes(i) = 0.d0
         end do
      else if (limiter == 1) then
         call mc_slopes(f,f_slopes,Nx,sym,phys_bdy)
      end if
      
      do i = 1, Nx
         if (i==1) then
            fr(i) = f(i)   - dxo2 * f_slopes(i)
            fl(i) = sym*fr(i)  
         else
            fr(i) = f(i)   - dxo2 * f_slopes(i)
            fl(i) = f(i-1) + dxo2 * f_slopes(i-1)
         end if
      end do
      
      ! Make the r and l states smooth near the center.
c$$$      do i = 1, 4
c$$$         xloc = (i-1)*dx
c$$$
c$$$         if (i==1 .or. i==2 .or. i==3) then
c$$$            xa(1) = 0.5*dx
c$$$            xa(2) = 1.5*dx
c$$$            xa(3) = 2.5*dx
c$$$            xa(4) = 3.5*dx
c$$$            ya(1) = f(1)
c$$$            ya(2) = f(2)
c$$$            ya(3) = f(3)
c$$$            ya(4) = f(4)
c$$$         else
c$$$            xa(1) = 1.5*dx
c$$$            xa(2) = 2.5*dx
c$$$            xa(3) = 3.5*dx
c$$$            xa(4) = 4.5*dx
c$$$            ya(1) = f(2)
c$$$            ya(2) = f(3)
c$$$            ya(3) = f(4)
c$$$            ya(4) = f(5)
c$$$         end if            
c$$$         
c$$$         call polint(xa,ya,4,xloc,interpolant,dy)
c$$$         fr(i) = interpolant
c$$$         fl(i) = interpolant
c$$$      end do

      return
      end 

c-----------------------------------------------------------------------
c     Initial data with polytropic or Shen EOS
c-----------------------------------------------------------------------
      subroutine id_tov(gamma,rho0_c,atm_frac,dn,srn,taun,rho,v,u,
     &     P,T,jr,jr_v,rhops,rhops_v,tau_v,ham_source,ham_source_v,
     &     alpha,beta,psi,alpha_c,
     &     beta_c,psi_c,x,x_c,Nx,ham,mask,phys_bdy,ghost_width,
     &     p_deplete,Ye,temp,rho0_scale,eos_flag,alpha_out,rho_atm,
     &     p_atm,u_atm,t_atm,deltat1,deltat2,v_pert,lorentz_max,
     &     psi0_out)
      implicit none
      real*8 gamma, rho0_c, atm_frac
      integer Nx
      real*8 dn(Nx), srn(Nx), taun(Nx), rho(Nx), v(Nx)
      real*8 rhops(Nx), jr(Nx), ham_source(Nx)
      real*8 u(Nx), P(Nx), T(Nx), alpha(Nx+1), beta(Nx+1), psi(Nx+1)
      real*8 alpha_c(Nx), beta_c(Nx), psi_c(Nx)
      real*8 x(Nx+1), x_c(Nx)
      real*8 jr_v(Nx+1), rhops_v(Nx+1), tau_v(Nx+1), ham_source_v(Nx+1)
      real*8 ham(Nx+1), mask(Nx+1)
      real*8 Ye, temp, rho0_scale, deltat1,deltat2, psi0_out
      integer phys_bdy(6)
      integer ghost_width(6)
      real*8 p_deplete, alpha_out, rho_atm, p_atm, u_atm, t_atm
      real*8 v_pert
      
c     Local variables.
      integer NRS, NRS_ex
      real*8 r_s(12*Nx+1), P_s(12*Nx+1), rho0_s(12*Nx+1)
      real*8 phi_s(12*Nx+1), m_s(12*Nx+1), r_iso_s(12*Nx+1)
      real*8 u_s(12*Nx+1), T_s(12*Nx+1)
      real*8 alpha_s(12*Nx+1), psi_s(12*Nx+1)
      real*8 r_s_ex(36*Nx+1), r_iso_s_ex(36*Nx+1), m_s_ex(36*Nx+1)
      real*8 rho_c, P_c, eps_c, phi_c, x_last
      real*8 r_iso_c
      integer n, n2, i, j, kill_flag, i_last, jup
      integer j1, j2, j3, j4, n_int, iloc,jloc,kloc
      parameter (n=3)
      parameter (n2=1)
      parameter (n_int=4)
      external derivs, derivs2
      real*8 y(n), dydx(n), h, yout(n), xloc
      real*8 y2(n2), dydx2(n2), yout2(n2)
      real*8 xa(n_int), ya(n_int), dy
      real*8 mtot, rsmax, rsmax_ex
      real*8 rloc, interpolant, ploc
      real*8 rholoc, uloc, vloc, alphaloc, betaloc, psiloc
      real*8 dloc, tauloc, srloc, cs_c, epsloc, csloc, tloc, mloc
      real*8 dphi,d2phi
      logical reached_surface, cap, floor
      integer eos_flag, istart
      real*8 lorentz_max
      real*8 w, rhoh, ham_norm, radius, L_scale, radius_iso
      real*8 G_newt, c
      integer il,jl,kl
      parameter (G_newt = 6.67259d-8)
      parameter (c = 2.99792458d10)
      character(len=30) out_name
      integer write_profile

      L_scale = sqrt(c**2/(G_newt*rho0_scale))
      
c     Decide on a grid for your tov solution in Schwarzshild 
c     coordinates
      NRS = 12*Nx+1
      NRS_ex = 36*Nx+1
      rsmax = x(Nx+1)
      rsmax_ex = (NRS_ex-1)/(NRS-1)*rsmax

      if (eos_flag==0) then
         P_c     = rho0_c**gamma
         eps_c   = P_c/((gamma-1.d0)*rho0_c)
      else 
c     You need to call the EOS for a known temperature here.
         il = 0
         jl = 0
         kl = 0
         call eos_wrapper_temp(P_c,rho0_c,eps_c,Ye,temp,rho0_scale,cs_c,
     &        il,jl,kl)
      end if
      rho_c   = rho0_c*(1.d0+eps_c)
      phi_c   = -1.d0
      r_iso_c = 0.d0
      p_atm   = atm_frac * P_c
      if (eos_flag == 0) then
         rho_atm = p_atm**(1.d0/gamma)
         u_atm = p_atm / rho_atm
         t_atm = 0.d0
      end if

      r_s(1)     = 0.d0
      P_s(1)     = P_c
      u_s(1)     = rho0_c*eps_c
      rho0_s(1)  = rho0_c
      phi_s(1)   = phi_c
      m_s(1)     = 0.d0
      if (eos_flag == 0) then
         T_s(1) = 0.d0
      else
         T_s(1)     = temp
      end if

c     Initialize
      y(1) = m_s(1)
      y(2) = P_s(1)
      y(3) = phi_s(1)
      dydx(1) = 0.d0
      dydx(2) = 0.d0
      dydx(3) = 0.d0
      reached_surface = .false.

c     set up radial grid.
      h = rsmax/(NRS-1)
      do i = 1, NRS-1
         r_s(i+1) = r_s(i) + h
      end do

      i_last = NRS
c     Integrate until you reach P=0.
      do i = 1, NRS-1
         xloc = r_s(i)
         iloc = i+1
         if (.not. reached_surface) then
            kill_flag = 0
            rholoc = rho0_s(i)
            call rk4(y,dydx,n,xloc,h,yout,derivs,kill_flag,gamma,
     &           Ye,temp,rholoc,rho0_scale,iloc,eos_flag,p_atm)
            
            if (kill_flag .ne. 0 .or. yout(2) < p_atm) then
               reached_surface = .true.
               mtot = m_s(i)
               P_s(i+1)  = p_atm
               m_s(i+1)  = mtot
               phi_s(i+1) = 0.5d0*log(1.d0-2.d0*mtot/r_s(i+1))
               radius = r_s(i)
               i_last = i

               if (eos_flag == 0) then
                  write(*,*) 'mtot = ', mtot
                  write(*,*) 'areal radius = ', radius
               else
c     Write out the mass in solar masses.
                  write(*,*) 'mtot = ', mtot*rho0_scale*
     &                 L_scale**3/(2.d33)
                  write(*,*) 'unscaled mtot = ', mtot
                  write(*,*) 'unscaled areal radius = ', radius
               end if
            else
               m_s(i+1)     = yout(1)
               P_s(i+1)     = yout(2)
               phi_s(i+1)   = yout(3)

               y(1) = m_s(i+1)
               y(2) = P_s(i+1)
               y(3) = phi_s(i+1)

               call derivs(r_s(i+1),y,dydx,kill_flag,gamma,Ye,temp,
     &              rholoc,rho0_scale,iloc,eos_flag,p_atm)
               if (kill_flag .ne. 0) then
                  reached_surface = .true.
                  mtot = m_s(i+1)
                  i_last = i+1
                  radius = r_s(i+1)
               end if
            end if
         else
            P_s(i+1)   = p_atm
            m_s(i+1)   = mtot
            phi_s(i+1) = 0.5d0*log(1.d0-2.d0*mtot/r_s(i+1))
         end if
         
         if (eos_flag == 0) then
            rho0_s(i+1) = (P_s(i+1))**(1.d0/gamma)
            u_s(i+1)    = P_s(i+1)/(gamma-1.d0)
            T_s(i+1)    = 0.d0
         else
            ploc = P_s(i+1)
            jloc = 1
            kloc = 1
c     Use previous value of rho as an initial guess.
            call finddens_wrapper(ploc,rholoc,epsloc,Ye,temp,
     &           rho0_scale,csloc,iloc,jloc,kloc)
            rho0_s(i+1) = rholoc
            u_s(i+1)  = epsloc*rholoc
            T_s(i+1)  = temp
         end if
      end do

c     Set the atmosphere density in the Shen eos case.
c     This assumes that the last point of your TOV integration 
c     lies in the atmosphere.  If that's not the case, then 
c     you're probably not too concerned about having an atmosphere.
      if (eos_flag==1) then
         rho_atm = rho0_s(NRS)
         u_atm = u_s(NRS)
         t_atm = temp
      end if

      ! Shift phi so that you match onto Schwarzshild exterior.
      x_last = r_s(i_last)
      do i=1,i_last
         phi_s(i) = phi_s(i) + 0.5d0*log(1.d0-2.d0*mtot/x_last) - 
     &        phi_s(i_last)
      end do

c     Determine isotropic radius by integrating from large distance.
c     To do this properly, we set up extended grids.
c     set up radial grid.
      h = rsmax/(NRS-1)
      r_s_ex(1) = 0.d0
      do i = 1, NRS_ex-1
         r_s_ex(i+1) = r_s_ex(i) + h
      end do

      do i = 1, NRS_ex
         if (i<i_last) then
            m_s_ex(i) = m_s(i)
         else
            m_s_ex(i) = mtot
         end if
      end do
            
      y2(1) = r_s_ex(NRS_ex)
      dydx2(1) = 1.d0
      r_iso_s_ex(NRS_ex) = r_s_ex(NRS_ex)
      xloc = r_s_ex(NRS_ex)
      mloc = m_s_ex(NRS_ex)
      
      do i = NRS_ex, 2, -1
         call rk4_1real(y2,dydx2,n2,xloc,-h,yout2,m_s_ex,r_s_ex,
     &        NRS_ex,derivs2)
         r_iso_s_ex(i-1) = yout2(1)
         y2(1) = yout2(1)
         xloc = r_s_ex(i-1)
         call derivs2(xloc,y2,dydx2,m_s_ex,r_s_ex,NRS_ex)
      end do

c     Now fill up the shorter array
      do i = 1, NRS
         r_iso_s(i) = r_iso_s_ex(i)
      end do

c     Calculate conformal factor and lapse.
      alpha_s(1) = exp(phi_s(1))
      do i = 2, NRS
         alpha_s(i) = exp(phi_s(i))
         psi_s(i) = sqrt(r_s(i)/r_iso_s(i))
      end do

c     First derivative of psi vanishes at origin, so, to second order,
      psi_s(1) = (4.d0*psi_s(2)-psi_s(3))/3.d0
c      psi_s(1) = psi_s(2)

c     Calculate the radius in isotropic coordinates
      radius_iso = radius * (psi_s(i_last))**2

c     Before interpolating onto the numerical grid,
c     multiply by the pressure depletion factor
      if (p_deplete .ne. 1.0) then
         do i = 1, NRS
            if (eos_flag==0) then
               P_s(i) = P_s(i) * p_deplete
               u_s(i) = u_s(i) * p_deplete
            else
               if (P_s(i) > p_atm) then
                  ploc = P_s(i)*p_deplete
                  rholoc = rho0_s(i)
                  epsloc = u_s(i)/rholoc
                  tloc = temp
                  iloc = i
                  jloc = 0
                  kloc = 0
                  call eos_wrapper_press(ploc,rholoc,epsloc,Ye,tloc,
     &                 rho0_scale,csloc,iloc,jloc,kloc,cap,floor,
     &              deltat1,deltat2)
c     Assign outgoing values.
                  P_s(i) = ploc
                  u_s(i) = epsloc*rholoc
                  T_s(i) = tloc
               end if
            end if
         end do
      end if

c     Interpolate to a grid in isotropic radius.
c     Do the cell centered stuff first.
      do i = 1, Nx
         rloc = x_c(i)
         
c     Find surrouding points for the polynomial interp
         do j = 1, NRS
            if (r_iso_s(j) .gt. rloc) then
               jup = j
               exit
            end if
         end do

         if (jup == 1 .or. jup == 2) then
            j1 = 1
            j2 = 2
            j3 = 3
            j4 = 4
         else if (jup == NRS .or. jup == NRS - 1) then
            j1 = NRS-3
            j2 = NRS-2
            j3 = NRS-1
            j4 = NRS
         else
            j1 = jup - 2
            j2 = jup - 1
            j3 = jup 
            j4 = jup + 1
         end if
      
         xa(1) = r_iso_s(j1)
         xa(2) = r_iso_s(j2)
         xa(3) = r_iso_s(j3)
         xa(4) = r_iso_s(j4)

c     Interpolate rho0
         ya(1) = rho0_s(j1)
         ya(2) = rho0_s(j2)
         ya(3) = rho0_s(j3)
         ya(4) = rho0_s(j4)

         call polint(xa,ya,n_int,rloc,interpolant,dy)
         rho(i) = interpolant

c     Also, assign velocity
         if (rloc .le. radius_iso) then
            v(i)   = - v_pert * rloc / radius_iso
         else
            v(i)   = 0.d0
         end if

c     Interpolate P
         ya(1) = P_s(j1)
         ya(2) = P_s(j2)
         ya(3) = P_s(j3)
         ya(4) = P_s(j4)

         call polint(xa,ya,n_int,rloc,interpolant,dy)
         ploc = interpolant
         P(i) = interpolant

         if (eos_flag == 0) then
            u(i) = ploc/(gamma-1.d0)
            T(i) = 0.d0
         else
c     You need to interpolate u
            ya(1) = u_s(j1)
            ya(2) = u_s(j2)
            ya(3) = u_s(j3)
            ya(4) = u_s(j4)

            call polint(xa,ya,n_int,rloc,interpolant,dy)
            u(i) = interpolant

c     You also need to interpolate T in case you did pressure depletion.
            ya(1) = T_s(j1)
            ya(2) = T_s(j2)
            ya(3) = T_s(j3)
            ya(4) = T_s(j4)

            call polint(xa,ya,n_int,rloc,interpolant,dy)
            T(i) = interpolant
         end if

c     Interpolate alpha
         ya(1) = alpha_s(j1)
         ya(2) = alpha_s(j2)
         ya(3) = alpha_s(j3)
         ya(4) = alpha_s(j4)

         call polint(xa,ya,n_int,rloc,interpolant,dy)
         alpha_c(i) = interpolant
         
c     Assign shift
         beta_c(i)  = 0.d0

c     Interpolate psi
         ya(1) = psi_s(j1)
         ya(2) = psi_s(j2)
         ya(3) = psi_s(j3)
         ya(4) = psi_s(j4)

         call polint(xa,ya,n_int,rloc,interpolant,dy)
         psi_c(i) = interpolant

         rholoc   = rho(i)
         uloc     = u(i)
         vloc     = v(i)
         alphaloc = alpha_c(i)
         betaloc  = beta_c(i)
         psiloc   = psi_c(i)

         call compute_conservatives(dloc,tauloc,srloc,rholoc,ploc,
     &        uloc,vloc,alphaloc,betaloc,psiloc,rloc,gamma,
     &        eos_flag,lorentz_max)

         dn(i)   = dloc
         taun(i) = tauloc
         srn(i)  = srloc

c     Calculate sources for elliptics.
         rhoh = rholoc + uloc + ploc
         w = 1.d0
         jr(i) = srloc/(rloc**2*psiloc**6)
         rhops(i) = rhoh*(-1.d0+2.d0*w**2)+2.d0*ploc
         ham_source(i) = rhoh*w**2-ploc
      end do

      do i=1, Nx+1
         rloc = x(i)

c     Find surrouding points for the polynomial interp
         do j = 1, NRS
            if (r_iso_s(j) .ge. rloc) then
               jup = j
               exit
            end if
         end do

         if (jup == 1 .or. jup == 2) then
            j1 = 1
            j2 = 2
            j3 = 3
            j4 = 4
         else if (jup == NRS .or. jup == NRS - 1) then
            j1 = NRS-3
            j2 = NRS-2
            j3 = NRS-1
            j4 = NRS
         else
            j1 = jup - 2
            j2 = jup - 1
            j3 = jup 
            j4 = jup + 1
         end if
      
         xa(1) = r_iso_s(j1)
         xa(2) = r_iso_s(j2)
         xa(3) = r_iso_s(j3)
         xa(4) = r_iso_s(j4)

c     Interpolate rho0
         ya(1) = rho0_s(j1)
         ya(2) = rho0_s(j2)
         ya(3) = rho0_s(j3)
         ya(4) = rho0_s(j4)

         call polint(xa,ya,n_int,rloc,interpolant,dy)
         rholoc = interpolant

c     Also, assign velocity
         if (rloc .le. radius_iso) then
            vloc = -v_pert*rloc/radius_iso
         else
            vloc = 0.d0
         end if

c     Interpolate P
         ya(1) = P_s(j1)
         ya(2) = P_s(j2)
         ya(3) = P_s(j3)
         ya(4) = P_s(j4)

         call polint(xa,ya,n_int,rloc,interpolant,dy)
         ploc = interpolant

         if (eos_flag == 0) then
            uloc = interpolant/(gamma-1.d0)
         else
c     Interpolate u
            ya(1) = u_s(j1)
            ya(2) = u_s(j2)
            ya(3) = u_s(j3)
            ya(4) = u_s(j4)
            
            call polint(xa,ya,n_int,rloc,interpolant,dy)
            uloc = interpolant
         end if

c     Interpolate alpha
         ya(1) = alpha_s(j1)
         ya(2) = alpha_s(j2)
         ya(3) = alpha_s(j3)
         ya(4) = alpha_s(j4)

         call polint(xa,ya,n_int,rloc,interpolant,dy)
         alpha(i) = interpolant
         
c     Assign shift
         beta(i)  = 0.d0

c     Interpolate psi
         ya(1) = psi_s(j1)
         ya(2) = psi_s(j2)
         ya(3) = psi_s(j3)
         ya(4) = psi_s(j4)

         call polint(xa,ya,n_int,rloc,interpolant,dy)
         psi(i) = interpolant
         
         alphaloc = alpha(i)
         betaloc  = beta(i)
         psiloc   = psi(i)
         
         call compute_conservatives(dloc,tauloc,srloc,rholoc,ploc,
     &        uloc,vloc,alphaloc,betaloc,psiloc,rloc,gamma,
     &        eos_flag,lorentz_max)

c     Calculate sources on vertices!  
         rhoh = rholoc + uloc + ploc
         if (rloc > 0.d0) then
            jr_v(i) = srloc/(rloc**2*psiloc**6)
            w = dloc/(rloc**2*psiloc**6*rholoc)
         else
            jr_v(i) = 0.d0
            w = 1.d0
         end if
         rhops_v(i) = rhoh*(-1.d0+2.d0*w**2)+2.d0*ploc
         tau_v(i) = tauloc
         ham_source_v(i) = rhoh*w**2 - ploc
c     Also initialize the mask
         mask(i) = 1.d0

      end do

c     Set the value of the lapse at the outer boundary.  This will be 
c     used for a Dirichlet boundary condition on the lapse.
      alpha_out = alpha(Nx+1)
      psi0_out = psi(Nx+1)

c     Evaluate the hamiltonian constraint and write the norm.
      write_profile = 0;
      out_name = 'id_ham.dat'
      call ham_const(ham,alpha,beta,psi,tau_v,ham_source_v,mask,
     &     phys_bdy,x,ham_norm,Nx,ghost_width,out_name,write_profile)

      write(*,*) 'inside id_tov, ham_norm = ', ham_norm

      end subroutine

c-----------------------------------------------------------------------
c     Find primitive variables.
c
c     Metric source terms:  Definitions
c     j_r \equiv -\gamma^{\mu}_i n^{\nu} T_{\mu\nu}
c     rhops \equiv \rho + S
c        where
c     rho = n_{\mu}n_{\nu}T^{\mu\nu}
c     S = g^{\mu\nu}S_{\mu\nu}
c     S_{\mu\nu} = \gamma^{\rho}_{\mu} \gamma^{\sigma}_{\nu}T_{\rho\sigma}
c-----------------------------------------------------------------------
      subroutine find_primitives(dn,srn,taun,rho,v,u,jr,rhops,
     &     ham_source,gamma,alpha_c,beta_c,psi_c,x_c,Nx,rho_atm,
     &     lorentz_max) 
      implicit none
      integer Nx
      real*8 dn(Nx),srn(Nx),taun(Nx)
      real*8 alpha_c(Nx),beta_c(Nx),psi_c(Nx), x_c(Nx)
      real*8 v(Nx),u(Nx),rho(Nx), jr(Nx), rhops(Nx), ham_source(Nx)
      real*8 gamma,rho_atm
      real*8 P(Nx)

      real*8 g_tt,   g_tx,   g_ty,   g_tz,   g_xx
      real*8 g_xy,   g_xz,   g_yy,   g_yz,   g_zz
      real*8 gup_tt, gup_tx, gup_ty, gup_tz, gup_xx
      real*8 gup_xy, gup_xz, gup_yy, gup_yz, gup_zz
      real*8 detg
      real*8 rholoc,   ploc,   vxloc,  vyloc,  vzloc
      real*8 dloc, eloc, sxloc, syloc, szloc, tauloc
      real*8 alphaloc,psiloc,betaloc,psi4loc, uloc
      real*8 rhoh, w, sg, xloc, lorentz_max
      
      integer i,j,k, retval, eos_flag
      
c     If you're inside this subroutine, then eos_flag = 0
c     But you still need a legit space in memory so you 
c     can pass it to compute_conservatives.
      eos_flag = 0
      
      do i=1,Nx

         alphaloc = alpha_c(i)
         psiloc   = psi_c(i)
         betaloc  = beta_c(i)
         psi4loc  = psiloc**4
         xloc     = x_c(i)

         ! Set metric variables once and for all
         g_tt = -alphaloc**2 + psi4loc*betaloc**2
         g_xx = psi4loc
         g_yy = 1.d0
         g_zz = 1.d0
         
         g_tx = psi4loc*betaloc
         g_ty = 0.d0
         g_tz = 0.d0
         g_xy = 0.d0
         g_xz = 0.d0
         g_yz = 0.d0
         
         gup_tt = -1.d0/alphaloc**2
         gup_xx = 1.d0/psi4loc - betaloc**2/alphaloc**2
         gup_yy = 1.d0
         gup_zz = 1.d0
         
         gup_tx = betaloc/alphaloc**2
         gup_ty = 0.d0
         gup_tz = 0.d0
         gup_xy = 0.d0
         gup_xz = 0.d0
         gup_yz = 0.d0
         
         detg = alphaloc*psi4loc*psiloc**2*xloc**2
               
         dloc    = dn(i)
         eloc  = -alphaloc*taun(i)+betaloc*srn(i)
         sxloc = srn(i)
         syloc = 0.d0
         szloc = 0.d0
         
         rholoc = rho(i)
         ploc = u(i)*(gamma-1.d0)
         vxloc = v(i)
         vyloc = 0.d0
         vzloc = 0.d0

c     This is the wrapper for the Noble et al. solver.
         call inversion_interface(retval,
     &        rholoc,ploc,vxloc,vyloc,vzloc, 
     &        dloc,eloc,sxloc,syloc,szloc, 
     &        g_tt,g_tx,g_ty,g_tz,g_xx,g_xy,g_xz,g_yy,g_yz,g_zz, 
     &        gup_tt,gup_tx,gup_ty,gup_tz,gup_xx,gup_xy,gup_xz,
     &        gup_yy,gup_yz,gup_zz,detg,gamma)
               
         dn(i) = dloc
         taun(i) = (-eloc+betaloc*sxloc)/alphaloc
         srn(i) = sxloc
         
         rho(i) = rholoc
         u(i)   = ploc/(gamma-1.d0)
         v(i)   = vxloc

c     Impose density floor, leaving the velocity unchanged.
         if (rho(i) .lt. rho_atm) then
            rholoc = rho_atm
            ploc = rholoc**gamma
            uloc = ploc/(gamma-1.d0)

c     Recalculate conserved variables.
            call compute_conservatives(dloc,tauloc,sxloc,rholoc,ploc,
     &           uloc,vxloc,alphaloc,betaloc,psiloc,xloc,gamma,eos_flag,
     &           lorentz_max)
            rho(i) = rholoc
            u(i) = uloc
            dn(i) = dloc
            taun(i) = tauloc
            srn(i) = sxloc
         end if

c     Calculate the matter source terms for the metric.
         sg = detg/alphaloc
         rhoh = rholoc + u(i) + ploc
         w = dloc/(sg*rholoc)
         jr(i) = sxloc/sg
         rhops(i) = rhoh*(-1.d0+2.d0*w**2)+2.d0*ploc
         ham_source(i) = rhoh*w**2 - ploc
      end do
      
      return
      end

c-----------------------------------------------------------------------
c     Find primitive variables.  For a tabulated EOS.
c-----------------------------------------------------------------------
      subroutine find_primitives_table(dn,srn,taun,rho,v,u,jr,rhops,
     &     ham_source,P,T,cs,Ye,rho0,alpha_c,beta_c,psi_c,x,Nx,
     &     deltat1,deltat2,rho_atm,p_atm,u_atm,t_atm,lorentz_max)
      implicit none
      integer Nx
      real*8 dn(Nx),srn(Nx),taun(Nx)
      real*8 alpha_c(Nx),beta_c(Nx),psi_c(Nx)
      real*8 v(Nx),u(Nx),rho(Nx),T(Nx)
      real*8 jr(Nx),rhops(Nx), x(Nx), ham_source(Nx)
      real*8 P(Nx), cs(Nx)
      real*8 Ye, rho0,rho_atm
      real*8 rholoc,   ploc,   vloc
      real*8 dloc, tauloc, srloc, tloc, csloc, uloc
      real*8 p_atm, u_atm, t_atm
      real*8 lapseloc,shiftloc,psiloc,xloc
      real*8 gamma, lorentz_max
      real*8 sg, rhoh, w
      real*8 deltat1,deltat2
      integer eos_flag
      logical cap, floor
      integer i,retval

      eos_flag = 1
c     We won't be using gamma here.  If something goes wrong
c     and you do use it, let's set it to zero so that it'll throw
c     a nan somewhere.
      gamma = 0.d0

      do i=1,Nx
         dloc     = dn(i)
         tauloc   = taun(i)
         srloc    = srn(i)
         rholoc   = rho(i)
         ploc     = P(i)
         vloc     = v(i)
         tloc     = T(i)
         uloc     = u(i)
         xloc     = x(i)
         lapseloc = alpha_c(i)
         psiloc   = psi_c(i)
         shiftloc = beta_c(i)
         sg = psiloc**6*xloc**2
         cap = .false.
         floor = .false.

         call table_inversion_1d(retval,rholoc,ploc,uloc,tloc,
     &        vloc,csloc,dloc,tauloc,srloc,Ye,rho0,psiloc,lapseloc,
     &        shiftloc,xloc,i,cap,floor,deltat1,deltat2)
         
         if (cap .eqv. .true.) then
            write(*,*) 'warning: capping temperature at i = ', i
         end if

         if (floor .eqv. .true.) then
            write(*,*) 'warning: flooring temperature at i = ', i
         end if

c     Impose atmosphere density floors if necessary.
         if (rholoc < rho_atm) then
            rholoc = rho_atm
            tloc   = t_atm
            ploc   = p_atm
            uloc   = u_atm
c     Recalculate the conserved variables.
            call compute_conservatives(dloc,tauloc,srloc,rholoc,ploc,
     &           uloc,vloc,lapseloc,shiftloc,psiloc,xloc,gamma,eos_flag,
     &           lorentz_max)
         end if
        
c     Assign outgoing values.
         dn(i)   = dloc
         taun(i) = tauloc
         srn(i)  = srloc
         rho(i)  = rholoc
         u(i)    = uloc
         v(i)    = vloc
         T(i)    = tloc
         P(i)    = ploc
c     Achtung.  This cs will be wrong if you just floored it.
         cs(i)   = csloc

         ! Calculate the matter source terms for the metric.
         rhoh = rholoc + uloc + ploc
         w = dloc/(sg*rholoc)
         jr(i) = srloc/sg
         rhops(i) = rhoh*(-1.d0+2.d0*w**2)+2.d0*ploc
         ham_source(i) = rhoh*w**2 - ploc
      end do
      
      return
      end

c-------------------------------------------------------------
c     Compute slopes using the Monotonized-Centered limiter
c-------------------------------------------------------------
      subroutine mc_slopes(f, f_slope, Nx, sym, phys_bdy)
      implicit none
      integer Nx,i
      real*8 f(Nx), f_slope(Nx)
      real*8  temp
      real*8 dfp, dfm, dfc
      real*8 sym
      integer phys_bdy(6)
      
      do i = 2, Nx - 1
         dfp = f(i+1) - f(i)
         dfm = f(i) - f(i-1)
         dfc = 0.5d0 * (dfp + dfm)
         
         if (dfp*dfm < 0.d0) then 
            f_slope(i) = 0.d0
         else
            if (abs(2.d0*dfp) .lt. abs(2.d0*dfm)) then
               temp = 2.d0*dfp
            else
               temp = 2.d0*dfm
            end if
            
            if (abs(temp) .lt. abs(dfc)) then
               f_slope(i) = temp
            else
               f_slope(i) = dfc
            end if
         end if
      end do
c     Set slope values at boundaries. 
      
      if (phys_bdy(1) .ne. 0) then
         dfp = f(2) - f(1)
         dfm = f(1) - sym*f(1)
         dfc = 0.5d0 * (dfp + dfm)
         
         if (dfp*dfm < 0.d0) then 
            f_slope(1) = 0.d0
         else
            if (abs(2.d0*dfp) .lt. abs(2.d0*dfm)) then
               temp = 2.d0*dfp
            else
               temp = 2.d0*dfm
            end if
            
            if (abs(temp) .lt. abs(dfc)) then
               f_slope(1) = temp
            else
               f_slope(1) = dfc
            end if
         end if
      else
         f_slope(1)  = 0.d0
      end if
      f_slope(Nx) = 0.d0
      
      end subroutine mc_slopes

c-----------------------------------------------------------------
c     I'm not sure whether it's actually necessary to initialize
c     the advanced time level like this.
c-----------------------------------------------------------------
      subroutine init_rest(dn,dnp1,srn,srnp1,taun,taunp1,Nx)
      implicit none
      integer Nx
      real*8 dn(Nx),dnp1(Nx)
      real*8 srn(Nx),srnp1(Nx)
      real*8 taun(Nx),taunp1(Nx)
      integer i,j,k

      do i=1,Nx
         dnp1(i)   = dn(i)
         srnp1(i)  = srn(i)
         taunp1(i) = taun(i)
      end do
      
      return
      end

c-----------------------------------------------------------------
c     Subroutine for zeroing fcs variables.  The flag determines
c     whether type a or type b cells should be zeroed.     
c-----------------------------------------------------------------
      subroutine zero_fcs_vars(type_flag,fc_mask,
     &     d_fcs,sr_fcs,tau_fcs,Nx)
      implicit none
      integer Nx,type_flag
      integer fc_mask(Nx)
      integer lfc_mask
      real*8 d_fcs(Nx), tau_fcs(Nx), sr_fcs(Nx)
      include 'fc_mask.inc'
      integer i,j,k

      do i=1,Nx
         lfc_mask = fc_mask(i)
         if (iand(lfc_mask,SIGN_FLAG) == type_flag) then 
            if ((iand(lfc_mask,PXFLAG) .ne. 0) .or. 
     &           (iand(lfc_mask,MXFLAG) .ne. 0) .or. 
     &           (iand(lfc_mask,PYFLAG) .ne. 0) .or. 
     &           (iand(lfc_mask,MYFLAG) .ne. 0) .or. 
     &           (iand(lfc_mask,PZFLAG) .ne. 0) .or. 
     &           (iand(lfc_mask,MZFLAG) .ne. 0)) then 
               d_fcs(i)   = 0.d0
               sr_fcs(i)  = 0.d0
               tau_fcs(i) = 0.d0
            end if
         end if
      end do
      
      end subroutine

c-----------------------------------------------------------------
c     Subroutine for applying the flux correction.  This
c     probably could have been done in C.     
c-----------------------------------------------------------------
      subroutine apply_fc(dn,srn,taun,d_fcs,sr_fcs,tau_fcs,Nx)
      implicit none
      integer Nx
      real*8 d_fcs(Nx), tau_fcs(Nx), sr_fcs(Nx)
      real*8 dn(Nx), taun(Nx), srn(Nx)
      ! Locals
      integer i

      do i=1,Nx
         dn(i)   = dn(i)   + d_fcs(i)
         srn(i)  = srn(i)  + sr_fcs(i)
         taun(i) = taun(i) + tau_fcs(i)
         
      end do
      
      end subroutine

c--------------------------------------------------------------------
c     Compute the conserved (evolution) variables given the 
c     primitive variables.
c--------------------------------------------------------------------
      subroutine compute_conservatives(rho_s,tau,s_r,rho,P,u,
     &     v,alpha,beta,psi,x,GAMMA,eos_flag,lorentz_max)
      implicit none

      real*8 rho_s, tau, s_r   ! The conservative variables
      real*8 rho, P, v,  u     ! The primitive variables
      real*8 alpha, beta, psi  ! Metric quantities.
      real*8 x                 ! where am I?
      real*8 GAMMA
      integer eos_flag
      
c     Variables for intermediate quantities.
      real*8 ut, ur, u_r, rhoh
      real*8 lorentz_max, factor
      real*8 vmp, vmm
      real*8 sg                ! sqrt(-g)
      real*8 psi4

c     Set metric variables once and for all
      psi4 = psi**4

      ut = -alpha**2 + psi4*(beta**2 + 
     &     2.d0*beta*v + v**2)
      if (ut .ge. 0.d0) then
         write(*,*) 'superluminality '
         factor = (1.d0-1.d0/lorentz_max**2)
         factor = sqrt(factor)
         vmp =  alpha/psi**2*factor - beta
         vmm = -alpha/psi**2*factor - beta
         factor = abs(vmp - v)/abs(vmm - v)
         if (factor > 1) then
            v = vmm
         else
            v = vmp
         end if
         ut = lorentz_max/alpha
      else
         ut = 1.d0/sqrt(-ut)
      end if
      ur = ut*v
      u_r = ut*psi4*(v+beta)

      if (eos_flag == 0) then
         rhoh = rho + GAMMA * P / (GAMMA - 1.d0)
      else if (eos_flag == 1) then
         rhoh = rho + P + u
      end if

      sg    = x**2*psi**6*alpha
      rho_s = sg*rho*ut
      tau   = sg*alpha*(rhoh*ut**2-P/alpha**2)
      s_r   = sg*rhoh*ut*u_r

      end subroutine compute_conservatives

c--------------------------------------------------------------------
c     For a given fluid state, find the maximum left- and right-going
c     wave speeds.  (These are only equal in the comoving frame 
c     of the fluid.  We are calculating grid speeds here.)
c--------------------------------------------------------------------
      subroutine find_cp_cm(cp,cm,cs,alpha,beta,psi,ut,ur)
      implicit none

      real*8 cp,cm,cs
      real*8 alpha,beta,psi
      real*8 gup_tt, gup_tr
      real*8 gup_rr
      real*8 ut, ur
      real*8 a, b, c, d, e, temp, desc

      gup_tt = -1.d0/alpha**2
      gup_tr = beta/alpha**2
      gup_rr = psi**(-4.d0) + beta**2/alpha**2

      d = ut
      a = d**2*(1.d0-cs**2) - cs**2*gup_tt
      e = ur
      b = -2.d0*(d*e*(1.d0-cs**2) - cs**2*gup_tr)
      c = -gup_rr*cs**2 + e**2*(1.d0-cs**2)
      desc = b**2 - 4.d0 * a * c

      if (desc .lt. 0.d0) then
         write(*,*) 'descriminant less than zero! complexity!'
         stop
      end if

      cp = (-b + sqrt(desc))/2.d0/a
      cm = (-b - sqrt(desc))/2.d0/a

      end subroutine find_cp_cm

c-------------------------------------------------------------
c      The refinement grid function.
c-------------------------------------------------------------
      subroutine find_tre_hydro(f_tre,rho,Nx)
      implicit none
      integer Nx,i
      real*8 f_tre(Nx), rho(Nx)
      
      do i = 1, Nx
         if (i==1) then
            f_tre(i) = 0.d0
         else if (i==Nx) then
            f_tre(i) = 0.d0
         else
c     2nd derivative in the x-direction.
            f_tre(i) = rho(i+1) - 2.d0*rho(i) + 
     &           rho(i-1)
         end if
      end do
      return
      end 

      subroutine total_mass(m,rho,x,Nx,ghost_width,my_rank)
      implicit none
      integer Nx
      real*8 m
      real*8 rho(Nx)
      real*8 x(Nx+1)
      real*8 pi
      integer i, j, k
      integer is,ie,js,je,ks,ke
      real*8 dx, dy, dz
      integer ghost_width(6)
      integer my_rank
      
      pi = acos(-1.d0)
      dx = x(2)-x(1)
      is = ghost_width(1) + 1
      ie = Nx - ghost_width(2) 

      m = 0.d0
      do i=is,ie
         m = m + rho(i)
      end do
      
      m = m * dx * 4.d0*pi

      return
      end

c-----------------------------------------------------------------------
c differential operator L, computed where cmask=CMASK_ON (saved in LV)
c
c L_{\alpha}(\alpha) = 0
c L_{\beta}(\beta) = 0
c 
c =>   \alpha'' + [2/r + 2\psi'/\psi]\alpha' + 1/\alpha*(2\psi^4/3)*
c      [2\beta\beta'/r - \beta^2/r^2 - (\beta')^2] 
c       - 4\pi\alpha\psi^4(\rho+S) = 0
c
c      \beta'' + [2/r+6\psi'/\psi-\alpha'/\alpha](\beta' - \beta/r)
c      - 2\pi\alpha\psi^4 j_r = 0
c
c where primes denote differentiation by r.
c We assume V=0 on the physical boundary.
c Notice that you need to plug the vertex centered jr and rhops into here.
c-----------------------------------------------------------------------
        subroutine lop(alpha,Lalpha,beta,Lbeta,psi,jr,rhops,cmask,x,Nx,
     &     alpha_out,phys_bdy)
        implicit none
        integer Nx
        integer phys_bdy(6)
        real*8 Lalpha(Nx+1),alpha(Nx+1),Lbeta(Nx+1),beta(Nx+1)
        real*8 psi(Nx+1),jr(Nx+1),rhops(Nx+1)
        real*8 cmask(Nx+1)
        real*8 x(Nx+1)
        real*8 pi
        include 'cmask.inc'
        integer i,j,k
        real*8 dx,dx2,dalpha,dbeta,dpsi,twodx
        real*8 d2alpha, d2beta, psi4, r, dalphaor
        real*8 alpha_temp, beta_temp, alpha_out

        pi = acos(-1.d0)
        dx=(x(2)-x(1))
        dx2=dx**2
        twodx = 2.d0*dx

        do i=1,Nx+1
           if (cmask(i).eq.CMASK_ON) then
c     Calculate some intermediate quantities.
              r = x(i)
              psi4 = psi(i)**4
              
              if (i==1 .and. phys_bdy(1) .ne. 0) then
c     First order, one sided derivatives
                 dalpha = (alpha(i+1)-alpha(i))/dx
              else if (i==Nx+1 .and. phys_bdy(2) .ne. 0) then
c     Second order, one sided derivatives.
c$$$                 dalpha = 3.d0*alpha(i)-4.d0*alpha(i-1)+alpha(i-2)
c$$$                 dalpha = dalpha / twodx
c$$$                 dbeta  = 3.d0*beta(i) -4.d0*beta(i-1) +beta(i-2)
c$$$                 dbeta  = dbeta / twodx
c     First order, one sided derivatives.
                 dalpha = (alpha(i)-alpha(i-1))/dx
                 dbeta  = (beta(i) -beta(i-1) )/dx
              else
                 dalpha = (alpha(i+1)-alpha(i-1))/twodx
                 dbeta  = (beta(i+1) -beta(i-1)) /twodx
                 dpsi   = (psi(i+1)  -psi(i-1))  /twodx
                 d2alpha = (alpha(i+1) - 2.d0*alpha(i) + alpha(i-1))/dx2
                 d2beta  = (beta(i+1)  - 2.d0*beta(i)  + beta(i-1)) /dx2
              end if

              if (i==1 .and. (phys_bdy(1) .ne. 0)) then
c     For the second order method.
c                Lalpha(i)= -3.d0*alpha(i)+4.d0*alpha(i+1)-alpha(i+2)
c     For the standard method.
                 Lalpha(i) = dalpha
                 Lbeta(i) = beta(i)
              else if (i==Nx+1 .and. (phys_bdy(2) .ne. 0)) then
c     This is for Neumann
c                 Lalpha(i)= dalpha + (alpha(i)-1.d0)/r
c     This is for Dirichlet
                 Lalpha(i)= alpha(i)-alpha_out
                 Lbeta(i) = dbeta  + beta(i)/r
              else                 
                 Lalpha(i)= d2alpha + (2.d0/r + 2.d0*dpsi/psi(i))*dalpha
     &                + 1.d0/alpha(i)*(2.d0*psi4/3.d0)*
     &                (2.d0*beta(i)*dbeta/r-beta(i)**2/r**2-dbeta**2)
     &                - 4.d0*pi*alpha(i)*psi4*rhops(i)

                 Lbeta(i) = d2beta + (2.d0/r + 6.d0*dpsi/psi(i) - 
     &                dalpha/alpha(i))*(dbeta - beta(i)/r) - 
     &                12.d0*pi*alpha(i)*jr(i)
              end if
           end if
        end do

        return
        end

c-----------------------------------------------------------------------
c differential operator L, computed where cmask=CMASK_ON (saved in LV)
c
c L_{\psi}(\psi) = 0
c 
c =>   \psi'' + 2/r\psi' + \frac{\psi^5}{12} [(\beta' -\beta/r)/\alpha]^2
c        + 2\pi\psi^5 \rho = 0
c
c where primes denote differentiation by r.
c Notice that you need to plug the vertex centered ham_source into here.
c-----------------------------------------------------------------------
      subroutine lop_t0(psi0,Lpsi0,alpha,beta,ham_source,cmask,x,Nx,
     &     psi0_out,phys_bdy)
      implicit none
      integer Nx
      integer phys_bdy(6)
      real*8 Lpsi0(Nx+1),psi0(Nx+1),alpha(Nx+1),beta(Nx+1)
      real*8 ham_source(Nx+1)
      real*8 cmask(Nx+1)
      real*8 x(Nx+1)
      real*8 pi, psi0_out
      include 'cmask.inc'
      integer i,j,k
      real*8 dx,dx2,dpsi0,d2psi0,dbeta,twodx
      real*8 psi5, r
      
      pi = acos(-1.d0)
      dx=(x(2)-x(1))
      dx2=dx**2
      twodx = 2.d0*dx

      do i=1,Nx+1
         if (cmask(i).eq.CMASK_ON) then
            r = x(i)
            psi5 = psi0(i)**5
            
            if (i==1 .and. phys_bdy(1) .ne. 0) then
c     First order, one sided derivatives
               dpsi0 = (psi0(i+1)-psi0(i))/dx
c               dpsi0 = (-3.d0*alpha(i)+4.d0*alpha(i+1)-alpha(i+2))/dx
               d2psi0 = 2.d0*(psi0(i+1)-psi0(i))/dx2
            else if (i==Nx+1 .and. phys_bdy(2) .ne. 0) then
c     Second order, one sided derivatives.
c$$$             dpsi0 = 3.d0*psi0(i)-4.d0*psi0(i-1)+psi0(i-2)
c$$$             dpsi0 = dpsi0 / twodx
c     First order, one sided derivatives.
               dpsi0 = (psi0(i)-psi0(i-1))/dx
            else
               dbeta  = (beta(i+1)-beta(i-1))/twodx
               dpsi0  = (psi0(i+1)-psi0(i-1))/twodx
               d2psi0 = (psi0(i+1) - 2.d0*psi0(i) + psi0(i-1))/dx2
            end if

            if (i==1 .and. (phys_bdy(1) .ne. 0)) then
c              Lpsi0(i) = dpsi0
               Lpsi0(i) = 3.d0*d2psi0 + 2.d0*pi*psi5*ham_source(i)
            else if (i==Nx+1 .and. (phys_bdy(2) .ne. 0)) then
c     This is for Neumann
c                 Lpsi0(i)= dpsi0 + (psi0(i)-1.d0)/r
c     This is for Dirichlet
               Lpsi0(i)= psi0(i)-psi0_out
            else                 
               Lpsi0(i)= d2psi0 + 2.d0/r*dpsi0
     &              + psi5/12.d0*((dbeta - beta(i)/r)/alpha(i))**2
     &              + 2.d0*pi*psi5*ham_source(i)
            end if
         end if
      end do
      
      return
      end

c-----------------------------------------------------------------------
c differential operator L, computed where cmask=CMASK_ON (saved in LV)
c This version is to be used for constrained evolution.
c
c L_{\psi}(\psi) = 0
c L_{\alpha}(\alpha) = 0
c L_{\beta}(\beta) = 0
c 
c =>   \psi'' + 2/r\psi' + \frac{\psi^5}{12} [(\beta' -\beta/r)/\alpha]^2
c        + 2\pi\psi^5 \rho = 0
cc
c      \alpha'' + [2/r + 2\psi'/\psi]\alpha' + 1/\alpha*(2\psi^4/3)*
c      [2\beta\beta'/r - \beta^2/r^2 - (\beta')^2] 
c       - 4\pi\alpha\psi^4(\rho+S) = 0
c
c      \beta'' + [2/r+6\psi'/\psi-\alpha'/\alpha](\beta' - \beta/r)
c      - 2\pi\alpha\psi^4 j_r = 0
c
c where primes denote differentiation by r.
c-----------------------------------------------------------------------
      subroutine lop_const_ev(psi,Lpsi,alpha,Lalpha,beta,Lbeta,jr,rhops,
     &     ham_source,cmask,x,Nx,alpha_out,psi_out,phys_bdy)
      implicit none
      integer Nx
      integer phys_bdy(6)
      real*8 Lpsi(Nx+1),psi(Nx+1)
      real*8 Lalpha(Nx+1),alpha(Nx+1),Lbeta(Nx+1),beta(Nx+1)
      real*8 jr(Nx+1),rhops(Nx+1),ham_source(Nx+1)
      real*8 cmask(Nx+1)
      real*8 x(Nx+1)
      real*8 pi
      include 'cmask.inc'
      integer i,j,k
      real*8 dx,dx2,dalpha,dbeta,twodx
      real*8 d2alpha, d2beta, psi4, r, dalphaor, psi5
      real*8 alpha_temp, beta_temp, alpha_out, psi_out
      real*8 dpsi, d2psi
      
      pi = acos(-1.d0)
      dx=(x(2)-x(1))
      dx2=dx**2
      twodx = 2.d0*dx
      
      do i=1,Nx+1
         if (cmask(i).eq.CMASK_ON) then
c     Calculate some intermediate quantities.
            r = x(i)
            psi4 = psi(i)**4
            psi5 = psi(i)*psi4
            
            if (i==1 .and. phys_bdy(1) .ne. 0) then
c     First order, one sided derivatives
               dalpha = (alpha(i+1)-alpha(i))/dx
c                dpsi   = (psi(i+1)-psi(i))/dx
c                dpsi = (-3.d0*psi(i)+4.d0*psi(i+1)-psi(i+2))/dx
               d2psi  = 2.d0*(psi(i+1)-psi(i))/dx2
            else if (i==Nx+1 .and. phys_bdy(2) .ne. 0) then
c     Second order, one sided derivatives.
c$$$             dalpha = 3.d0*alpha(i)-4.d0*alpha(i-1)+alpha(i-2)
c$$$             dalpha = dalpha / twodx
c$$$             dbeta  = 3.d0*beta(i) -4.d0*beta(i-1) +beta(i-2)
c$$$             dbeta  = dbeta / twodx
c$$$             dpsi   = 3.d0*psi(i)-4.d0*psi(i-1)+psi(i-2)
c$$$             dpsi   = dpsi / twodx
c     First order, one sided derivatives.
               dpsi   = (psi(i)-psi(i-1))/dx
               dalpha = (alpha(i)-alpha(i-1))/dx
               dbeta  = (beta(i) -beta(i-1) )/dx
            else
               dpsi    = (psi(i+1)-psi(i-1))/twodx
               dalpha  = (alpha(i+1)-alpha(i-1))/twodx
               dbeta   = (beta(i+1) -beta(i-1)) /twodx
               d2psi   = (psi(i+1) - 2.d0*psi(i) + psi(i-1))/dx2
               d2alpha = (alpha(i+1) - 2.d0*alpha(i) + alpha(i-1))/dx2
               d2beta  = (beta(i+1)  - 2.d0*beta(i)  + beta(i-1)) /dx2
            end if
            
            if (i==1 .and. (phys_bdy(1) .ne. 0)) then
c     For the second order method.
c                Lalpha(i)= -3.d0*alpha(i)+4.d0*alpha(i+1)-alpha(i+2)
c     For the standard method.
               Lpsi(i)   = 3.d0*d2psi + 2.d0*pi*psi5*ham_source(i)
               Lalpha(i) = dalpha
               Lbeta(i)  = beta(i)
            else if (i==Nx+1 .and. (phys_bdy(2) .ne. 0)) then
c     This is for Neumann
c                Lalpha(i) = dalpha + (alpha(i)-1.d0)/r
c                Lpsi(i)   = dpsi +   (psi(i)-1.d0)/r
c     This is for Dirichlet
               Lpsi(i)   = psi(i)-psi_out
               Lalpha(i) = alpha(i)-alpha_out
               Lbeta(i)  = dbeta  + beta(i)/r
            else                 
               Lpsi(i)= d2psi + 2.d0/r*dpsi
     &              + psi5/12.d0*((dbeta - beta(i)/r)/alpha(i))**2
     &              + 2.d0*pi*psi5*ham_source(i)
               
               Lalpha(i)= d2alpha + (2.d0/r + 2.d0*dpsi/psi(i))*dalpha
     &              + 1.d0/alpha(i)*(2.d0*psi4/3.d0)*
     &              (2.d0*beta(i)*dbeta/r-beta(i)**2/r**2-dbeta**2)
     &              - 4.d0*pi*alpha(i)*psi4*rhops(i)
               
               Lbeta(i) = d2beta + (2.d0/r + 6.d0*dpsi/psi(i) - 
     &              dalpha/alpha(i))*(dbeta - beta(i)/r) - 
     &              12.d0*pi*alpha(i)*jr(i)
            end if
         end if
      end do
      
      return
      end
      
c-----------------------------------------------------------------------
c computes the residual L[V]-rhs L, computed where cmask=CMASK_ON 
c (saved in res) ... the L2 norm of the residual is stored in norm
c-----------------------------------------------------------------------
        subroutine residual(alpha,alpha_res,alpha_rhs,beta,beta_res,
     &     beta_rhs,psi,jr,rhops,cmask,x,norm,Nx,alpha_out,phys_bdy)
        implicit none
        integer Nx
        integer phys_bdy(6)
        real*8 alpha(Nx+1),alpha_res(Nx+1),alpha_rhs(Nx+1)
        real*8 beta(Nx+1),beta_res(Nx+1),beta_rhs(Nx+1)
        real*8 psi(Nx+1),jr(Nx+1),rhops(Nx+1),cmask(Nx+1)
        real*8 x(Nx+1),norm,alpha_out

        integer i,j,k,sum
        include 'cmask.inc'

        call lop(alpha,alpha_res,beta,beta_res,psi,jr,rhops,cmask,x,Nx,
     &       alpha_out,phys_bdy)

        norm=0
        sum=0

        do i=1,Nx+1
           if (cmask(i).eq.CMASK_ON) then
              alpha_res(i)=alpha_res(i)-alpha_rhs(i)
              beta_res(i) =beta_res(i) -beta_rhs(i)
              norm=norm+alpha_res(i)**2+beta_res(i)**2
              sum=sum+1
           end if
        end do

        norm=sqrt(norm/sum)
        return
        end

c-----------------------------------------------------------------------
c computes the residual L[V]-rhs L, computed where cmask=CMASK_ON 
c (saved in res) ... the L2 norm of the residual is stored in norm
c-----------------------------------------------------------------------
        subroutine residual_t0(psi0,psi0_res,psi0_rhs,alpha,beta,
     &     ham_source,cmask,x,norm,Nx,psi0_out,phys_bdy)
        implicit none
        integer Nx
        integer phys_bdy(6)
        real*8 psi0(Nx+1), psi0_res(Nx+1), psi0_rhs(Nx+1)
        real*8 alpha(Nx+1)
        real*8 beta(Nx+1)
        real*8 ham_source(Nx+1),cmask(Nx+1)
        real*8 x(Nx+1),norm,psi0_out

        integer i,j,k,sum
        include 'cmask.inc'

        call lop_t0(psi0,psi0_res,alpha,beta,ham_source,cmask,x,Nx,
     &       psi0_out,phys_bdy)

        norm=0
        sum=0

        do i=1,Nx+1
           if (cmask(i).eq.CMASK_ON) then
              psi0_res(i)=psi0_res(i)-psi0_rhs(i)
              norm=norm+psi0_res(i)**2
              sum=sum+1
           end if
        end do

        norm=sqrt(norm/sum)
        return
        end

c-----------------------------------------------------------------------
c computes the residual L[V]-rhs L, computed where cmask=CMASK_ON 
c (saved in res) ... the L2 norm of the residual is stored in norm
c-----------------------------------------------------------------------
        subroutine residual_const_ev(psi,psi_res,psi_rhs,alpha,
     &     alpha_res,alpha_rhs,beta,beta_res,beta_rhs,jr,rhops,
     &     ham_source,cmask,x,norm,Nx,alpha_out,psi_out,phys_bdy)
        implicit none
        integer Nx
        integer phys_bdy(6)
        real*8 psi(Nx+1),psi_res(Nx+1),psi_rhs(Nx+1)
        real*8 alpha(Nx+1),alpha_res(Nx+1),alpha_rhs(Nx+1)
        real*8 beta(Nx+1),beta_res(Nx+1),beta_rhs(Nx+1)
        real*8 jr(Nx+1),rhops(Nx+1),cmask(Nx+1),ham_source(Nx+1)
        real*8 x(Nx+1),norm,alpha_out,psi_out

        integer i,j,k,sum
        include 'cmask.inc'

        call lop_const_ev(psi,psi_res,alpha,alpha_res,beta,beta_res,jr,
     &       rhops,ham_source,cmask,x,Nx,alpha_out,psi_out,phys_bdy)

        norm=0
        sum=0

        do i=1,Nx+1
           if (cmask(i).eq.CMASK_ON) then
              psi_res(i)   = psi_res(i)  - psi_rhs(i)
              alpha_res(i) = alpha_res(i)- alpha_rhs(i)
              beta_res(i)  = beta_res(i) - beta_rhs(i)
              norm=norm+alpha_res(i)**2+beta_res(i)**2+psi_res(i)**2
              sum=sum+1
           end if
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
        subroutine relax(alpha,alpha_rhs,beta,beta_rhs,psi,jr,rhops,
     &     cmask,phys_bdy,x,norm,Nx,ghost_width,alpha_out)
        implicit none
        integer Nx,Ny,Nz
        real*8 alpha(Nx+1),alpha_rhs(Nx+1)
        real*8 beta(Nx+1),beta_rhs(Nx+1)
        real*8 psi(Nx+1),jr(Nx+1),rhops(Nx+1)
        real*8 cmask(Nx+1)
        real*8 x(Nx+1),norm
        real*8 alpha_out  
        integer phys_bdy(6)
        integer ghost_width(6)

        integer i,j,k,N,pass,red,sum
        real*8 pi
        real*8 dx,dx2,twodx
        real*8 res,Jac,r, psi4
        real*8 resalpha, resbeta, jacalpha, jacbeta
        real*8 dalpha, dbeta, dpsi
        real*8 d2alpha, d2beta
        real*8 alpha_temp, beta_temp
        integer have_point
        real*8 norm_alpha, norm_beta ! for debugging
c     for debugging
        real*8 term1, term2, term4, term3

        include 'cmask.inc'

        pi    = acos(-1.d0)
        dx    = (x(2)-x(1))
        dx2   = dx**2
        twodx = 2.d0*dx
        norm  = 0
        sum   = 0
        norm_alpha = 0.d0
        norm_beta  = 0.d0

        do pass = 0,1
           do i=1,Nx+1
              if (mod(i+pass,2).eq.0.and.cmask(i).eq.CMASK_ON) then
                 r = x(i)
                 psi4 = psi(i)**4
                 if (i==1 .and. (phys_bdy(1) .ne. 0)) then
                    dbeta  = beta(i+1)/dx
                    dpsi   = 0.d0
                    d2alpha = 2.d0*(alpha(i+1) - alpha(i))/dx2
                    d2beta  = 2.d0*(beta(i+1)  - beta(i) )/dx2
                    dalpha = (alpha(i+1)-alpha(i))/dx
c     This is the second order method
c$$$                    resalpha = -3.d0*alpha(i) + 4.d0*alpha(i+1) 
c$$$     &                   - alpha(i+2) - alpha_rhs(i)
c$$$                    jacalpha = -3.d0
c     This is the first order method.
                    resalpha = dalpha - alpha_rhs(i)
                    jacalpha = -1.d0/dx
                    
                    resbeta  = beta(i) - beta_rhs(i)
                    jacbeta  = 1.d0
                 else if (i==Nx+1 .and. (phys_bdy(2) .ne. 0)) then
c     Second order in space.
c$$$                    dalpha = 3.d0*alpha(i)-4.d0*alpha(i-1)+alpha(i-2)
c$$$                    dalpha = dalpha / twodx
c$$$                    dbeta  = 3.d0*beta(i) -4.d0*beta(i-1) +beta(i-2)
c$$$                    dbeta  = dbeta / twodx
c     First order in space
                    dalpha = (alpha(i)-alpha(i-1))/dx
                    dbeta  = (beta(i) -beta(i-1) )/dx
c     For Neumann BCs
c                    resalpha = dalpha + (alpha(i)-1.d0)/r - alpha_rhs(i) 
c                    jacalpha = 3.d0/twodx + 1.d0/r
c                    jacalpha = 1.d0/dx + 1.d0/r
c     For Dirichlet BCs
                    resalpha = alpha(i)-alpha_out - alpha_rhs(i)
                    jacalpha = 1.d0

                    resbeta  = dbeta + beta(i)/r - beta_rhs(i)
c                    jacbeta  = 3.d0/twodx + 1.d0/r
                    jacbeta  = 1.d0/dx + 1.d0/r
                 else
                    dalpha = (alpha(i+1)-alpha(i-1))/twodx
                    dbeta  = (beta(i+1) -beta(i-1)) /twodx
                    dpsi   = (psi(i+1)  -psi(i-1))  /twodx
                    d2alpha = (alpha(i+1)-2.d0*alpha(i)+alpha(i-1))/dx2
                    d2beta  = (beta(i+1) -2.d0*beta(i) +beta(i-1)) /dx2

                    resalpha = d2alpha + 
     &                   (2.d0/r + 2.d0*dpsi/psi(i))*dalpha
     &                   + 1.d0/alpha(i)*(2.d0*psi4/3.d0)*(2.d0*beta(i)*
     &                   dbeta/r - beta(i)**2/r**2 - dbeta**2) - 
     &                   4.d0*alpha(i)*pi*psi4*rhops(i) - alpha_rhs(i)

                    resbeta  = d2beta + (2.d0/r + 6.d0*dpsi/psi(i) - 
     &                   dalpha/alpha(i))*(dbeta - beta(i)/r) - 
     &                   12.d0*pi*alpha(i)*jr(i) - beta_rhs(i)
                    
                    jacalpha = -2.d0/dx2 - 1.d0/alpha(i)**2 * 
     &                   (2.d0*psi4/3.d0)*(2.d0*beta(i)*
     &                   dbeta/r - beta(i)**2/r**2 - dbeta**2) - 
     &                   4.d0*pi*psi4*rhops(i)
                    
                    jacbeta  = -2.d0/dx2 - 1.d0/r*(2.d0/r + 
     &                   6.d0*dpsi/psi(i) - dalpha/alpha(i))
                 end if

                 alpha(i) = alpha(i) - resalpha/jacalpha
                 beta(i)  = beta(i)  - resbeta/jacbeta

                 norm_alpha = norm_alpha + resalpha**2
                 norm_beta  = norm_beta  + resbeta**2
                 norm = norm + resalpha**2 + resbeta**2
                 sum = sum + 1
              end if
           end do
        end do

        norm=sqrt(norm/sum)
        norm_alpha = sqrt(norm_alpha/sum)
        norm_beta  = sqrt(norm_beta/sum)

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
      subroutine relax_t0(psi0,psi0_rhs,alpha,beta,ham_source,
     &     cmask,phys_bdy,x,norm,Nx,ghost_width,psi0_out)
      implicit none
      integer Nx,Ny,Nz
      real*8 psi0(Nx+1),psi0_rhs(Nx+1)
      real*8 alpha(Nx+1)
      real*8 beta(Nx+1)
      real*8 ham_source(Nx+1)
      real*8 cmask(Nx+1)
      real*8 x(Nx+1),norm
      real*8 psi0_out  
      integer phys_bdy(6)
      integer ghost_width(6)
      
      integer i,j,k,N,pass,red,sum
      real*8 pi
      real*8 dx,dx2,twodx
      real*8 res,Jac,r, psi4, psi5
      real*8 respsi0, jacpsi0
      real*8 dbeta, dpsi0
      real*8 d2psi0
      integer have_point
c     for debugging
      real*8 term1, term2, term4, term3
      
      include 'cmask.inc'
      
      pi    = acos(-1.d0)
      dx    = (x(2)-x(1))
      dx2   = dx**2
      twodx = 2.d0*dx
      norm  = 0
      sum   = 0
      
      do pass = 0,1
         do i=1,Nx+1
            if (mod(i+pass,2).eq.0.and.cmask(i).eq.CMASK_ON) then
               r = x(i)
               psi4 = psi0(i)**4
               psi5 = psi0(i)*psi4
               
               if (i==1 .and. (phys_bdy(1) .ne. 0)) then
                  dpsi0   = (psi0(i+1)-psi0(i))/dx
                  d2psi0  = 2.d0*(psi0(i+1)-psi0(i))/dx2
c     This is the second order method
c$$$              respsi0 = (-3.d0*psi0(i) + 4.d0*psi0(i+1) 
c$$$  &                   - psi0(i+2))/dx - psi0_rhs(i)
c$$$              jacpsi0 = -3.d0/dx
c     This is the first order method.
c                  respsi0 = dpsi0 - psi0_rhs(i)
c                  jacpsi0 = -1.d0/dx
                  respsi0 = 3.d0*d2psi0 + 2.d0*pi*psi5*ham_source(i) 
     &                 - psi0_rhs(i)
                  jacpsi0  = -6.d0/dx2 + 10.d0*pi*psi4*ham_source(i)
               else if (i==Nx+1 .and. (phys_bdy(2) .ne. 0)) then
c     Second order in space.
c$$$                    dpsi0 = 3.d0*psi0(i)-4.d0*psi0(i-1)+psi0(i-2)
c$$$                    dpsi0 = dpsi0 / twodx
c     First order in space
                  dpsi0 = (psi0(i)-psi0(i-1))/dx
c     For Neumann BCs
c                    respsi0 = dpsi0 + (psi0(i)-1.d0)/r - psi0_rhs(i) 
c                    jacpsi0 = 3.d0/twodx + 1.d0/r
c     For Dirichlet BCs
                  respsi0 = psi0(i)-psi0_out - psi0_rhs(i)
                  jacpsi0 = 1.d0
               else
                  dpsi0  = (psi0(i+1)-psi0(i-1))/twodx
                  dbeta  = (beta(i+1)-beta(i-1)) /twodx
                  d2psi0 = (psi0(i+1)-2.d0*psi0(i)+psi0(i-1))/dx2
                  
                  respsi0 = d2psi0 + 2.d0/r*dpsi0
     &                 + psi5/12.d0*((dbeta - beta(i)/r)/alpha(i))**2
     &                 + 2.d0*pi*psi5*ham_source(i) - psi0_rhs(i)
                  
                  jacpsi0 = -2.d0/dx2 
     &                 + 5.d0*psi4/12.d0
     &                 *((dbeta - beta(i)/r)/alpha(i))**2
     &                 + 10.d0*pi*psi4*ham_source(i)
               end if

               psi0(i) = psi0(i) - respsi0/jacpsi0
               norm = norm + respsi0**2 
               sum = sum + 1
            end if
         end do
      end do
      
      norm=sqrt(norm/sum)

      return
      end

c-----------------------------------------------------------------------
c A version of relaxation for constrained evolution.
c
c the 'approximate' L2 norm of the residual prior to the sweep is
c returned
c
c phys_bdy=1 when corresponding boundary is physical
c-----------------------------------------------------------------------
      subroutine relax_const_ev(psi,psi_rhs,alpha,alpha_rhs,
     &     beta,beta_rhs,jr,rhops,ham_source,cmask,phys_bdy,x,
     &     norm,Nx,ghost_width,alpha_out,psi_out)
      implicit none
      integer Nx,Ny,Nz
      real*8 psi(Nx+1),psi_rhs(Nx+1)
      real*8 alpha(Nx+1),alpha_rhs(Nx+1)
      real*8 beta(Nx+1),beta_rhs(Nx+1)
      real*8 jr(Nx+1),rhops(Nx+1),ham_source(Nx+1)
      real*8 cmask(Nx+1)
      real*8 x(Nx+1),norm
      real*8 alpha_out,psi_out  
      integer phys_bdy(6)
      integer ghost_width(6)
      
      integer i,j,k,N,pass,red,sum
      real*8 pi
      real*8 dx,dx2,twodx
      real*8 res,Jac,r, psi4, psi5
      real*8 respsi, jacpsi
      real*8 resalpha, resbeta, jacalpha, jacbeta
      real*8 dalpha, dbeta, dpsi
      real*8 d2alpha, d2beta, d2psi
      real*8 alpha_temp, beta_temp
      integer have_point
      real*8 norm_alpha, norm_beta, norm_psi 
      real*8 term1, term2, term3, term4
      
      include 'cmask.inc'
      
      pi    = acos(-1.d0)
      dx    = (x(2)-x(1))
      dx2   = dx**2
      twodx = 2.d0*dx
      norm  = 0
      sum   = 0
      norm_alpha = 0.d0
      norm_beta  = 0.d0
      norm_psi   = 0.d0
      
      do pass = 0,1
         do i=1,Nx+1
            if (mod(i+pass,2).eq.0.and.cmask(i).eq.CMASK_ON) then
               r = x(i)
               psi4 = psi(i)**4
               psi5 = psi4*psi(i)

               if (i==1 .and. (phys_bdy(1) .ne. 0)) then
                  dalpha   = (alpha(i+1)-alpha(i))/dx
                  dbeta    = beta(i+1)/dx
                  dpsi     = 0.d0
                  d2psi    = 2.d0*(psi(i+1)   - psi(i)  )/dx2
                  d2alpha  = 2.d0*(alpha(i+1) - alpha(i))/dx2
                  d2beta   = 2.d0*(beta(i+1)  - beta(i) )/dx2

                  respsi = 3.d0*d2psi + 2.d0*pi*psi5*ham_source(i) 
     &                 - psi_rhs(i)
                  jacpsi  = -6.d0/dx2 + 10.d0*pi*psi4*ham_source(i)
                 
                  resalpha = dalpha - alpha_rhs(i)
                  jacalpha = -1.d0/dx
                  
                  resbeta  = beta(i) - beta_rhs(i)
                  jacbeta  = 1.d0
               else if (i==Nx+1 .and. (phys_bdy(2) .ne. 0)) then
                  dpsi   = (psi(i)  - psi(i-1))/dx
                  dalpha = (alpha(i)- alpha(i-1))/dx
                  dbeta  = (beta(i) - beta(i-1) )/dx
c     For Neumann BCs
c                  respsi   = dpsi + (psi(i)-1.d0)/r - psi_rhs(i) 
c                  jacpsi   = 1.d0/dx + 1.d0/r
c                  resalpha = dalpha + (alpha(i)-1.d0)/r - alpha_rhs(i) 
c                  jacalpha = 1.d0/dx + 1.d0/r
c     For Dirichlet BCs
                  respsi = psi(i)-psi_out - psi_rhs(i)
                  jacpsi = 1.d0

                  resalpha = alpha(i)-alpha_out - alpha_rhs(i)
                  jacalpha = 1.d0

                  resbeta  = dbeta + beta(i)/r - beta_rhs(i)
                  jacbeta  = 1.d0/dx + 1.d0/r
               else
                  dpsi    = (psi(i+1)  -psi(i-1)  )/twodx
                  dalpha  = (alpha(i+1)-alpha(i-1))/twodx
                  dbeta   = (beta(i+1) -beta(i-1) )/twodx
                  d2psi   = (psi(i+1)  -2.d0*psi(i)  +psi(i-1)  )/dx2
                  d2alpha = (alpha(i+1)-2.d0*alpha(i)+alpha(i-1))/dx2
                  d2beta  = (beta(i+1) -2.d0*beta(i) +beta(i-1) )/dx2
                  
                  respsi = d2psi + 2.d0/r*dpsi
     &                 + psi5/12.d0*((dbeta - beta(i)/r)/alpha(i))**2
     &                 + 2.d0*pi*psi5*ham_source(i) - psi_rhs(i)
                  
                  resalpha = d2alpha + 
     &                 (2.d0/r + 2.d0*dpsi/psi(i))*dalpha
     &                 + 1.d0/alpha(i)*(2.d0*psi4/3.d0)*(2.d0*beta(i)*
     &                 dbeta/r - beta(i)**2/r**2 - dbeta**2) - 
     &                 4.d0*alpha(i)*pi*psi4*rhops(i) - alpha_rhs(i)
                  
                  resbeta  = d2beta + (2.d0/r + 6.d0*dpsi/psi(i) - 
     &                 dalpha/alpha(i))*(dbeta - beta(i)/r) - 
     &                 12.d0*pi*alpha(i)*jr(i) - beta_rhs(i)

                  jacpsi = -2.d0/dx2 
     &                 + 5.d0*psi4/12.d0
     &                 *((dbeta - beta(i)/r)/alpha(i))**2
     &                 + 10.d0*pi*psi4*ham_source(i)
                  
                  jacalpha = -2.d0/dx2 - 1.d0/alpha(i)**2 * 
     &                 (2.d0*psi4/3.d0)*(2.d0*beta(i)*
     &                 dbeta/r - beta(i)**2/r**2 - dbeta**2) - 
     &                 4.d0*pi*psi4*rhops(i)
                  
                  jacbeta  = -2.d0/dx2 - 1.d0/r*(2.d0/r + 
     &                 6.d0*dpsi/psi(i) - dalpha/alpha(i))
               end if
                    
               psi(i)   = psi(i)   - respsi/jacpsi
               alpha(i) = alpha(i) - resalpha/jacalpha
               beta(i)  = beta(i)  - resbeta/jacbeta
               
               norm_psi   = norm_psi   + respsi**2
               norm_alpha = norm_alpha + resalpha**2
               norm_beta  = norm_beta  + resbeta**2
               norm = norm + respsi**2 + resalpha**2 + resbeta**2
               sum = sum + 1
            end if
         end do
      end do
      
      norm       = sqrt(norm/sum)
      norm_psi   = sqrt(norm_psi/sum)
      norm_alpha = sqrt(norm_alpha/sum)
      norm_beta  = sqrt(norm_beta/sum)

      return
      end

c----------------------------------------------------------      
c     Derivs subroutine for the TOV solution.
c----------------------------------------------------------      
      subroutine derivs(x,y,dydx,kill_flag,gamma,Ye,temp,
     &     rholoc,rho0_scale,il,eos_flag,p_atm)
      integer NMAX,kill_flag
      parameter (NMAX = 50)
      real*8 x, y(NMAX), dydx(NMAX)
      real*8 pi, m, P, phi, rho, rho0, gamma, eps
      real*8 r_iso
      real*8 Ye,temp,rholoc,rho0_scale,p_atm
      real*8 ploc,epsloc,csloc
      integer il,jl,kl, eos_flag

      pi = acos(-1.d0)
      m     = y(1)
      P     = y(2)
      phi   = y(3)

      jl = 1
      kl = 1

      if (P < p_atm) then
         dydx(1) = 0.d0
         dydx(2) = 0.d0
         dydx(3) = 0.d0
         kill_flag = 1
      end if

      if (kill_flag .eq. 0) then

         if (eos_flag == 0) then
            rho0 = P**(1.d0/gamma)
            eps  = P/((gamma-1.d0)*rho0)
         else
            ploc = P
            call finddens_wrapper(ploc,rholoc,epsloc,Ye,temp,
     &           rho0_scale,csloc,il,jl,kl)
            rho0 = rholoc
            eps  = epsloc
         end if
         rho  = rho0*(1.d0+eps)
         
         dydx(1) = 4.d0*pi*x**2*rho
         
         if (m==0.d0) then
            dydx(2) = 0.d0
            dydx(3) = 0.d0
         else
            dydx(2) = -rho*m/x**2*(1.d0+P/rho)
            dydx(2) = dydx(2)*(1+(4.d0*pi*P*x**3)/m)
            dydx(2) = dydx(2)/(1-2.d0*m/x)
            
            dydx(3) = -dydx(2)/rho
            dydx(3) = dydx(3)/(1.d0+P/rho)
         end if
      end if

      end subroutine

      subroutine derivs2(x,y,dydx,m_s_ex,r_s_ex,NRS_ex)
      integer NMAX,kill_flag,NRS_ex
      parameter (NMAX = 50)
      real*8 x, y(NMAX), dydx(NMAX), m, r_iso, xloc
      real*8 m_s_ex(NRS_ex), r_s_ex(NRS_ex)
      real*8 ya(4), xa(4), mloc, dy
      integer iloc
      xloc = x
c     First, interpolate to find the local mass.
      do i = 1, NRS_ex
         if (xloc .lt. r_s_ex(i)) exit
      end do
      iloc = i

      if (iloc == 1 .or. iloc == 2) then
         ya(1) = m_s_ex(1)
         ya(2) = m_s_ex(2)
         ya(3) = m_s_ex(3)
         ya(4) = m_s_ex(4)
         xa(1) = r_s_ex(1)
         xa(2) = r_s_ex(2)
         xa(3) = r_s_ex(3)
         xa(4) = r_s_ex(4)
      else if (iloc == NRS_ex) then
         ya(1) = m_s_ex(NRS_ex-3)
         ya(2) = m_s_ex(NRS_ex-2)
         ya(3) = m_s_ex(NRS_ex-1)
         ya(4) = m_s_ex(NRS_ex)
         xa(1) = r_s_ex(NRS_ex-3)
         xa(2) = r_s_ex(NRS_ex-2)
         xa(3) = r_s_ex(NRS_ex-1)
         xa(4) = r_s_ex(NRS_ex)
      else
         ya(1) = m_s_ex(iloc-2)
         ya(2) = m_s_ex(iloc-1)
         ya(3) = m_s_ex(iloc)
         ya(4) = m_s_ex(iloc+1)
         xa(1) = r_s_ex(iloc-2)
         xa(2) = r_s_ex(iloc-1)
         xa(3) = r_s_ex(iloc)
         xa(4) = r_s_ex(iloc+1)
      end if         

      call polint(xa,ya,4,xloc,mloc,dy)

      r_iso = y(1)
      if (x .eq. 0.d0) then
         dydx(1) = 1.d0
      else
         dydx(1) = r_iso/(x*sqrt(1.d0-2.d0*mloc/x))
      end if

      end subroutine

      subroutine rk4(y,dydx,n,x,h,yout,derivs,kill_flag,gamma,
     &     Ye,temp,rholoc,rho0_scale,iloc,eos_flag,p_atm)
      integer n, NMAX
      real*8 h, x, dydx(n), y(n), yout(n)
      external derivs
      parameter (NMAX = 50)
      integer i, kill_flag
      real*8 gamma,Ye,temp,rholoc,rho0_scale,p_atm
      integer iloc, eos_flag
      real*8 h6, hh, xh, dym(NMAX), dyt(NMAX), yt(NMAX)
      
      hh = h*0.5d0
      h6 = h / 6.d0
      xh = x + hh

      do i = 1,n
         yt(i) = y(i)+hh*dydx(i)
      end do

      call derivs(xh,yt,dyt,kill_flag,gamma,Ye,temp,
     &     rholoc,rho0_scale,iloc,eos_flag,p_atm)
      if (kill_flag .ne. 0) return

      do i = 1,n
         yt(i) = y(i)+hh*dyt(i)
      end do

      call derivs(xh,yt,dym,kill_flag,gamma,Ye,temp,
     &     rholoc,rho0_scale,iloc,eos_flag,p_atm)
      if (kill_flag .ne. 0) return

      do i = 1,n
         yt(i) = y(i)+h*dym(i)
         dym(i) = dyt(i)+dym(i)
      end do

      call derivs(x+h,yt,dyt,kill_flag,gamma,Ye,temp,
     &     rholoc,rho0_scale,iloc,eos_flag,p_atm)
      if (kill_flag .ne. 0) return
      
      do i = 1,n
         yout(i) = y(i) + h6*(dydx(i)+dyt(i)+2.d0*dym(i))
      end do
      
      return 
      end

      subroutine rk4_1real(y,dydx,n,x,h,yout,m_s_ex,r_s_ex,NRS_ex,
     &     derivs)
      integer n, NMAX, NRS_ex
      real*8 h, x, dydx(n), y(n), yout(n), m
      real*8 m_s_ex(NRS_ex), r_s_ex(NRS_ex)
      external derivs
      parameter (NMAX = 50)
      integer i
      real*8 h6, hh, xh, dym(NMAX), dyt(NMAX), yt(NMAX)
      
      hh = h*0.5d0
      h6 = h / 6.d0
      xh = x + hh

      do i = 1,n
         yt(i) = y(i)+hh*dydx(i)
      end do

      call derivs(xh,yt,dyt,m_s_ex,r_s_ex,NRS_ex)

      do i = 1,n
         yt(i) = y(i)+hh*dyt(i)
      end do

      call derivs(xh,yt,dym,m_s_ex,r_s_ex,NRS_ex)

      do i = 1,n
         yt(i) = y(i)+h*dym(i)
         dym(i) = dyt(i)+dym(i)
      end do

      call derivs(x+h,yt,dyt,m_s_ex,r_s_ex,NRS_ex)
      
      do i = 1,n
         yout(i) = y(i) + h6*(dydx(i)+dyt(i)+2.d0*dym(i))
      end do
      
      return 
      end

      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END

c-----------------------------------------------------------------------
c     Evaluates the Hamiltonian constraint and returns a norm.
c-----------------------------------------------------------------------
      subroutine ham_const(ham,alpha,beta,psi,tau_v,ham_source_v,
     &     cmask,phys_bdy,x,norm,Nx,ghost_width,out_name,write_profile)
      implicit none
      integer Nx
      real*8 ham(Nx+1)
      real*8 ham_source_v(Nx+1)
      real*8 alpha(Nx+1)
      real*8 beta(Nx+1)
      real*8 psi(Nx+1),tau_v(Nx+1)
      real*8 cmask(Nx+1)
      real*8 x(Nx+1),norm
      integer phys_bdy(6)
      integer ghost_width(6)
      character(len=30) out_name
      integer write_profile
      
      integer i,j,k,N,sum
      real*8 pi
      real*8 dx,dx2,twodx
      real*8 res,r, psi4
      real*8 dalpha, dbeta, dpsi
      real*8 d2psi
      real*8 term1, term2, term3, term4
      real*8 loc_norm
      integer write_step
      include 'cmask.inc'
      
      pi    = acos(-1.d0)
      dx    = (x(2)-x(1))
      dx2   = dx**2
      twodx = 2.d0*dx
      norm  = 0
      sum   = 0

      write_step = Nx/256

      if (write_profile .ne. 0) then
         open (unit=1,file=out_name,status='unknown')
      end if

      do i=1,Nx
         if (cmask(i).eq.CMASK_ON) then
            r = x(i)
            psi4 = psi(i)**4

            if (r == 0.d0) then
               d2psi  = (2.d0*psi(i+1)-2.d0*psi(i))/dx2
               term2  = 2.d0*d2psi
               term3  = 0.d0
            else
               dalpha = (alpha(i+1)-alpha(i-1))/twodx
               dbeta  = (beta(i+1) -beta(i-1)) /twodx
               dpsi   = (psi(i+1)  -psi(i-1))  /twodx
               d2psi  = (psi(i+1) - 2.d0*psi(i) + psi(i-1))/dx2

               term2 = 2.d0*dpsi/r
               term3 = (dbeta-beta(i)/r)/alpha(i)
               term3 = term3*term3
               term3 = term3*psi4*psi(i)/12.d0
            end if
            
            term1 = d2psi
            term4 = 2.d0*pi*psi4*psi(i)*ham_source_v(i)
            ham(i) = term1 + term2 + term3 + term4
            if (write_profile .ne. 0) then
               write(1,101) x(i),ham(i),term1,term2,term3,term4
            end if
            loc_norm = term1**2 + term2**2 + term3**2 + term4**2
            loc_norm = sqrt(loc_norm)
            
            ham(i) = ham(i)/loc_norm
            norm = norm + ham(i)**2
            sum = sum + 1
            
         end if
      end do
      norm=sqrt(norm/sum)
 101  format (6e17.5)
      return
      end


#ifndef _NUM
#define _NUM 1
/*======================================================================================== */
/* num.h                                                                                   */
/*                                                                                         */
/* c-prototypes for the fortran routines                                                   */
/*======================================================================================== */
void hydro_1step_(real *dn, real *srn, real *taun, 
		  real *d_p, real *sr_p, real *tau_p, real *rho, real *v, real *u,
		  real *P, real *T, real *cs, real *jr, real *rhops, real *ham_source,
		  real *x, real *x_c, real *gamma, int *limiter, real *Ye, real *dt, real *rho0, 
		  int *Nx, real *psi, real *psi_c, real *alpha, real *alpha_c,
		  real *beta, real *beta_c, int *fc_mask, real *d_fcs, real *sr_fcs, 
		  real *tau_fcs, int *save_fcs, int *eos_flag, int *phys_bdy,
		  real *deltat1, real *deltat2, real *rho_atm, real *p_atm,
		  real *u_atm, real *t_atm, real *lorentz_max);

void psi_1step_(real *psin, real *psi_p, real *beta, real *x, int *Nx, int *phys_bdy, int *iter);

//void psi_save_(real *psinp1, real *psinp1_c, real *psi_v_sav, real *psi_sav, int *Nx);

void id_tov_(real *gamma, real *rho_c, real *atm_frac, real *dn,
	     real *srn, real *taun, real *rho, real *v, real *u, real *P, real *T, 
	     real *jr, real *jr_v, real *rhops, real *rhops_v, real *tau_v, 
	     real *ham_source, real *ham_source_v,
	     real *alpha, real *beta, real *psi, real *alpha_c,
	     real *beta_c, real *psi_c, real *x, real *x_c, int *Nx,
	     real *ham, real *mask, int *phys_bdy, int *ghost_width, real *p_deplete,
	     real *Ye, real *temperature, real *rho0, int *eos_flag, real *alpha_out, 
	     real *rho_atm, real *p_atm, real *u_atm, real *t_atm, real *deltat1, 
	     real *deltat2, real *v_pert, real *lorentz_max, real *psi0_out);

void find_primitives_(real *dn, real *srn, real *taun, real *rho, 
		      real *v, real *u, real *jr, real *rhops, real *ham_source, 
		      real *gamma, real *alpha_c, real *beta_c, real *psi_c, real *x, 
		      int *Nx, real *rho_atm, real *lorentz_max);

void find_primitives_table_(real *dn, real *srn, real *taun, real *rho, real *v, real *u, 
			    real *jr, real *rhops, real *ham_source, real *P, real *T, real *cs, 
			    real *Ye, real *rho0, real *alpha_c, real *beta_c, real *psi_c, 
			    real *x, int *Nx, real *deltat1, real* deltat2, real *rho_atm,
			    real *p_atm, real *u_atm, real *t_atm, real *lorentz_max);

void init_rest_(real *dn, real *dnp1, real *srn, real *srnp1, real *taun, 
		real *taunp1, int *Nx);

void zero_fcs_vars_(int *type_flag, int *fc_mask, real *d_fcs, real *sr_fcs,
		    real *tau_fcs, int *Nx);

void apply_fc_(real *dn, real *srn, real *taun, real *d_fcs,
	      real *sr_fcs, real *tau_fcs, int *Nx);

void find_tre_hydro_(real *f_tre, real *rho, int *Nx);

void total_mass_(real *m, real *rho, real *x, int *Nx, int *ghost_width, int *my_rank);

void residual_(real *alpha, real *alpha_res, real *alpha_rhs, real *beta,
	       real *beta_res, real *beta_rhs, real *psi, real *jr, real *rhops,
	       real *cmask, real *x, real *norm, int *Nx, real *alpha_out, int *phys_bdy);

void residual_t0_(real *psi0, real *psi0_res, real *psi0_rhs, real *alpha,
		  real *beta, real *ham_source, real *cmask, real *x, 
		  real *norm, int *Nx, real *psi0_out, int *phys_bdy);

void residual_const_ev_(real *psi, real *psi_res, real *psi_rhs, 
			real *alpha, real *alpha_res, real *alpha_rhs, 
			real *beta, real *beta_res, real *beta_rhs, 
			real *jr, real *rhops, real *ham_source,
			real *cmask, real *x, real *norm, int *Nx, 
			real *alpha_out, real *psi_out, int *phys_bdy);

void lop_(real *alpha, real *Lalpha, real *beta, real *Lbeta, real *psi, real *jr,
	  real *rhops, real *cmask, real *x, int *Nx, real *alpha_out, int *phys_bdy);

void lop_t0_(real *psi0, real *Lpsi0, real *alpha, real *beta, real *ham_source,
	     real *cmask, real *x, int *Nx, real *psi_out, int *phys_bdy);

void lop_const_ev_(real *psi, real *Lpsi, real *alpha, real *Lalpha, 
		   real *beta, real *Lbeta, real *jr, real *rhops, real *ham_source, 
		   real *cmask, real *x, int *Nx, real *alpha_out, real *psi_out,
		   int *phys_bdy);

void relax_(real *alpha, real *alpha_rhs, real *beta, real *beta_rhs, real *psi,
	    real *jr, real *rhops, real *cmask, int *phys_bdy, real *x, real *norm,
	    int *Nx, int *ghost_width, real *alpha_out);

void relax_t0_(real *psi0, real *psi0_rhs, real *alpha, real *beta, real *ham_source,
	       real *cmask, int *phys_bdy, real *x, real *norm,int *Nx, 
	       int *ghost_width, real *psi0_out);

void relax_const_ev_(real *psi, real *psi_rhs, real *alpha, real *alpha_rhs, 
		     real *beta, real *beta_rhs, real *jr, real *rhops, real *ham_source, 
		     real *cmask, int *phys_bdy, real *x, real *norm,
		     int *Nx, int *ghost_width, real *alpha_out, real *psi_out);

void ham_const_(real *ham, real *alpha, real *beta, real *psi, real *tau_v, real *ham_source_v,
		real *cmask, int *phys_bdy, real *x, real *norm, int *Nx,
		int *ghost_width, char *out4_name, int *write_profile);

#endif /*_NUM */

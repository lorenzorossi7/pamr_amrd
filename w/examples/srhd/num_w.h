#ifndef _NUM
#define _NUM 1
/*======================================================================================== */
/* num.h                                                                                   */
/*                                                                                         */
/* c-prototypes for the fortran routines                                                   */
/*======================================================================================== */
void hydro_1step_(real *dnp1, real *sxnp1, real *synp1, real *sznp1, real *enp1, 
		  real *d_p, real *sx_p, real *sy_p, real *sz_p, real *e_p,
		  real *rho, real *vx, real *vy, real *vz, real *u,
		  real *P, real *T, real *cs,
		  real *x, real *y, real *z, real *gamma, int *limiter,
		  real *dt, real *Ye, real *rho0, int *Nx, int *Ny, int *Nz, 
		  int *fc_mask, real *d_fcs, real *sx_fcs, real *sy_fcs,
		  real *sz_fcs, real *e_fcs, int *save_fcs, int *eos_flag, int *phys_bdy,
		  int *bc_type, int *N_bound, real *lorentz_max);

void init_riemann_(real *gamma, real *dn, real *sxn, real *syn, real *szn,
		   real *en, real *rho, real *vx, real *vy, 
		   real *vz, real *u, real *P, real *T, real *cs,
		   real *x_c, real *y_c, real *z_c,
		   int *Nx, int *Ny, int *Nz, real *x0, real *rho_l, real *rho_r, 
		   real *u_l, real *u_r, real *ut_l, real *ux_l, real *uy_l, real *uz_l,
		   real *ut_r, real *ux_r, real *uy_r, real *uz_r, real *T_r, real *T_l,
		   real *Ye, real *rho0, int *eos_flag, int *riemann_prob, real *lorentz_max);

void find_primitives_(real *dn, real *sxn, real *syn, real *szn,   
		      real *en, real *rho, 
		      real *vx, real *vy, real *vz, real *u,
		      real *gamma, int *Nx, int *Ny, int *Nz);

void find_primitives_table_(real *dn, real *sxn, real *syn, real *szn,   
			    real *en, real *rho, 
			    real *vx, real *vy, real *vz, real *u, 
			    real *P, real *T, real *cs, real *Ye, real *rho0,
			    int *Nx, int *Ny, int *Nz);

void init_rest_(real *dn, real *dnp1, real *sxn, real *sxnp1, real *syn, real *synp1,
		real *szn, real *sznp1, real *en, real *enp1, int *Nx, int *Ny, int *Nz);

void zero_fcs_vars_(int *type_flag, int *fc_mask, real *d_fcs, real *sx_fcs,
		    real *sy_fcs, real *sz_fcs, real *e_fcs, int *Nx, int *Ny, int *Nz);

void apply_fc_(real *dn, real *sxn, real *syn, real *szn, real *en, real *d_fcs,
	      real *sx_fcs, real *sy_fcs, real *sz_fcs, real *e_fcs, 
	      int *Nx, int *Ny, int *Nz);

void find_tre_hydro_(real *f_tre, real *rho, int *Nx, int *Ny, int *Nz);

void total_mass_(real *m, real *rho, real *x, real *y, real *z, 
		 int *Nx, int *Ny, int *Nz, int *ghost_width, int *my_rank);

void cons_bc_(int *bc_type, real *dn, real *en, real *sxn, real *syn, real *szn, 
		 int *Nx, int *Ny, int *Nz, int *N_bound, int *phys_bdy);

#endif /*_NUM */

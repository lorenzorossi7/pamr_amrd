#ifndef _NUM
#define _NUM 1
/*======================================================================================== */
/* num.h                                                                                   */
/*                                                                                         */
/* c-prototypes for the fortran routines                                                   */
/*======================================================================================== */
void lop_(real *LV, real *V, real *rho, real *cmask, 
          real *x, real *y, real *z, int *Nx, int *Ny, int *Nz);

void residual_(real *res, real *rhs, real *V, real *rho, real *cmask, 
               real *x, real *y, real *z, real *norm, int *Nx, int *Ny, int *Nz);

void residual_diagnostic_(real *res, real *cmask, 
               real *x, real *y, real *z, int *Nx, int *Ny, int *Nz, int *ghost_width);

void diagnose_res_(real *res, real *cmask, 
		   real *x, real *y, real *z, int *Nx, int *Ny, int *Nz, int *ghost_width);

void relax_(real *V, real *rhs, real *rho, real *cmask, int *phys_bdy,
            real *x, real *y, real *z, real *norm, int *Nx, int *Ny, int *Nz, int *ghost_width);

void hydro_1step_(int *step_flag, real *d, real *sx, real *sy, real *sz, 
		  real *en, real *phi, real *d_p, real *sx_p, real *sy_p, real *sz_p, real *e_p,
		  real *rho, real *vx, real *vy, real *vz, real *u,
		  real *x, real *y, real *z, real *gamma, int *limiter,
		  real *rho_atm, real *dt, int *Nx, int *Ny, int *Nz, real *p_deplete,
		  int *fc_mask, real *d_fcs, real *sx_fcs, real *sy_fcs,
		  real *sz_fcs, real *e_fcs);

void polytrope_id_(real *gamma, real *dn, real *sxn, real *syn, real *szn,
		   real *en, real *rho, real *rho_v, real *vx, real *vy, 
		   real *vz, real *u, real *phiv, real *p_deplete, 
		   real *x, real *y, real *z, 
		   real *x_c, real *y_c, real *z_c,
		   int *Nx, int *Ny, int *Nz, real *rho_atm);

void find_primitives_(real *dn, real *sxn, real *syn, real *szn,   
		      real *en, real *phin, real *rho, 
		      real *vx, real *vy, real *vz, real *u,
		      real *gamma, int *Nx, int *Ny, int *Nz,
		      real *rho_atm, real *p_deplete);

void total_mass_(real *m, real *rho, real *x, real *y, real *z, 
		 int *Nx, int *Ny, int *Nz, int *ghost_width, int *my_rank);

void get_rho_c_(int *has_origin, real *rho_c, real *rho, real *x, real *y, real *z, 
		int *Nx, int *Ny, int *Nz, int *ghost_width);

void init_rest_(real *dn, real *dnp1, real *sxn, real *sxnp1, real *syn, real *synp1,
		real *szn, real *sznp1, real *en, real *enp1, real *phiv, 
		real *phiv_np1, real *phi_extrap_tm1, 
		real *phi_extrap_tm2, int *Nx, int *Ny, int *Nz);

void zero_fcs_vars_(int *type_flag, int *fc_mask, real *d_fcs, real *sx_fcs,
		    real *sy_fcs, real *sz_fcs, real *e_fcs, int *Nx, int *Ny, int *Nz);

void apply_fc_(real *dn, real *sxn, real *syn, real *szn, real *en, real *d_fcs,
	      real *sx_fcs, real *sy_fcs, real *sz_fcs, real *e_fcs, 
	      int *Nx, int *Ny, int *Nz);

#endif /*_NUM */

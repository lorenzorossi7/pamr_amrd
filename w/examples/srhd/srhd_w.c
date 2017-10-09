//=============================================================================
// application functions/variables for Newtonian Boson Stars
//=============================================================================

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pamr_w.h"
#include "amrd_w.h"
#include "num_w.h"
#include "shen.h"
#include <mpi.h>
#define INTSIZE 10
//=============================================================================
// some convenient, "local" global variables
//=============================================================================
char *itoa(int value);
double gamma;
int limiter;
int bc_type;
int N_bound;
double lorentz_max;

double *rho, *u, *vx, *vy, *vz, *T, *P, *cs;
double *d, *Sx, *Sy, *Sz, *E;
double *d_p, *Sx_p, *Sy_p, *Sz_p, *E_p;
double *d_fcs, *Sx_fcs, *Sy_fcs, *Sz_fcs, *E_fcs;
double *d_tre;

double *d_n, *Sx_n, *Sy_n, *Sz_n, *E_n;
double *d_np1, *Sx_np1, *Sy_np1, *Sz_np1, *E_np1;

double *x,*y,*z, *x_c, *y_c, *z_c;
int shape[3],shape_c[3],ghost_width[6],Nx,Ny,Nz,phys_bdy[6],size,size_c;
double base_bbox[6],bbox[6],dx,dy,dz,dt,dto2;
int g_L;
int loc_lev;

int rho_gfn, u_gfn, vx_gfn, vy_gfn, vz_gfn, T_gfn, P_gfn, cs_gfn;
int d_gfn, Sx_gfn, Sy_gfn, Sz_gfn, E_gfn;
int d_p_gfn, Sx_p_gfn, Sy_p_gfn, Sz_p_gfn, E_p_gfn;
int d_fcs_gfn, Sx_fcs_gfn, Sy_fcs_gfn, Sz_fcs_gfn, E_fcs_gfn;
int d_tre_gfn;

int d_n_gfn, Sx_n_gfn, Sy_n_gfn, Sz_n_gfn, E_n_gfn;
int d_np1_gfn, Sx_np1_gfn, Sy_np1_gfn, Sz_np1_gfn, E_np1_gfn;

int max_lev;
int init_depth;
int dim;
double t0;
int ID_DV_trace;
int eos_flag;
double Ye;  // For now, we'll set this to a constant.
double rho0; // The density scale to pass into the eos.

int num_TRE_max_hydro;
real *TRE_max_hydro;

// Parameters for the Riemann problem.
double x0;
double rho_l, rho_r;
double u_l, u_r;
double ut_l,ux_l,uy_l,uz_l;
double ut_r,ux_r,uy_r,uz_r;
double T_l, T_r;
int riemann_prob;

void srhd_free_data(void);
void srhd_init_rest(void);
void srhd_find_primitives(void);
void srhd_fcs_var_clear(int type_flag, int *ifc_mask);
void srhd_flux_correct(void);
void srhd_post_flux_correct(void);
void srhd_pre_tstep(int L);

//=============================================================================
// call after variables have been defined
//=============================================================================
void set_gfns(void)
{
    if ((d_gfn=PAMR_get_gfn("d",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((d_np1_gfn=PAMR_get_gfn("d",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((d_n_gfn=PAMR_get_gfn("d",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((Sx_gfn=PAMR_get_gfn("Sx",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((Sx_np1_gfn=PAMR_get_gfn("Sx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Sx_n_gfn=PAMR_get_gfn("Sx",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((Sy_gfn=PAMR_get_gfn("Sy",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((Sy_np1_gfn=PAMR_get_gfn("Sy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Sy_n_gfn=PAMR_get_gfn("Sy",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((Sz_gfn=PAMR_get_gfn("Sz",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((Sz_np1_gfn=PAMR_get_gfn("Sz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Sz_n_gfn=PAMR_get_gfn("Sz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((E_gfn=PAMR_get_gfn("E",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((E_np1_gfn=PAMR_get_gfn("E",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((E_n_gfn=PAMR_get_gfn("E",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    // Stuff for flux correction.  They are all work variables in the AMR hierarchy.
    if ((d_fcs_gfn   = PAMR_get_gfn("d_fcs"  , PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((Sx_fcs_gfn  = PAMR_get_gfn("Sx_fcs" , PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((Sy_fcs_gfn  = PAMR_get_gfn("Sy_fcs" , PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((Sz_fcs_gfn  = PAMR_get_gfn("Sz_fcs" , PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((E_fcs_gfn   = PAMR_get_gfn("E_fcs"  , PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 

    // Some work variables
    if ((d_p_gfn=PAMR_get_gfn("d_p",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((Sx_p_gfn=PAMR_get_gfn("Sx_p",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((Sy_p_gfn=PAMR_get_gfn("Sy_p",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((Sz_p_gfn=PAMR_get_gfn("Sz_p",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((E_p_gfn=PAMR_get_gfn("E_p",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((rho_gfn=PAMR_get_gfn("rho",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((vx_gfn=PAMR_get_gfn("vx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((vy_gfn=PAMR_get_gfn("vy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((vz_gfn=PAMR_get_gfn("vz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((u_gfn=PAMR_get_gfn("u",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((T_gfn=PAMR_get_gfn("T",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((P_gfn=PAMR_get_gfn("P",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((cs_gfn=PAMR_get_gfn("cs",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    
    if ((d_tre_gfn=PAMR_get_gfn("d_tre",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 

}

//=============================================================================
// call with valid iter to set up globals:
//=============================================================================
void ldptr(void)
{
   int rank,dim,ngfs,lev;
   double t,*x0[3],*x0_c[3],*gfs[100],dx0[3];
   static int first=1;

   if (first) 
   {
      first=0; 
      set_gfns();
      PAMR_get_global_bbox(base_bbox);
   }

   if (!(PAMR_get_g_attribs(&rank,&dim,shape,shape_c,bbox,ghost_width,&t,&ngfs,x0,x0_c,gfs))) 
      AMRD_stop("ldptr: PAMR_get_g_attribs failed\n","");

   d=gfs[d_gfn-1];
   d_fcs=gfs[d_fcs_gfn-1];
   d_n=gfs[d_n_gfn-1];
   d_np1=gfs[d_np1_gfn-1];

   Sx=gfs[Sx_gfn-1];
   Sx_fcs=gfs[Sx_fcs_gfn-1];
   Sx_n=gfs[Sx_n_gfn-1];
   Sx_np1=gfs[Sx_np1_gfn-1];

   Sy=gfs[Sy_gfn-1];
   Sy_fcs=gfs[Sy_fcs_gfn-1];
   Sy_n=gfs[Sy_n_gfn-1];
   Sy_np1=gfs[Sy_np1_gfn-1];

   Sz=gfs[Sz_gfn-1];
   Sz_fcs=gfs[Sz_fcs_gfn-1];
   Sz_n=gfs[Sz_n_gfn-1];
   Sz_np1=gfs[Sz_np1_gfn-1];

   E=gfs[E_gfn-1];
   E_fcs=gfs[E_fcs_gfn-1];
   E_n=gfs[E_n_gfn-1];
   E_np1=gfs[E_np1_gfn-1];

   rho   = gfs[rho_gfn-1];
   vx    = gfs[vx_gfn-1];
   vy    = gfs[vy_gfn-1];
   vz    = gfs[vz_gfn-1];
   u     = gfs[u_gfn-1];
   T     = gfs[T_gfn-1];
   P     = gfs[P_gfn-1];
   cs    = gfs[cs_gfn-1];

   d_tre = gfs[d_tre_gfn-1];

   x=x0[0]; dx=x[1]-x[0];
   y=x0[1]; dy=y[1]-y[0];
   z=x0[2]; dz=z[1]-z[0];
   x_c=x0_c[0]; 
   y_c=x0_c[1]; 
   z_c=x0_c[2]; 
   PAMR_get_g_level(&lev);
   loc_lev = lev;
   PAMR_get_dxdt(lev,dx0,&dt);
   dto2 =dt/2.0;

   if ((bbox[0]-base_bbox[0])<dx/2) phys_bdy[0]=1; else phys_bdy[0]=0;
   if ((base_bbox[1]-bbox[1])<dx/2) phys_bdy[1]=1; else phys_bdy[1]=0;
   if ((bbox[2]-base_bbox[2])<dy/2) phys_bdy[2]=1; else phys_bdy[2]=0;
   if ((base_bbox[3]-bbox[3])<dy/2) phys_bdy[3]=1; else phys_bdy[3]=0;
   if ((bbox[4]-base_bbox[4])<dz/2) phys_bdy[4]=1; else phys_bdy[4]=0;
   if ((base_bbox[5]-bbox[5])<dz/2) phys_bdy[5]=1; else phys_bdy[5]=0;

   Nx=shape_c[0];
   Ny=shape_c[1];
   Nz=shape_c[2];

   d_p  = gfs[d_p_gfn-1];
   Sx_p = gfs[Sx_p_gfn-1];
   Sy_p = gfs[Sy_p_gfn-1];
   Sz_p = gfs[Sz_p_gfn-1];
   E_p  = gfs[E_p_gfn-1];

   size_c=Nx*Ny*Nz;
}

//=============================================================================
// utility routines
//=============================================================================
void zero_c(double *f)
{
   int i;
   for (i=0; i<shape_c[0]*shape_c[1]*shape_c[2]; i++) f[i]=0;
}

//=============================================================================
// Routines required by amrd:
//=============================================================================

//=============================================================================
// Returns 0 
//=============================================================================
int srhd_id(void)
{
  return 0;
}

//=============================================================================
// Sets custom parameters, variables, etc. Split up into two segments,
// one called before the pamr context is initialized and standard
// parameters are read, and the other afterwards
//=============================================================================
void srhd_var_pre_init(char *pfile)
{
  // Read in some global parameters which are needed by srhd_id
  max_lev = 1;
  init_depth = 2;
  dim = 3;
  t0 = 0.0;
  base_bbox[0] = base_bbox[2] = base_bbox[4] = -1.0;
  base_bbox[1] = base_bbox[3] = base_bbox[5] = 1.0;
  ID_DV_trace = 0;
  AMRD_int_param(pfile,"dim",&dim,1);
  AMRD_int_param(pfile,"max_lev",&max_lev,1);
  AMRD_int_param(pfile,"init_depth",&init_depth,1);
  AMRD_real_param(pfile,"base_abbox",base_bbox,2*dim);
  AMRD_real_param(pfile,"t0",&t0,1);
  AMRD_int_param(pfile,"ID_DV_trace",&ID_DV_trace,1);

  return;
}

void srhd_var_post_init(char *pfile)
{

   if (my_rank==0)
   {
      printf("===================================================================\n");
      printf("Reading SRHD parameters:\n\n");
   }

   gamma = 1.333333333333333333;
   limiter = 1;
   x0 = 1.0;
   rho_l = 1.0;
   rho_r = 1.0;
   u_l = 1.0;
   u_r = 1.0;
   ut_l = 1.0;
   ut_r = 1.0;
   ux_l = 0.0;
   ux_r = 0.0;
   uy_l = 0.0;
   uy_r = 0.0;
   uz_l = 0.0;
   uz_r = 0.0;
   Ye = 0.1;      // Default to neutron rich matter.
   T_l = 0.0;
   T_r = 0.0;
   rho0 = 1.e14;
   eos_flag = 0;  // Defaults to polytrope.
   num_TRE_max_hydro = 0;
   riemann_prob = 0;
   bc_type = 0;
   N_bound = 0;
   lorentz_max = 10.0;
  
   AMRD_real_param(pfile,"gamma",&gamma,1); 
   AMRD_int_param(pfile,"limiter",&limiter,1); 
   AMRD_real_param(pfile,"x0",&x0,1);
   AMRD_real_param(pfile,"rho_l",&rho_l,1);
   AMRD_real_param(pfile,"rho_r",&rho_r,1);
   AMRD_real_param(pfile,"u_l",&u_l,1);
   AMRD_real_param(pfile,"u_r",&u_r,1);
   AMRD_real_param(pfile,"ut_l",&ut_l,1);
   AMRD_real_param(pfile,"ut_r",&ut_r,1);
   AMRD_real_param(pfile,"ux_l",&ux_l,1);
   AMRD_real_param(pfile,"ux_r",&ux_r,1);
   AMRD_real_param(pfile,"uy_l",&uy_l,1);
   AMRD_real_param(pfile,"uy_r",&uy_r,1);
   AMRD_real_param(pfile,"uz_l",&uz_l,1);
   AMRD_real_param(pfile,"uz_r",&uz_r,1);
   AMRD_int_param(pfile,"eos_flag",&eos_flag,1);
   AMRD_real_param(pfile,"Ye",&Ye,1);
   AMRD_real_param(pfile,"T_l",&T_l,1);
   AMRD_real_param(pfile,"T_r",&T_r,1);
   AMRD_real_param(pfile,"rho0",&rho0,1);
   AMRD_real_param_v(pfile,"TRE_max_hydro",&TRE_max_hydro,&num_TRE_max_hydro); 
   AMRD_int_param(pfile,"riemann_prob",&riemann_prob,1); 
   AMRD_int_param(pfile,"bc_type",&bc_type,1); 
   AMRD_int_param(pfile,"N_bound",&N_bound,1); 
   AMRD_real_param(pfile,"lorentz_max",&lorentz_max,1);

   if (my_rank==0) printf("===================================================================\n");
   return;
}

//=============================================================================
// Sets all variables to their 'zero' values:
//=============================================================================
void srhd_AMRH_var_clear(void)
{
   ldptr();

   zero_c(d_n);  zero_c(d_np1);
   zero_c(Sx_n); zero_c(Sx_np1);
   zero_c(Sy_n); zero_c(Sy_np1);
   zero_c(Sz_n); zero_c(Sz_np1);
      
   return;
}

//=============================================================================
// Sets fcs variables to their 'zero' values.
// typeflag = 0:  corresponds to negative sign flag in mask.  zero cells of type A
// typeflag = 1:  corresponds to positive sign flag in mask.  zero cells of type B
//=============================================================================
void srhd_fcs_var_clear(int type_flag, int *ifc_mask)
{
   ldptr();
   zero_fcs_vars_(&type_flag,ifc_mask,d_fcs,Sx_fcs,Sy_fcs,Sz_fcs,E_fcs,&Nx,&Ny,&Nz);
   return;
}

//=============================================================================
// Initial data for flat space hydro tests
// Problem type:
// 0 -- linear shock tube 
// 1 -- 2d Riemann problem
// 2 -- 2d cylindrical blast wave
// 3 -- shock at 45 degrees.
//=============================================================================
void srhd_free_data(void)
{
   int i;

   ldptr();

   if (eos_flag == 1) shen_initialize_();

   init_riemann_(&gamma,d_n,Sx_n,Sy_n,Sz_n,E_n,rho,vx,vy,vz,
		 u,P,T,cs,x_c,y_c,z_c,&Nx,&Ny,&Nz,
		 &x0,&rho_l,&rho_r,&u_l,&u_r,&ut_l,&ux_l,&uy_l,&uz_l,
		 &ut_r,&ux_r,&uy_r,&uz_r,&T_r,&T_l,&Ye,&rho0,&eos_flag,
		 &riemann_prob,&lorentz_max);

   for (i=0; i<size_c; i++) {
     d_fcs[i]  = 0;
     Sx_fcs[i] = 0;
     Sy_fcs[i] = 0;
     Sz_fcs[i] = 0;
     E_fcs[i]  = 0;
   }
   
   return;
}  

void srhd_init_rest(void)
{
   ldptr();
   
   init_rest_(d_n,d_np1,Sx_n,Sx_np1,Sy_n,Sy_np1,Sz_n,Sz_np1,E_n,E_np1,
	      &Nx,&Ny,&Nz);

   return;
}  

//=============================================================================
// Initial constraint data --- called after each MG iteration:
//=============================================================================
void srhd_t0_cnst_data(void)
{
   return;
}

//=============================================================================
// Returns some norm of the residual for the evolution variables,
// called after an evolution iteration:
//=============================================================================
double srhd_evo_residual(void)
{
   ldptr();

   return 0.0;
}

//=============================================================================
// Returns some norm of the residual for the MG variables, *AND* 
// stores the point-wise residual in "f_res" for each MG variable "f" (for
// computing new RHS's)
//=============================================================================
double srhd_MG_residual(void)
{
   return 0.0;
}

//=============================================================================
// Performs 1 iteration of the evolution equations 
//=============================================================================
void srhd_evolve(int iter, int *ifc_mask)
{
   ldptr();
   int i,j,k,ind;
   if (iter == 1) {
     hydro_1step_(d_np1,Sx_np1,Sy_np1,Sz_np1,E_np1,d_p,Sx_p,Sy_p,Sz_p,E_p,
		  rho,vx,vy,vz,u,P,T,cs,x_c,y_c,z_c,&gamma,&limiter,
		  &dto2,&Ye,&rho0,&Nx,&Ny,&Nz,ifc_mask,d_fcs,Sx_fcs,Sy_fcs,
		  Sz_fcs,E_fcs,&iter,&eos_flag,phys_bdy,&bc_type,&N_bound,
		  &lorentz_max);

     for (j=0; j<size_c; j++) {
       d_np1[j]  = d_n[j]  + dto2*d_p[j];
       Sx_np1[j] = Sx_n[j] + dto2*Sx_p[j];
       Sy_np1[j] = Sy_n[j] + dto2*Sy_p[j];
       Sz_np1[j] = Sz_n[j] + dto2*Sz_p[j];
       E_np1[j]  = E_n[j]  + dto2*E_p[j];
     }

   } else if (iter == 2) {
     // Perform full step.
     hydro_1step_(d_np1,Sx_np1,Sy_np1,Sz_np1,E_np1,d_p,Sx_p,Sy_p,Sz_p,E_p,
		  rho,vx,vy,vz,u,P,T,cs,x_c,y_c,z_c,&gamma,&limiter,
		  &dt,&Ye,&rho0,&Nx,&Ny,&Nz,ifc_mask,d_fcs,Sx_fcs,Sy_fcs,
		  Sz_fcs,E_fcs,&iter,&eos_flag,phys_bdy,&bc_type,&N_bound,
		  &lorentz_max);

     for (j=0; j<size_c; j++) {
       d_np1[j]  = d_n[j]  + dt*d_p[j];
       Sx_np1[j] = Sx_n[j] + dt*Sx_p[j];
       Sy_np1[j] = Sy_n[j] + dt*Sy_p[j];
       Sz_np1[j] = Sz_n[j] + dt*Sz_p[j];
       E_np1[j]  = E_n[j]  + dt*E_p[j];

     }
     // Perform final inversion and boundary conditions.
     cons_bc_(&bc_type,d_np1,E_np1,Sx_np1,Sy_np1,Sz_np1,&Nx,&Ny,&Nz,&N_bound,phys_bdy);

     // First call primitive variable inversion
     if (eos_flag == 0) {
       find_primitives_(d_np1,Sx_np1,Sy_np1,Sz_np1,E_np1,rho,vx,vy,vz,u,&gamma,&Nx,&Ny,&Nz);
     } else if (eos_flag == 1) {
       find_primitives_table_(d_np1,Sx_np1,Sy_np1,Sz_np1,E_np1,rho,vx,vy,vz,u,P,T,cs,&Ye,&rho0,&Nx,&Ny,&Nz);
     }
   }

   return;
}

//=============================================================================
// Performs 1 relaxation sweep of the MG equations, and returns an estimate
// of the norm of the residual.
//=============================================================================
double srhd_MG_relax(void)
{
   return 0.0;
}

//=============================================================================
// Computes the differential operator L acting on the current grid,
// in a region specified by the grid function "mask". Stores the result
// in "f_lop" for each MG variable "f"
//=============================================================================
void srhd_L_op(void)
{
   return;
}

//=============================================================================
// Calculations prior to saving info to disk.
//=============================================================================
void srhd_pre_io_calc(void)
{
   return;
}

//=============================================================================
// Called after calculating the TRE for all variables
//=============================================================================
void srhd_scale_tre(void)
{
  ldptr();
  find_tre_hydro_(d_tre,rho,&Nx,&Ny,&Nz);

  int i,j,k, ind;
  
  for (i=0; i<Nx; i++) {
    for (j=0; j<Ny; j++) {
      for (k=0; k<Nz; k++) {
	ind = i + Nx * j + Nx * Ny * k;

	// Branson:  Note, I needed information about the grid level I'm on.
	// This is probably against the spirit of the application subroutines
	// doing the identical operation on every grid.  

	if (d_tre[ind] > TRE_max_hydro[loc_lev-1]) d_tre[ind] = 1.0;
      }
    }
  }
  
}

//=============================================================================
// no post-regrid/tstep stuff
//=============================================================================
void srhd_post_regrid(void)
{
}

void srhd_post_tstep(int L)
{
  return;
}

//=============================================================================
void srhd_pre_tstep(int L)
{

  int valid;
  static int local_first=1;
  static int step_counter = 1;
  char name[256];
  static FILE *out1;
  static FILE *out2;
  char out1_name[256];
  char out2_name[256];
  char buffer[20];
  int mpi_size;
  int Lf,Lc,write_lev,ict;
  int ind, i, j, k;
  real m, lm,ct;
  real rho_c, lrho_c;

  Lf=PAMR_get_max_lev(PAMR_AMRH);
  Lc=PAMR_get_min_lev(PAMR_AMRH);
  
  if (Lf==1) write_lev = Lf;
  if (Lf>1)  write_lev = Lc + 1; // not the shadow level

  // Next, write some diagnostics to a file.
  if (L==write_lev) {
    
    if (my_rank == 0) {
      ct=PAMR_get_time(L);
      ict = ct * 1000;
      //sprintf(out1_name,"%s.dat",itoa(ict));
      sprintf(out1_name,"out.dat");
      if (!(out1=fopen(out1_name,"w"))) {
	AMRD_stop("srhd_post_tstep: error opening out.dat file\n","");
      }

      valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
      
      while(valid)
	{
	  ldptr();
	  valid=PAMR_next_g();
	  
	  for (i=0; i<Nx; i++) {
	    for (j=0; j<Ny; j++) {
	      for (k=0; k<Nz; k++) {
		ind = i + Nx * j + Nx * Ny * k;
		if (j==2 && k==2) fprintf(out1,"%15.6E %15.6E %15.6E %15.6E \n",x_c[i],rho[ind],P[ind],vx[ind]);
	      }
	    }
	  }
	}

      fflush(out1);

      if (AMRD_num_f_tre_vars > 0) {
	// Have a look at your refinement grid function.
	sprintf(out2_name,"tre.dat");
	if (!(out2=fopen(out2_name,"w"))) {
	  AMRD_stop("srhd_post_tstep: error opening tre.dat file\n","");
	}
	valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
	
	while(valid)
	  {
	    ldptr();
	    find_tre_hydro_(d_tre,rho,&Nx,&Ny,&Nz);
	    
	    valid=PAMR_next_g();
	    
	    for (i=0; i<Nx; i++) {
	      for (j=0; j<Ny; j++) {
		for (k=0; k<Nz; k++) {
		  ind = i + Nx * j + Nx * Ny * k;
		  if (j==2 && k==2) fprintf(out2,"%15.6E %15.6E\n",x_c[i],d_tre[ind]);
		}
	      }
	    }
	  }

	fflush(out2);
      }
    }
      
    valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
    while(valid)
      {
	ldptr();
	total_mass_(&lm,d_n,x,y,z,&Nx,&Ny,&Nz,ghost_width,&my_rank);
	valid=PAMR_next_g();
      }
    
    MPI_Allreduce(&lm,&m,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if (my_rank == 0) ct=PAMR_get_time(L);
  }

  step_counter++;
}

void srhd_find_primitives()
{
  ldptr();
  if (eos_flag == 0) {
    find_primitives_(d_n,Sx_n,Sy_n,Sz_n,E_n,rho,vx,vy,vz,u,&gamma,&Nx,&Ny,&Nz);
  } else if (eos_flag == 1) {
    find_primitives_table_(d_n,Sx_n,Sy_n,Sz_n,E_n,rho,vx,vy,vz,u,P,T,cs,&Ye,&rho0,&Nx,&Ny,&Nz);
  }
}

void srhd_fill_ex_mask(double *mask, double *mask_c, int dim, int *shape, int *shape_c, double *bbox, double excised)
{
}

//=============================================================================
int main(int argc, char **argv)
{
  amrd_set_app_fcs_var_clear_hook(srhd_fcs_var_clear);
  amrd_set_app_flux_correct_hook(srhd_flux_correct);
  amrd_set_app_post_flux_correct_hook(srhd_post_flux_correct);
  amrd_set_app_pre_tstep_hook(srhd_pre_tstep);
  amrd(argc,argv,&srhd_id,&srhd_var_pre_init,
       &srhd_var_post_init, &srhd_AMRH_var_clear,
       &srhd_free_data, &srhd_t0_cnst_data,
       &srhd_evo_residual, &srhd_MG_residual,
       &srhd_evolve, &srhd_MG_relax, &srhd_L_op, 
       &srhd_pre_io_calc,&srhd_scale_tre,
       &srhd_post_regrid,&srhd_post_tstep,
       &srhd_fill_ex_mask,0);
}

// A hook function to apply the flux correction.

void srhd_flux_correct(void)
{
  ldptr();
  
  apply_fc_(d_n,Sx_n,Sy_n,Sz_n,E_n,d_fcs,Sx_fcs,Sy_fcs,Sz_fcs,E_fcs,&Nx,&Ny,&Nz);
  
  return;
}  

// A hook function to take care of any unfinished business
// following the flux correction step.

void srhd_post_flux_correct(void)
{
  ldptr();
  if (eos_flag==0) {
    find_primitives_(d_n,Sx_n,Sy_n,Sz_n,E_n,rho,vx,vy,vz,u,&gamma,&Nx,&Ny,&Nz);
  } else if (eos_flag==1) {
    find_primitives_table_(d_n,Sx_n,Sy_n,Sz_n,E_n,rho,vx,vy,vz,u,P,T,cs,&Ye,&rho0,&Nx,&Ny,&Nz);
  }
}

//--------------------------------------------------------------------
// This is an integer to string converter I found on the internet.
//--------------------------------------------------------------------
char *itoa(int value)
{
int count,                   /* number of characters in string       */
    i,                       /* loop control variable                */
    sign;                    /* determine if the value is negative   */
char *ptr,                   /* temporary pointer, index into string */
     *string,                /* return value                         */
     *temp;                  /* temporary string array               */

count = 0;
if ((sign = value) < 0)      /* assign value to sign, if negative    */
   {                         /* keep track and invert value          */
   value = -value;
   count++;                  /* increment count                      */
   }

/* allocate INTSIZE plus 2 bytes (sign and NULL)                     */
temp = (char *) malloc(INTSIZE + 2);
if (temp == NULL)
   {
   return(NULL);
   }
memset(temp,'\0', INTSIZE + 2);

string = (char *) malloc(INTSIZE + 2);
if (string == NULL)
   {
   return(NULL);
   }
memset(string,'\0', INTSIZE + 2);
ptr = string;                /* set temporary ptr to string          */

/*--------------------------------------------------------------------+
| NOTE: This process reverses the order of an integer, ie:            |
|       value = -1234 equates to: char [4321-]                        |
|       Reorder the values using for {} loop below                    |
+--------------------------------------------------------------------*/
do {
   *temp++ = value % 10 + '0';   /* obtain modulus and or with '0'   */
   count++;                      /* increment count, track iterations*/
   }  while (( value /= 10) >0);

if (sign < 0)                /* add '-' when sign is negative        */
   *temp++ = '-';

*temp-- = '\0';              /* ensure null terminated and point     */
                             /* to last char in array                */

/*--------------------------------------------------------------------+
| reorder the resulting char *string:                                 |
| temp - points to the last char in the temporary array               |
| ptr  - points to the first element in the string array              |
+--------------------------------------------------------------------*/
for (i = 0; i < count; i++, temp--, ptr++)
   {
   memcpy(ptr,temp,sizeof(char));
   }

return(string);
}

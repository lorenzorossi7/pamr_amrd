//=============================================================================
// application functions/variables for GRHD in one dimension
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

int limiter;

double *rho, *u, *v, *T, *P, *cs, *jr, *rhops, *ham_source;
double *d, *Sr, *tau;
double *d_p, *Sr_p, *tau_p, *psi_p_v;
double *d_fcs, *Sr_fcs, *tau_fcs;
double *d_tre;
double *psi_v, *psi_lop_v, *psi_res_v, *psi_rhs_v;
double *psi_extrap_tm1_v, *psi_extrap_tm2_v;
double *alpha_v, *alpha_lop_v, *alpha_res_v, *alpha_rhs_v;
double *alpha_extrap_tm1_v, *alpha_extrap_tm2_v;
double *beta_v, *beta_lop_v, *beta_res_v, *beta_rhs_v;
double *beta_extrap_tm1_v, *beta_extrap_tm2_v;
double *psi0_v, *psi0_lop_v, *psi0_res_v, *psi0_rhs_v;
double *jr_v, *rhops_v, *ham_source_v, *alpha, *beta, *psi;
double *jr_v_mg, *rhops_v_mg, *ham_source_v_mg;
double *ham_v, *tau_v;

double *d_n, *Sr_n, *tau_n;
double *psi_v_n, *alpha_v_n, *beta_v_n;
double *d_np1, *Sr_np1, *tau_np1;
double *psi_v_np1, *alpha_v_np1, *beta_v_np1;
double *mask, *mask_mg;

double *x, *x_c;
int shape[3],shape_c[3],ghost_width[6],Nx,phys_bdy[6],size,size_c;
double base_bbox[6],bbox[6],dx,dt,dto2;
int g_L;
double alpha_out, psi0_out;
// Note that these are not independently specifiable parameters, 
// but are to be calculated using atm_frac.
double rho_atm, p_atm, u_atm, t_atm;

int rho_gfn, u_gfn, v_gfn, T_gfn, P_gfn, cs_gfn;
int d_gfn, Sr_gfn, tau_gfn;
int d_p_gfn, Sr_p_gfn, tau_p_gfn, psi_p_v_gfn;
int psi_v_gfn, psi_lop_v_gfn, psi_res_v_gfn, psi_rhs_v_gfn;
int alpha_v_gfn, alpha_lop_v_gfn, alpha_res_v_gfn, alpha_rhs_v_gfn;
int beta_v_gfn, beta_lop_v_gfn, beta_res_v_gfn, beta_rhs_v_gfn;
int psi0_v_gfn, psi0_lop_v_gfn, psi0_res_v_gfn, psi0_rhs_v_gfn;
int psi_gfn, alpha_gfn, beta_gfn;
int d_fcs_gfn, Sr_fcs_gfn, tau_fcs_gfn;
int d_tre_gfn;

int d_n_gfn, Sr_n_gfn, tau_n_gfn, psi_v_n_gfn;
int alpha_v_n_gfn, beta_v_n_gfn;
int d_np1_gfn, Sr_np1_gfn, tau_np1_gfn, psi_v_np1_gfn;
int alpha_v_np1_gfn, beta_v_np1_gfn;
int psi_extrap_tm1_v_gfn, psi_extrap_tm2_v_gfn;
int alpha_extrap_tm1_v_gfn, alpha_extrap_tm2_v_gfn;
int beta_extrap_tm1_v_gfn, beta_extrap_tm2_v_gfn;
int mask_gfn, mask_mg_gfn;
int jr_gfn, jr_v_gfn, jr_v_mg_gfn;
int rhops_gfn, rhops_v_gfn, rhops_v_mg_gfn;
int ham_v_gfn, tau_v_gfn;
int ham_source_gfn, ham_source_v_gfn, ham_source_v_mg_gfn;

int max_lev;
int init_depth;
int dim;
double t0;
int ID_DV_trace;
int eos_flag;
double Ye;  // For now, we'll set this to a constant.
double rho0; // The density scale to pass into the eos.
double temperature;
double deltat1, deltat2;
int const_ev; // Whether to do constrained evolution.

// Parameters for Initial Data
double rho0_c;
double atm_frac;
double gamma;
double x0;
double rho_l, rho_r;
double u_l, u_r;
double ut_l,ur_l;
double ut_r,ur_r;
double T_l, T_r;
double p_deplete;
double lorentz_max;
double v_pert;
int ascii_dump_int;

void grhd_1d_free_data(void);
void grhd_1d_init_rest(void);
void grhd_1d_find_primitives(void);
void grhd_1d_fcs_var_clear(int type_flag, int *ifc_mask);
void grhd_1d_flux_correct(void);
void grhd_1d_post_flux_correct(void);
void grhd_1d_find_tre_hydro(void);
void grhd_1d_pre_tstep(int L);
void grhd_1d_elliptic_vars_t0_init(void);
void grhd_1d_t0_cnst_data(void);

//=============================================================================
// call after variables have been defined
//=============================================================================
void set_gfns(void)
{
    if ((d_gfn=PAMR_get_gfn("d",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((d_np1_gfn=PAMR_get_gfn("d",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((d_n_gfn=PAMR_get_gfn("d",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((Sr_gfn=PAMR_get_gfn("Sr",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((Sr_np1_gfn=PAMR_get_gfn("Sr",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Sr_n_gfn=PAMR_get_gfn("Sr",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((tau_gfn=PAMR_get_gfn("tau",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((tau_np1_gfn=PAMR_get_gfn("tau",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((tau_n_gfn=PAMR_get_gfn("tau",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((d_fcs_gfn   = PAMR_get_gfn("d_fcs"  , PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((Sr_fcs_gfn  = PAMR_get_gfn("Sr_fcs" , PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((tau_fcs_gfn = PAMR_get_gfn("tau_fcs", PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 

    if ((psi_v_gfn=PAMR_get_gfn("psi_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_v_np1_gfn=PAMR_get_gfn("psi_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_v_n_gfn=PAMR_get_gfn("psi_v",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((alpha_v_gfn=PAMR_get_gfn("alpha_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((alpha_v_np1_gfn=PAMR_get_gfn("alpha_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((alpha_v_n_gfn=PAMR_get_gfn("alpha_v",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((beta_v_gfn=PAMR_get_gfn("beta_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((beta_v_np1_gfn=PAMR_get_gfn("beta_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((beta_v_n_gfn=PAMR_get_gfn("beta_v",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((psi_extrap_tm1_v_gfn=PAMR_get_gfn("psi_extrap_tm1_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_extrap_tm2_v_gfn=PAMR_get_gfn("psi_extrap_tm2_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((alpha_extrap_tm1_v_gfn=PAMR_get_gfn("alpha_extrap_tm1_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((alpha_extrap_tm2_v_gfn=PAMR_get_gfn("alpha_extrap_tm2_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((beta_extrap_tm1_v_gfn=PAMR_get_gfn("beta_extrap_tm1_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((beta_extrap_tm2_v_gfn=PAMR_get_gfn("beta_extrap_tm2_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((psi_res_v_gfn=PAMR_get_gfn("psi_res_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_lop_v_gfn=PAMR_get_gfn("psi_lop_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_rhs_v_gfn=PAMR_get_gfn("psi_rhs_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((alpha_res_v_gfn=PAMR_get_gfn("alpha_res_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((alpha_lop_v_gfn=PAMR_get_gfn("alpha_lop_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((alpha_rhs_v_gfn=PAMR_get_gfn("alpha_rhs_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((beta_res_v_gfn=PAMR_get_gfn("beta_res_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((beta_lop_v_gfn=PAMR_get_gfn("beta_lop_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((beta_rhs_v_gfn=PAMR_get_gfn("beta_rhs_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mask_mg_gfn=PAMR_get_gfn("cmask_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mask_gfn=PAMR_get_gfn("cmask_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((psi0_v_gfn=PAMR_get_gfn("psi0_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((psi0_res_v_gfn=PAMR_get_gfn("psi0_res_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((psi0_lop_v_gfn=PAMR_get_gfn("psi0_lop_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((psi0_rhs_v_gfn=PAMR_get_gfn("psi0_rhs_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);

    // Some work variables
    if ((d_p_gfn=PAMR_get_gfn("d_p",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((Sr_p_gfn=PAMR_get_gfn("Sr_p",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((tau_p_gfn=PAMR_get_gfn("tau_p",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((psi_p_v_gfn=PAMR_get_gfn("psi_p_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 

    if ((rho_gfn=PAMR_get_gfn("rho",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((v_gfn=PAMR_get_gfn("v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((u_gfn=PAMR_get_gfn("u",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((T_gfn=PAMR_get_gfn("T",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((P_gfn=PAMR_get_gfn("P",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((cs_gfn=PAMR_get_gfn("cs",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((jr_gfn=PAMR_get_gfn("jr",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((jr_v_gfn=PAMR_get_gfn("jr_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((jr_v_mg_gfn=PAMR_get_gfn("jr_v",PAMR_MGH,0) )<0) AMRD_stop("set_gnfs error",0);
    if ((rhops_gfn=PAMR_get_gfn("rhops",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((rhops_v_gfn=PAMR_get_gfn("rhops_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((rhops_v_mg_gfn=PAMR_get_gfn("rhops_v",PAMR_MGH,0) )<0) AMRD_stop("set_gnfs error",0);

    if ((ham_source_gfn=PAMR_get_gfn("ham_source",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((ham_source_v_gfn=PAMR_get_gfn("ham_source_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((ham_source_v_mg_gfn=PAMR_get_gfn("ham_source_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0); 

    if ((alpha_gfn=PAMR_get_gfn("alpha",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((beta_gfn=PAMR_get_gfn("beta",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((psi_gfn=PAMR_get_gfn("psi",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 

    if ((ham_v_gfn=PAMR_get_gfn("ham_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((tau_v_gfn=PAMR_get_gfn("tau_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    
    //    if ((d_tre_gfn=PAMR_get_gfn("d_tre",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 

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

   d    = gfs[d_gfn-1];
   d_fcs= gfs[d_fcs_gfn-1];
   d_n  = gfs[d_n_gfn-1];
   d_np1= gfs[d_np1_gfn-1];

   Sr     = gfs[Sr_gfn-1];
   Sr_fcs = gfs[Sr_fcs_gfn-1];
   Sr_n   = gfs[Sr_n_gfn-1];
   Sr_np1 = gfs[Sr_np1_gfn-1];

   tau     = gfs[tau_gfn-1];
   tau_fcs = gfs[tau_fcs_gfn-1];
   tau_n   = gfs[tau_n_gfn-1];
   tau_np1 = gfs[tau_np1_gfn-1];

   psi_v              = gfs[psi_v_gfn-1];
   psi_v_n            = gfs[psi_v_n_gfn-1];
   psi_v_np1          = gfs[psi_v_np1_gfn-1];
   psi_res_v          = gfs[psi_res_v_gfn-1];
   psi_lop_v          = gfs[psi_lop_v_gfn-1];
   psi_rhs_v          = gfs[psi_rhs_v_gfn-1];
   psi_extrap_tm1_v   = gfs[psi_extrap_tm1_v_gfn-1];
   psi_extrap_tm2_v   = gfs[psi_extrap_tm2_v_gfn-1];

   alpha_v            = gfs[alpha_v_gfn-1];
   alpha_v_n          = gfs[alpha_v_n_gfn-1];
   alpha_v_np1        = gfs[alpha_v_np1_gfn-1];
   alpha_res_v        = gfs[alpha_res_v_gfn-1];
   alpha_lop_v        = gfs[alpha_lop_v_gfn-1];
   alpha_rhs_v        = gfs[alpha_rhs_v_gfn-1];
   alpha_extrap_tm1_v = gfs[alpha_extrap_tm1_v_gfn-1];
   alpha_extrap_tm2_v = gfs[alpha_extrap_tm2_v_gfn-1];

   beta_v             = gfs[beta_v_gfn-1];
   beta_v_n           = gfs[beta_v_n_gfn-1];
   beta_v_np1         = gfs[beta_v_np1_gfn-1];
   beta_res_v         = gfs[beta_res_v_gfn-1];
   beta_lop_v         = gfs[beta_lop_v_gfn-1];
   beta_rhs_v         = gfs[beta_rhs_v_gfn-1];
   beta_extrap_tm1_v  = gfs[beta_extrap_tm1_v_gfn-1];
   beta_extrap_tm2_v  = gfs[beta_extrap_tm2_v_gfn-1];

   psi0_v             = gfs[psi0_v_gfn-1];
   psi0_res_v         = gfs[psi0_res_v_gfn-1];
   psi0_lop_v         = gfs[psi0_lop_v_gfn-1];
   psi0_rhs_v         = gfs[psi0_rhs_v_gfn-1];

   mask    = gfs[mask_gfn-1];
   mask_mg = gfs[mask_mg_gfn-1];

   d_p     = gfs[d_p_gfn-1];
   Sr_p    = gfs[Sr_p_gfn-1];
   tau_p   = gfs[tau_p_gfn-1];
   psi_p_v = gfs[psi_p_v_gfn-1];

   rho   = gfs[rho_gfn-1];
   v     = gfs[v_gfn-1];
   u     = gfs[u_gfn-1];
   T     = gfs[T_gfn-1];
   P     = gfs[P_gfn-1];
   cs    = gfs[cs_gfn-1];

   jr         = gfs[jr_gfn-1];
   jr_v       = gfs[jr_v_gfn-1];
   jr_v_mg    = gfs[jr_v_mg_gfn-1];
   rhops      = gfs[rhops_gfn-1];
   rhops_v    = gfs[rhops_v_gfn-1];
   rhops_v_mg = gfs[rhops_v_mg_gfn-1];

   ham_source      = gfs[ham_source_gfn-1];
   ham_source_v    = gfs[ham_source_v_gfn-1];
   ham_source_v_mg = gfs[ham_source_v_mg_gfn-1];

   alpha  = gfs[alpha_gfn-1];
   beta   = gfs[beta_gfn-1];
   psi    = gfs[psi_gfn-1];

   ham_v       = gfs[ham_v_gfn-1];
   tau_v       = gfs[tau_v_gfn-1];

   //   d_tre = gfs[d_tre_gfn-1];

   x=x0[0]; dx=x[1]-x[0];
   x_c=x0_c[0]; 
   PAMR_get_g_level(&lev);
   PAMR_get_dxdt(lev,dx0,&dt);
   dto2 =dt/2.0;

   if ((bbox[0]-base_bbox[0])<dx/2) phys_bdy[0]=1; else phys_bdy[0]=0;
   if ((base_bbox[1]-bbox[1])<dx/2) phys_bdy[1]=1; else phys_bdy[1]=0;

   Nx=shape_c[0];

   size  = shape[0];
   size_c=Nx;
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
// Returns 0 to use default mechanism, or is expected to calculate
// the correct initial hierarchy and return 1:
//=============================================================================
int grhd_1d_id(void)
{
  int valid;
  int Lf,L,Lc;
  int gnum;
  real base_length_x,base_length_y;
  real new_length_x,new_length_y;
  int glev[PAMR_MAX_LEVS];
  real gbbox[PAMR_MAX_LEVS*PAMR_MAX_DIM*2];
  int i,j;
  // Achtung.  This assumes that, if you want more than one level, 
  // you'll be specifying them yourself, a la FMR.

  if (max_lev > 1) {
    // Devise your boxes (for levels higher than the base level)
    gnum = max_lev - init_depth;
    if (gnum < 1) printf("problem inside grhd_1d_id\n");

    base_length_x = base_bbox[1]-base_bbox[0];

    for (i=1; i<=gnum; i++) {
      glev[i-1] = init_depth + i;
      if (i==1) {
	new_length_x = base_length_x / 2.0;
	gbbox[(i-1)*2*dim ] = base_bbox[0] + 0.25*base_length_x;
	gbbox[(i-1)*2*dim ] = base_bbox[0] + 0.15*base_length_x;
      } else {
	new_length_x = new_length_x / 2.0;
	gbbox[(i-1)*2*dim ] = gbbox[(i-2)*2*dim ] + 0.5*new_length_x;
	gbbox[(i-1)*2*dim+1] = gbbox[(i-1)*2*dim ] + new_length_x;
      }
    }

    // Compose the hierarchy
    t0 = 0.0;

    if (!(PAMR_compose_hierarchy(init_depth+1,max_lev,gnum,glev,gbbox,t0))) printf("ouch\n");
    // Fill grid functions on all levels with the initial data.
    Lf=PAMR_get_max_lev(PAMR_AMRH);
    Lc=PAMR_get_min_lev(PAMR_AMRH);

    // Branson messing around
    for (L=1; L<=Lf; L++) {
      valid = PAMR_init_s_iter(L,PAMR_AMRH,0);
      while(valid) {
	grhd_1d_free_data();
	valid = PAMR_next_g();
      }
      PAMR_sync(L,2,PAMR_AMRH,0);
    }

    // Branson messing around
    for (L=Lf; L>1; L--) {
      PAMR_inject(L,2,PAMR_AMRH);
      valid = PAMR_init_s_iter(L-1,PAMR_AMRH,0);
      while(valid) {
	//	grhd_1d_find_primitives();
	grhd_1d_post_flux_correct();
	valid = PAMR_next_g();
      }
    }

    // Dump the initial data stuff if desired.
    if (ID_DV_trace) {
      PAMR_save_gfn("d",PAMR_AMRH,2,-1,-1.0,"ID_","");
      PAMR_save_gfn("Sr",PAMR_AMRH,2,-1,-1.0,"ID_","");
      PAMR_save_gfn("tau",PAMR_AMRH,2,-1,-1.0,"ID_","");
    }

    // You need to do the stuff contained in init_rest
    for (L=1; L<=Lf; L++) {
      valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
      while(valid) {
	grhd_1d_init_rest();
	valid=PAMR_next_g();
      }
    }

    // Set the fc mask
    for (L=2; L<=Lf; L++) set_fc_mask(L);

    return 1;
  } else {
    return 0;
  }
}

//=============================================================================
// Sets custom parameters, variables, etc. Split up into two segments,
// one called before the pamr context is initialized and standard
// parameters are read, and the other afterwards
//=============================================================================
void grhd_1d_var_pre_init(char *pfile)
{
  // Read in some global parameters which are needed by grhd_1d_id
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

void grhd_1d_var_post_init(char *pfile)
{

   if (my_rank==0)
   {
      printf("===================================================================\n");
      printf("Reading GRHD_1D parameters:\n\n");
   }

   gamma = 1.333333333333333333;
   rho0_c = 0.1;
   atm_frac = 1.e-6;
   
   limiter = 1;
   x0 = 1.0;
   rho_l = 1.0;
   rho_r = 1.0;
   u_l = 1.0;
   u_r = 1.0;
   ut_l = 1.0;
   ut_r = 1.0;
   u_l = 0.0;
   u_r = 0.0;
   Ye = 0.1;               // Defaults to neutron rich matter.
   T_l = 0.0;
   T_r = 0.0;
   rho0 = 1.e14;
   ascii_dump_int = 1000;  // Defaults to a largish number.
   p_deplete = 1.0;        // Defaults to no pressure depletion.
   v_pert = 0.0;           // Defaults to zero veloctiy perturbation
   lorentz_max = 10.0;     // Defaults to kinda relativistic.
   eos_flag = 0;           // Defaults to polytrope.
   deltat1 = 0.1;          // Defaults to 10%
   deltat2 = 0.1;
   const_ev = 0;           // Defaults to unconstrained evolution.
  
   AMRD_real_param(pfile,"gamma",&gamma,1); 
   AMRD_real_param(pfile,"rho0_c",&rho0_c,1); 
   AMRD_real_param(pfile,"atm_frac",&atm_frac,1); 
   AMRD_int_param(pfile,"limiter",&limiter,1); 
   AMRD_real_param(pfile,"x0",&x0,1);
   AMRD_real_param(pfile,"rho_l",&rho_l,1);
   AMRD_real_param(pfile,"rho_r",&rho_r,1);
   AMRD_real_param(pfile,"u_l",&u_l,1);
   AMRD_real_param(pfile,"u_r",&u_r,1);
   AMRD_real_param(pfile,"ut_l",&ut_l,1);
   AMRD_real_param(pfile,"ut_r",&ut_r,1);
   AMRD_real_param(pfile,"ur_l",&ur_l,1);
   AMRD_real_param(pfile,"ur_r",&ur_r,1);
   AMRD_int_param(pfile,"eos_flag",&eos_flag,1);
   AMRD_real_param(pfile,"Ye",&Ye,1);
   AMRD_real_param(pfile,"temperature",&temperature,1);
   AMRD_real_param(pfile,"deltat1",&deltat1,1);
   AMRD_real_param(pfile,"deltat2",&deltat2,1);
   AMRD_real_param(pfile,"T_l",&T_l,1);
   AMRD_real_param(pfile,"T_r",&T_r,1);
   AMRD_real_param(pfile,"rho0",&rho0,1);
   AMRD_real_param(pfile,"p_deplete",&p_deplete,1);
   AMRD_real_param(pfile,"v_pert",&v_pert,1);
   AMRD_real_param(pfile,"lorentz_max",&lorentz_max,1);
   AMRD_int_param(pfile,"ascii_dump_int",&ascii_dump_int,1);
   AMRD_int_param(pfile,"const_ev",&const_ev,1);

   if (my_rank==0) printf("===================================================================\n");
   return;
}

//=============================================================================
// Sets all variables to their 'zero' values:
//=============================================================================
void grhd_1d_AMRH_var_clear(void)
{
   ldptr();
   zero_c(Sr_n); zero_c(Sr_np1);
   return;
}

//=============================================================================
// Sets fcs variables to their 'zero' values.
// typeflag = 0:  corresponds to negative sign flag in mask.  zero cells of type A
// typeflag = 1:  corresponds to positive sign flag in mask.  zero cells of type B
//=============================================================================
void grhd_1d_fcs_var_clear(int type_flag, int *ifc_mask)
{
   ldptr();
   zero_fcs_vars_(&type_flag,ifc_mask,d_fcs,Sr_fcs,tau_fcs,&Nx);
   return;
}

//=============================================================================
// Initial data for Newtonian polytrope: (at tn=2)
//=============================================================================
void grhd_1d_free_data(void)
{
   int i;

   ldptr();

   if (eos_flag == 1) { 
     shen_initialize_();
   }

   id_tov_(&gamma,&rho0_c,&atm_frac,d_n,Sr_n,tau_n,rho,v,u,P,T,jr,jr_v,
	   rhops,rhops_v,tau_v,ham_source,ham_source_v,
	   alpha_v_n,beta_v_n,psi_v_n,alpha,beta,psi,x,x_c,&Nx,
	   ham_v,mask,phys_bdy,ghost_width,&p_deplete,&Ye,
	   &temperature,&rho0,&eos_flag,&alpha_out,&rho_atm,
	   &p_atm,&u_atm,&t_atm,&deltat1,&deltat2,&v_pert,&lorentz_max,&psi0_out);

   // Testing.
   if (eos_flag == 0) {
    find_primitives_(d_n,Sr_n,tau_n,rho,v,u,jr,rhops,ham_source,&gamma,alpha,beta,psi,x_c,
		     &Nx,&rho_atm,&lorentz_max);
   } else if (eos_flag==1) {
     find_primitives_table_(d_n,Sr_n,tau_n,rho,v,u,jr,rhops,ham_source,P,T,cs,&Ye,&rho0,
			    alpha,beta,psi,x_c,&Nx,&deltat1,&deltat2,&rho_atm,
			    &p_atm,&u_atm,&t_atm,&lorentz_max);
   }

   return;
}  

void grhd_1d_init_rest(void)
{
   ldptr();
   
   init_rest_(d_n,d_np1,Sr_n,Sr_np1,tau_n,tau_np1,&Nx);

   return;
}  

//=============================================================================
// Returns some norm of the residual for the evolution variables,
// called after an evolution iteration:
//=============================================================================
double grhd_1d_evo_residual(void)
{
   ldptr();

   return 0.0;
}

//=============================================================================
// Returns some norm of the residual for the MG variables, *AND* 
// stores the point-wise residual in "f_res" for each MG variable "f" (for
// computing new RHS's)
//=============================================================================

//=============================================================================
// Returns some norm of the residual for the MG variables, *AND* 
// stores the point-wise residual in "f_res" for each MG variable "f" (for
// computing new RHS's)
//=============================================================================
double grhd_1d_MG_residual(void)
{
   double norm;
   double norm1, norm2;

   ldptr();

   if (const_ev == 0) {
     residual_(alpha_v,alpha_res_v,alpha_rhs_v,beta_v,beta_res_v,beta_rhs_v,psi_v,
	       jr_v_mg,rhops_v_mg,mask_mg,x,&norm1,&Nx,&alpha_out,phys_bdy);
   } else {
     residual_const_ev_(psi_v,psi_res_v,psi_rhs_v,alpha_v,alpha_res_v,alpha_rhs_v,
			beta_v,beta_res_v,beta_rhs_v,jr_v_mg,rhops_v_mg,ham_source_v_mg,
			mask_mg,x,&norm1,&Nx,&alpha_out,&psi0_out,phys_bdy);
   }

   norm2 = 0.0;
   if (psi0_v && const_ev == 0) {
     residual_t0_(psi0_v,psi0_res_v,psi0_rhs_v,alpha_v,beta_v,
		  ham_source_v_mg,mask_mg,x,&norm2,&Nx,&psi0_out,phys_bdy);
   }

   norm = sqrt(norm1*norm1 + norm2*norm2);

   return norm;
}

//=============================================================================
// Performs 1 relaxation sweep of the MG equations, and returns an estimate
// of the norm of the residual.
//=============================================================================
double grhd_1d_MG_relax(void)
{
   double norm;
   double norm1, norm2;
   ldptr();

   if (const_ev == 0) {
     relax_(alpha_v,alpha_rhs_v,beta_v,beta_rhs_v,psi_v,jr_v_mg,rhops_v_mg,mask_mg,
	    phys_bdy,x,&norm1,&Nx,ghost_width,&alpha_out);
   } else {
     relax_const_ev_(psi_v,psi_rhs_v,alpha_v,alpha_rhs_v,beta_v,beta_rhs_v,jr_v_mg,
		     rhops_v_mg,ham_source_v_mg,mask_mg,phys_bdy,x,&norm1,&Nx,
		     ghost_width,&alpha_out,&psi0_out);
   }

   norm2 = 0.0;
   if (psi0_v && const_ev == 0) {
     relax_t0_(psi0_v,psi0_rhs_v,alpha_v,beta_v,ham_source_v_mg,mask_mg,
	       phys_bdy,x,&norm2,&Nx,ghost_width,&psi0_out);
   }

   norm = sqrt(norm1*norm1+norm2*norm2);

   return norm;
}

//=============================================================================
// Computes the differential operator L acting on the current grid,
// in a region specified by the grid function "mask". Stores the result
// in "f_lop" for each MG variable "f"
//=============================================================================
void grhd_1d_L_op(void)
{
   ldptr();

   if (const_ev == 0) {
     lop_(alpha_v,alpha_lop_v,beta_v,beta_lop_v,psi_v,jr_v_mg,rhops_v_mg,mask_mg,
	  x,&Nx,&alpha_out,phys_bdy);
   } else {
     lop_const_ev_(psi_v,psi_lop_v,alpha_v,alpha_lop_v,beta_v,beta_lop_v,jr_v_mg,
		   rhops_v_mg,ham_source_v_mg,mask_mg,x,&Nx,&alpha_out,&psi0_out,
		   phys_bdy);
   }

   if (psi0_v && const_ev==0) {
     lop_t0_(psi0_v,psi0_lop_v,alpha_v,beta_v,ham_source_v_mg,mask_mg,x,&Nx,
	     &psi0_out,phys_bdy);
   }

   return;
}


void grhd_1d_evolve(int iter, int *ifc_mask)
{
   ldptr();
   int j;

   // Call local v_to_c to transfer interpolated metric quantities
   // to cell centers.
   PAMR_v_to_c_local(1,PAMR_AMRH);

   if (iter == 1) {
     if (const_ev == 0) {
       psi_1step_(psi_v_n,psi_p_v,beta_v_n,x,&Nx,phys_bdy,&iter);
     }
     hydro_1step_(d_np1,Sr_np1,tau_np1,d_p,Sr_p,tau_p,
		  rho,v,u,P,T,cs,jr,rhops,ham_source,x,x_c,&gamma,&limiter,
		  &Ye,&dto2,&rho0,&Nx,
		  psi_v_np1,psi,alpha_v_np1,alpha,beta_v_np1,beta,
		  ifc_mask,d_fcs,Sr_fcs,tau_fcs,&iter,&eos_flag,phys_bdy,
		  &deltat1,&deltat2,&rho_atm,&p_atm,&u_atm,&t_atm,&lorentz_max);
     for (j=0; j<size_c; j++) {
       d_np1[j]   = d_n[j]   + dto2 * d_p[j];
       Sr_np1[j]  = Sr_n[j]  + dto2 * Sr_p[j];
       tau_np1[j] = tau_n[j] + dto2 * tau_p[j];
     }

     if (const_ev == 0) {
       for (j=0; j<size; j++) {
	 psi_v_np1[j] = psi_v_n[j] + dto2 * psi_p_v[j];
       }
     }
   } else if (iter==2) {
     // Perform full step.
     if (const_ev == 0) {
       psi_1step_(psi_v_np1,psi_p_v,beta_v_np1,x,&Nx,phys_bdy,&iter);
     }

     hydro_1step_(d_np1,Sr_np1,tau_np1,d_p,Sr_p,tau_p,
		  rho,v,u,P,T,cs,jr,rhops,ham_source,x,x_c,&gamma,&limiter,
		  &Ye,&dt,&rho0,&Nx,
		  psi_v_np1,psi,alpha_v_np1,alpha,beta_v_np1,beta,
		  ifc_mask,d_fcs,Sr_fcs,tau_fcs,&iter,&eos_flag,phys_bdy,
		  &deltat1,&deltat2,&rho_atm,&p_atm,&u_atm,&t_atm,&lorentz_max);


     for (j=0; j<size_c; j++) {
       d_np1[j]   = d_n[j]   + dt * d_p[j];
       Sr_np1[j]  = Sr_n[j]  + dt * Sr_p[j];
       tau_np1[j] = tau_n[j] + dt * tau_p[j];
     }

     if (const_ev == 0) {
       for (j=0; j<size; j++) {
	 psi_v_np1[j] = psi_v_n[j] + dt * psi_p_v[j];
       }
     }

   } else if (iter == 3) {
     // Final inversion and boundary conditions.
     if (const_ev == 0) {
       psi_1step_(psi_v_np1,psi_p_v,beta_v_np1,x,&Nx,phys_bdy,&iter);
     }
     
     hydro_1step_(d_np1,Sr_np1,tau_np1,d_p,Sr_p,tau_p,
		  rho,v,u,P,T,cs,jr,rhops,ham_source,x,x_c,&gamma,&limiter,
		  &Ye,&dt,&rho0,&Nx,
		  psi_v_np1,psi,alpha_v_np1,alpha,beta_v_np1,beta,
		  ifc_mask,d_fcs,Sr_fcs,tau_fcs,&iter,&eos_flag,phys_bdy,
		  &deltat1,&deltat2,&rho_atm,&p_atm,&u_atm,&t_atm,&lorentz_max);

   }

   return;
}

//=============================================================================
// Calculations prior to saving info to disk.
//=============================================================================
void grhd_1d_pre_io_calc(void)
{
   return;
}

//=============================================================================
// Called after calculating the TRE for all variables
//=============================================================================
void grhd_1d_scale_tre(void)
{
   return;
}

//=============================================================================
// no post-regrid/tstep stuff
//=============================================================================
void grhd_1d_post_regrid(void)
{

  // We need to resolve for the primitives and calculate the metric source
  // terms!  We'll need them if we resolve the elliptics after regridding.
  ldptr();
  if (eos_flag==0) {
    find_primitives_(d_n,Sr_n,tau_n,rho,v,u,jr,rhops,ham_source,&gamma,alpha,beta,psi,x_c,
		     &Nx,&rho_atm,&lorentz_max);
  } else if (eos_flag==1) {
    find_primitives_table_(d_n,Sr_n,tau_n,rho,v,u,jr,rhops,ham_source,P,T,cs,&Ye,&rho0,
			   alpha,beta,psi,x_c,&Nx,&deltat1,&deltat2,&rho_atm,
			   &p_atm,&u_atm,&t_atm,&lorentz_max);
  }
}

// an empty post_tstep function.
void grhd_1d_post_tstep(int L)
{
  return;
}

// The pre-tstep function

void grhd_1d_pre_tstep(int L)
{

  int valid;
  static int local_first=1;
  static int step_counter = 1;
  char name[256];
  static FILE *out1;
  static FILE *out2;
  static FILE *out3;
  static FILE *out4;
  char out1_name[256];
  char out2_name[256];
  char out3_name[256];
  char out4_name[256];
  char buffer[20];
  int mpi_size;
  int Lf,Lc,write_lev,ict;
  int ind, i, j, k;
  real dt, dx0[3];
  real ham_norm;
  real m, lm,ct;
  real rho_c, lrho_c;
  real alpha_c, psi_c;
  int write_profile;

  Lf=PAMR_get_max_lev(PAMR_AMRH);
  Lc=PAMR_get_min_lev(PAMR_AMRH);
  
  if (Lf==1) write_lev = Lf;
  if (Lf>1)  write_lev = Lc + 1; // not the shadow level

  if (local_first && L==write_lev) {
    sprintf(out3_name,"grhd_1d_w.mon");
    if (my_rank==0) { if (!(out3=fopen(out3_name,"w"))) AMRD_stop("grhd_1d_pre_tstep: error opening out3 file\n",""); }
    local_first=0;
  }

  // Next, write some diagnostics to a file.
  if (L==write_lev) {
    PAMR_c_to_v(L,1,PAMR_AMRH,0);
    PAMR_c_to_v(L,2,PAMR_AMRH,0);

    if (my_rank == 0) {
      ct=PAMR_get_time(L);
      PAMR_get_dxdt(L,dx0,&dt);
      //      ict = ct * 1000 + 0.5;
      ict = ct/dt + 0.5;

      if (ict%ascii_dump_int == 0){
	sprintf(out1_name,"%s.hydro.dat",itoa(ict));
	//sprintf(out1_name,"out.dat");
	if (!(out1=fopen(out1_name,"w"))) {
	  AMRD_stop("grhd_1d_post_tstep: error opening output file\n","");
	}

	valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
      
	while(valid)
	  {
	    ldptr();
	    for (i=0; i<Nx; i++) {
	      if (eos_flag == 0) {
		fprintf(out1,"%15.6E %15.6E %15.6E %15.6E %15.6E\n",x_c[i],rho[i],v[i],(gamma-1.0)*u[i],d_n[i]);
	      } else {
		fprintf(out1,"%15.6E %15.6E %15.6E %15.6E %15.6E %15.6E\n",x_c[i],rho[i],v[i],P[i],d_n[i],T[i]);
	      }
	    }
	    valid=PAMR_next_g();
	  }
	
	fflush(out1);

	sprintf(out2_name,"%s.metric.dat",itoa(ict));
	if (!(out2=fopen(out2_name,"w"))) {
	  AMRD_stop("grhd_1d_post_tstep: error opening output file\n","");
	}
	
	valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
      
	while(valid)
	  {
	    ldptr();
	    for (i=0; i<Nx+1; i++) {
	      fprintf(out2,"%15.6E %15.6E %15.6E %15.6E\n",x[i],alpha_v_n[i],beta_v_n[i],psi_v_n[i]);
	    }
	    valid=PAMR_next_g();
	  }
	
	fflush(out2);
      }

      // Branson messing around.  Want to look at your refinement grid function.
/*       sprintf(out2_name,"tre.dat"); */
/*       if (!(out2=fopen(out2_name,"w"))) { */
/* 	AMRD_stop("grhd_1d_post_tstep: error opening tre.dat file\n",""); */
/*       } */
/*       valid=PAMR_init_s_iter(L,PAMR_AMRH,0); */
      
/*       while(valid) */
/* 	{ */
/* 	  ldptr(); */
/* 	  find_tre_hydro_(d_tre,rho,&Nx); */

/* 	  valid=PAMR_next_g(); */
	  
/* 	  for (i=0; i<Nx; i++) { */
/* 	    fprintf(out2,"%15.6E %15.6E\n",x_c[i],d_tre[i]); */
/* 	  } */
/* 	} */
      
/*       fflush(out2); */
    }
      
    valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
    m = 0.0;
    while(valid)
      {
	ldptr();
	lm = 0.0;
	total_mass_(&lm,d_n,x,&Nx,ghost_width,&my_rank);
	valid=PAMR_next_g();

	alpha_c = alpha_v_n[1];
	psi_c = psi_v_n[1];
	rho_c = rho[1];	
      }

    MPI_Allreduce(&lm,&m,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    // Calculate Hamiltonian constraint and write the norm.
    if (L==write_lev && my_rank == 0) {
      ct=PAMR_get_time(L);
      PAMR_get_dxdt(L,dx0,&dt);
      ict = ct/dt + 0.5;
      
      valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
      while(valid) { 
	ldptr();

	write_profile = 0;
	if (ict%ascii_dump_int == 0) {
	  sprintf(out4_name,"%s.ham.dat",itoa(ict));	  
	  write_profile = 1;
	}

	ham_const_(ham_v,alpha_v_n,beta_v_n,psi_v_n,tau_v,ham_source_v,mask,phys_bdy,x,
		   &ham_norm,&Nx,ghost_width,out4_name,&write_profile);

      	valid=PAMR_next_g();
      }
    }
    if (my_rank == 0) ct=PAMR_get_time(L);
    if (my_rank == 0) fprintf(out3,"%14.7E %22.15E %14.7E %14.7E %14.7E %14.7E\n ",
			      ct,m,ham_norm,rho_c,alpha_c,psi_c);
    fflush(out3);
    
  }

  step_counter++;
}

void grhd_1d_find_primitives()
{
  ldptr();
  if (eos_flag == 0) {
    find_primitives_(d_np1,Sr_np1,tau_np1,rho,v,u,jr,rhops,ham_source,&gamma,alpha,beta,psi,x_c,
		     &Nx,&rho_atm,&lorentz_max);
  } else if (eos_flag == 1) {
    find_primitives_table_(d_np1,Sr_np1,tau_np1,rho,v,u,jr,rhops,ham_source,P,T,cs,&Ye,&rho0,alpha,
			   beta,psi,x_c,&Nx,&deltat1,&deltat2,&rho_atm,&p_atm,&u_atm,&t_atm,
			   &lorentz_max);
  }
}

void grhd_1d_fill_ex_mask(double *mask, double *mask_c, int dim, int *shape, int *shape_c, 
			  double *bbox, double excised)
{
}

//=============================================================================
int main(int argc, char **argv)
{
  amrd_set_app_fcs_var_clear_hook(grhd_1d_fcs_var_clear);
  amrd_set_app_flux_correct_hook(grhd_1d_flux_correct);
  amrd_set_app_post_flux_correct_hook(grhd_1d_post_flux_correct);
  amrd_set_app_find_tre_hydro_hook(grhd_1d_find_tre_hydro);
  amrd_set_app_pre_tstep_hook(grhd_1d_pre_tstep);
  amrd_set_elliptic_vars_t0_init(grhd_1d_elliptic_vars_t0_init);
  amrd_set_app_find_primitives_hook(grhd_1d_find_primitives);
  amrd(argc,argv,&grhd_1d_id,&grhd_1d_var_pre_init,
       &grhd_1d_var_post_init, &grhd_1d_AMRH_var_clear,
       &grhd_1d_free_data, &grhd_1d_t0_cnst_data,
       &grhd_1d_evo_residual, &grhd_1d_MG_residual,
       &grhd_1d_evolve, &grhd_1d_MG_relax, &grhd_1d_L_op, 
       &grhd_1d_pre_io_calc,&grhd_1d_scale_tre,
       &grhd_1d_post_regrid,&grhd_1d_post_tstep,
       &grhd_1d_fill_ex_mask,0);
}

// A hook function to apply the flux correction.

void grhd_1d_flux_correct(void)
{
  ldptr();
  
  apply_fc_(d_n,Sr_n,tau_n,d_fcs,Sr_fcs,tau_fcs,&Nx);
  
  return;
}  

// A hook function to take care of any unfinished business
// following the flux correction step.

void grhd_1d_post_flux_correct(void)
{
  ldptr();
  if (eos_flag==0) {
    find_primitives_(d_n,Sr_n,tau_n,rho,v,u,jr,rhops,ham_source,&gamma,alpha,beta,psi,x_c,
		     &Nx,&rho_atm,&lorentz_max);
  } else if (eos_flag==1) {
    find_primitives_table_(d_n,Sr_n,tau_n,rho,v,u,jr,rhops,ham_source,P,T,cs,&Ye,&rho0,alpha,
			   beta,psi,x_c,&Nx,&deltat1,&deltat2,&rho_atm,&p_atm,&u_atm,&t_atm,
			   &lorentz_max);
  }
}


//----------------------------------------------------------------------
// A hook function to calculate the grid function(s) to be compared 
// with the refinement criterion.
//----------------------------------------------------------------------

void grhd_1d_find_tre_hydro(void)
{
  //  ldptr();
  //  find_tre_hydro_(d_tre,rho,&Nx);
  
}


//----------------------------------------------------------------------
// Hooks for the t0 resolving of the Hamiltonian constraint
//----------------------------------------------------------------------
void grhd_1d_elliptic_vars_t0_init(void)
{
   ldptr();

   int i;
   for (i=0; i<Nx+1; i++) {
     psi0_v[i] = psi_v[i];
   }

   return;
}  

void grhd_1d_t0_cnst_data(void)
{
   ldptr();

   int i;
   if (psi_v && const_ev==0) {
     for (i=0; i<Nx+1; i++) {
       psi_v[i] = psi0_v[i];
     }
   }
   
   return;
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

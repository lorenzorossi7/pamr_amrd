//=============================================================================
// application functions/variables for Newtonian Boson Stars
//=============================================================================

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pamr_w.h"
#include "amrd_w.h"
#include "num_w.h"
#include <mpi.h>

//=============================================================================
// some convenient, "local" global variables
//=============================================================================
double p_deplete;
double gamma;
int limiter;
double rho_atm;

double *P, *rho, *u, *vx, *vy, *vz;
double *d, *Sx, *Sy, *Sz, *E;
double *d_p, *Sx_p, *Sy_p, *Sz_p, *E_p;
double *d_fcs, *Sx_fcs, *Sy_fcs, *Sz_fcs, *E_fcs;
double *phi_v, *phi_lop_v, *phi_res_v, *phi_rhs_v;
double *rho_v, *rho_mg, *rho_v_mg;
double *phi_extrap_tm1_v, *phi_extrap_tm2_v;

double *d_n, *Sx_n, *Sy_n, *Sz_n, *E_n, *phi_n;
double *phi_v_n, *phi_v_np1;
double *d_np1, *Sx_np1, *Sy_np1, *Sz_np1, *E_np1, *phi_np1;

double *mask,*mask_mg,*fc_mask;

double *x,*y,*z, *x_c, *y_c, *z_c;
int shape[3],shape_c[3],ghost_width[6],Nx,Ny,Nz,phys_bdy[6],size,size_c;
double base_bbox[6],bbox[6],dx,dy,dz,dt,dto2;
int g_L;

int P_gfn, rho_gfn, u_gfn, vx_gfn, vy_gfn, vz_gfn;
int d_gfn, Sx_gfn, Sy_gfn, Sz_gfn, E_gfn;
int d_p_gfn, Sx_p_gfn, Sy_p_gfn, Sz_p_gfn, E_p_gfn;
int d_fcs_gfn, Sx_fcs_gfn, Sy_fcs_gfn, Sz_fcs_gfn, E_fcs_gfn;
int phi_v_gfn, phi_lop_v_gfn, phi_res_v_gfn, phi_rhs_v_gfn;
int rho_v_gfn, rho_mg_gfn, rho_v_mg_gfn;
int phi_extrap_tm1_v_gfn, phi_extrap_tm2_v_gfn;

int d_n_gfn, Sx_n_gfn, Sy_n_gfn, Sz_n_gfn, E_n_gfn, phi_v_n_gfn;
int d_np1_gfn, Sx_np1_gfn, Sy_np1_gfn, Sz_np1_gfn, E_np1_gfn, phi_v_np1_gfn;

int mask_gfn,mask_mg_gfn,fc_mask_gfn;

int max_lev;
int init_depth;
int dim;
double t0;
int ID_DV_trace;

void nps_free_data(void);
void nps_init_rest(void);
void nps_find_primitives(void);
void nps_fcs_var_clear(int type_flag, int *ifc_mask);
void nps_flux_correct(void);
void nps_post_flux_correct(void);
void nps_pre_tstep(int L);

//=============================================================================
// call after variables have been defined
//=============================================================================
void set_gfns(void)
{
    if ((d_gfn=PAMR_get_gfn("d",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((d_np1_gfn=PAMR_get_gfn("d",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((d_n_gfn=PAMR_get_gfn("d",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((d_p_gfn=PAMR_get_gfn("d_p",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((Sx_gfn=PAMR_get_gfn("Sx",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((Sx_np1_gfn=PAMR_get_gfn("Sx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Sx_n_gfn=PAMR_get_gfn("Sx",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Sx_p_gfn=PAMR_get_gfn("Sx_p",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((Sy_gfn=PAMR_get_gfn("Sy",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((Sy_np1_gfn=PAMR_get_gfn("Sy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Sy_n_gfn=PAMR_get_gfn("Sy",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Sy_p_gfn=PAMR_get_gfn("Sy_p",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((Sz_gfn=PAMR_get_gfn("Sz",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((Sz_np1_gfn=PAMR_get_gfn("Sz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Sz_n_gfn=PAMR_get_gfn("Sz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Sz_p_gfn=PAMR_get_gfn("Sz_p",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((E_gfn=PAMR_get_gfn("E",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((E_np1_gfn=PAMR_get_gfn("E",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((E_n_gfn=PAMR_get_gfn("E",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((E_p_gfn=PAMR_get_gfn("E_p",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((phi_v_gfn=PAMR_get_gfn("phi_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((phi_v_np1_gfn=PAMR_get_gfn("phi_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((phi_v_n_gfn=PAMR_get_gfn("phi_v",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((phi_extrap_tm1_v_gfn=PAMR_get_gfn("phi_extrap_tm1_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((phi_extrap_tm2_v_gfn=PAMR_get_gfn("phi_extrap_tm2_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((phi_res_v_gfn=PAMR_get_gfn("phi_res_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((phi_lop_v_gfn=PAMR_get_gfn("phi_lop_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((phi_rhs_v_gfn=PAMR_get_gfn("phi_rhs_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mask_mg_gfn=PAMR_get_gfn("cmask_v",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mask_gfn=PAMR_get_gfn("cmask_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    // Stuff for flux correction.  They are all work variables in the AMR hierarchy.
    if ((d_fcs_gfn   = PAMR_get_gfn("d_fcs"  , PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((Sx_fcs_gfn  = PAMR_get_gfn("Sx_fcs" , PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((Sy_fcs_gfn  = PAMR_get_gfn("Sy_fcs" , PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((Sz_fcs_gfn  = PAMR_get_gfn("Sz_fcs" , PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((E_fcs_gfn   = PAMR_get_gfn("E_fcs"  , PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 

    // Some work variables
    if ((rho_gfn=PAMR_get_gfn("rho",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((vx_gfn=PAMR_get_gfn("vx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((vy_gfn=PAMR_get_gfn("vy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((vz_gfn=PAMR_get_gfn("vz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((P_gfn=PAMR_get_gfn("P",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((u_gfn=PAMR_get_gfn("u",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((rho_mg_gfn=PAMR_get_gfn("rho",PAMR_MGH,0) )<0) AMRD_stop("set_gnfs error",0); 
    if ((rho_v_gfn=PAMR_get_gfn(   "rho_v",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0); 
    if ((rho_v_mg_gfn=PAMR_get_gfn("rho_v",PAMR_MGH,0) )<0) AMRD_stop("set_gnfs error",0);
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
   d_p=gfs[d_p_gfn-1];

   Sx=gfs[Sx_gfn-1];
   Sx_fcs=gfs[Sx_fcs_gfn-1];
   Sx_n=gfs[Sx_n_gfn-1];
   Sx_np1=gfs[Sx_np1_gfn-1];
   Sx_p=gfs[Sx_p_gfn-1];

   Sy=gfs[Sy_gfn-1];
   Sy_fcs=gfs[Sy_fcs_gfn-1];
   Sy_n=gfs[Sy_n_gfn-1];
   Sy_np1=gfs[Sy_np1_gfn-1];
   Sy_p=gfs[Sy_p_gfn-1];

   Sz=gfs[Sz_gfn-1];
   Sz_fcs=gfs[Sz_fcs_gfn-1];
   Sz_n=gfs[Sz_n_gfn-1];
   Sz_np1=gfs[Sz_np1_gfn-1];
   Sz_p=gfs[Sz_p_gfn-1];

   E=gfs[E_gfn-1];
   E_fcs=gfs[E_fcs_gfn-1];
   E_n=gfs[E_n_gfn-1];
   E_np1=gfs[E_np1_gfn-1];
   E_p=gfs[E_p_gfn-1];

   phi_v=gfs[phi_v_gfn-1];
   phi_v_n=gfs[phi_v_n_gfn-1];
   phi_v_np1=gfs[phi_v_np1_gfn-1];

   mask=gfs[mask_gfn-1];
   mask_mg=gfs[mask_mg_gfn-1];
   phi_res_v=gfs[phi_res_v_gfn-1];
   phi_rhs_v=gfs[phi_rhs_v_gfn-1];
   phi_lop_v=gfs[phi_lop_v_gfn-1];

   phi_extrap_tm1_v=gfs[phi_extrap_tm1_v_gfn-1];
   phi_extrap_tm2_v=gfs[phi_extrap_tm2_v_gfn-1];

   rho   = gfs[rho_gfn-1];
   rho_v = gfs[rho_v_gfn-1];
   rho_mg = gfs[rho_mg_gfn-1];
   rho_v_mg = gfs[rho_v_mg_gfn-1];
   vx    = gfs[vx_gfn-1];
   vy    = gfs[vy_gfn-1];
   vz    = gfs[vz_gfn-1];
   P     = gfs[P_gfn-1];
   u     = gfs[u_gfn-1];

   x=x0[0]; dx=x[1]-x[0];
   y=x0[1]; dy=y[1]-y[0];
   z=x0[2]; dz=z[1]-z[0];
   x_c=x0_c[0]; 
   y_c=x0_c[1]; 
   z_c=x0_c[2]; 
   PAMR_get_g_level(&lev);
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

   size_c=Nx*Ny*Nz;
}

//=============================================================================
// utility routines
//=============================================================================
void zero(double *f)
{
   int i;
   for (i=0; i<shape[0]*shape[1]*shape[2]; i++) f[i]=0;
}

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
int nps_id(void)
{
  int valid;
  int Lf,L,Lc;
  int gnum;
  int glev[PAMR_MAX_LEVS];
  real gbbox[PAMR_MAX_LEVS*PAMR_MAX_DIM*2];
  int i,j;
  real frac = 0.5;
  // Warning:  This assumes that, if you want more than one level, 
  // you'll be specifying them yourself, a la FMR.

  if (max_lev > 1) {
    // Devise your boxes (for levels higher than the base level)
    gnum = max_lev - init_depth;
    if (gnum < 1) printf("problem inside nps_id\n");

    for (i=1; i<=gnum; i++) {
      glev[i-1] = init_depth + i;
      for (j=0; j<2*dim; j++) gbbox[(i-1)*2*dim + j] = frac*base_bbox[j];
      frac = frac * frac;
    }

    // Compose the hierarchy
    t0 = 0.0;

    if (!(PAMR_compose_hierarchy(init_depth+1,max_lev,gnum,glev,gbbox,t0))) printf("ouch\n");
    // Fill grid functions on all levels with the initial data.
    Lf=PAMR_get_max_lev(PAMR_AMRH);
    Lc=PAMR_get_min_lev(PAMR_AMRH);

    for (L=1; L<=Lf; L++) {
      valid = PAMR_init_s_iter(L,PAMR_AMRH,0);
      while(valid) {
	nps_free_data();
	valid = PAMR_next_g();
      }
      PAMR_sync(L,2,PAMR_AMRH,0);
    }

    for (L=Lf; L>1; L--) {
      PAMR_inject(L,2,PAMR_AMRH);
      valid = PAMR_init_s_iter(L-1,PAMR_AMRH,0);
      while(valid) {
	nps_find_primitives();
	valid = PAMR_next_g();
      }
    }

    // Dump the initial data stuff if desired.
    if (ID_DV_trace) {
      PAMR_save_gfn("d",PAMR_AMRH,2,-1,-1.0,"ID_","");
      PAMR_save_gfn("Sx",PAMR_AMRH,2,-1,-1.0,"ID_","");
      PAMR_save_gfn("Sy",PAMR_AMRH,2,-1,-1.0,"ID_","");
      PAMR_save_gfn("Sz",PAMR_AMRH,2,-1,-1.0,"ID_","");
      PAMR_save_gfn("E",PAMR_AMRH,2,-1,-1.0,"ID_","");
      PAMR_save_gfn("phi_v",PAMR_AMRH,2,-1,-1.0,"ID_","");
    }

    // You need to do the stuff contained in init_rest
    for (L=1; L<=Lf; L++) {
      valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
      while(valid) {
	nps_init_rest();
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
void nps_var_pre_init(char *pfile)
{
  // Read in some global parameters which are needed by nps_id
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

void nps_var_post_init(char *pfile)
{

   if (my_rank==0)
   {
      printf("===================================================================\n");
      printf("Reading NPS parameters:\n\n");
   }

   p_deplete = 1.0;
   gamma = 2.0;
   limiter = 1;
   rho_atm = 1.e-5;

   AMRD_real_param(pfile,"p_deplete",&p_deplete,1);
   AMRD_real_param(pfile,"gamma",&gamma,1); 
   AMRD_int_param(pfile,"limiter",&limiter,1); 
   AMRD_real_param(pfile,"rho_atm",&rho_atm,1);

   if (my_rank==0) printf("===================================================================\n");
   return;
}

//=============================================================================
// Sets all variables to their 'zero' values:
//=============================================================================
void nps_AMRH_var_clear(void)
{
   ldptr();

   zero_c(d_n); zero_c(d_np1);
   zero_c(Sx_n); zero_c(Sx_np1);
   zero_c(Sy_n); zero_c(Sy_np1);
   zero_c(Sz_n); zero_c(Sz_np1);
   zero(phi_v_n); zero(phi_v_np1);
   
   return;
}

//=============================================================================
// Sets fcs variables to their 'zero' values.
// typeflag = 0:  corresponds to negative sign flag in mask.  zero cells of type A
// typeflag = 1:  corresponds to positive sign flag in mask.  zero cells of type B
//=============================================================================
void nps_fcs_var_clear(int type_flag, int *ifc_mask)
{
   ldptr();
   zero_fcs_vars_(&type_flag,ifc_mask,d_fcs,Sx_fcs,Sy_fcs,Sz_fcs,E_fcs,&Nx,&Ny,&Nz);
   return;
}

//=============================================================================
// Initial data for Newtonian polytrope: (at tn=2)
//=============================================================================
void nps_free_data(void)
{
   int i;

   ldptr();
   
   polytrope_id_(&gamma,d_n,Sx_n,Sy_n,Sz_n,E_n,rho,rho_v,vx,vy,vz,
		 u,phi_v_n,&p_deplete,x,y,z,x_c,y_c,z_c,&Nx,&Ny,&Nz,&rho_atm);

   return;
}  

void nps_init_rest(void)
{
   ldptr();
   
   init_rest_(d_n,d_np1,Sx_n,Sx_np1,Sy_n,Sy_np1,Sz_n,Sz_np1,E_n,E_np1,
	      phi_v_n,phi_v_np1,phi_extrap_tm1_v,phi_extrap_tm2_v,
	      &Nx,&Ny,&Nz);

   return;
}  

//=============================================================================
// Initial constraint data --- called after each MG iteration:
//=============================================================================
void nps_t0_cnst_data(void)
{
   return;
}

//=============================================================================
// Returns some norm of the residual for the evolution variables,
// called after an evolution iteration:
//=============================================================================
double nps_evo_residual(void)
{
   ldptr();

   return 0.0;
}

//=============================================================================
// Returns some norm of the residual for the MG variables, *AND* 
// stores the point-wise residual in "f_res" for each MG variable "f" (for
// computing new RHS's)
//=============================================================================
double nps_MG_residual(void)
{
   double norm;

   ldptr();

   residual_(phi_res_v,phi_rhs_v,phi_v,rho_v_mg,mask_mg,x,y,z,&norm,&Nx,&Ny,&Nz);
   return norm;
}


//=============================================================================
// Performs 1 iteration of the evolution equations 
//=============================================================================
void nps_evolve(int iter, int *ifc_mask)
{
   ldptr();
   int j;

   if (iter == 1) {
     // Perform half step.  
     hydro_1step_(&iter,d_np1,Sx_np1,Sy_np1,Sz_np1,E_np1,phi_v_np1,
		  d_p,Sx_p,Sy_p,Sz_p,E_p,
		  rho,vx,vy,vz,u,x_c,y_c,z_c,&gamma,&limiter,&rho_atm,
		  &dto2,&Nx,&Ny,&Nz,&p_deplete,ifc_mask,d_fcs,Sx_fcs,Sy_fcs,
		  Sz_fcs,E_fcs);

     for (j=0; j<size_c; j++) {
       d_np1[j]  = d_n[j]  + dto2*d_p[j];
       Sx_np1[j] = Sx_n[j] + dto2*Sx_p[j];
       Sy_np1[j] = Sy_n[j] + dto2*Sy_p[j];
       Sz_np1[j] = Sz_n[j] + dto2*Sz_p[j];
       E_np1[j]  = E_n[j]  + dto2*E_p[j];
     }

   } else if (iter==2) {
     // Perform full step.
     hydro_1step_(&iter,d_np1,Sx_np1,Sy_np1,Sz_np1,E_np1,phi_v_np1,
		  d_p,Sx_p,Sy_p,Sz_p,E_p,
		  rho,vx,vy,vz,u,x_c,y_c,z_c,&gamma,&limiter,&rho_atm,
		  &dt,&Nx,&Ny,&Nz,&p_deplete,ifc_mask,d_fcs,Sx_fcs,Sy_fcs,
		  Sz_fcs,E_fcs);

     for (j=0; j<size_c; j++) {
       d_np1[j]  = d_n[j]  + dt*d_p[j];
       Sx_np1[j] = Sx_n[j] + dt*Sx_p[j];
       Sy_np1[j] = Sy_n[j] + dt*Sy_p[j];
       Sz_np1[j] = Sz_n[j] + dt*Sz_p[j];
       E_np1[j]  = E_n[j]  + dt*E_p[j];
     }

   } else if (iter==3) {
     // Final primitive variable inversion.
     hydro_1step_(&iter,d_np1,Sx_np1,Sy_np1,Sz_np1,E_np1,phi_v_np1,
		  d_p,Sx_p,Sy_p,Sz_p,E_p,
		  rho,vx,vy,vz,u,x_c,y_c,z_c,&gamma,&limiter,&rho_atm,
		  &dt,&Nx,&Ny,&Nz,&p_deplete,ifc_mask,d_fcs,Sx_fcs,Sy_fcs,
		  Sz_fcs,E_fcs);
   }
   
   return;
}

//=============================================================================
// Performs 1 relaxation sweep of the MG equations, and returns an estimate
// of the norm of the residual.
//=============================================================================
double nps_MG_relax(void)
{
   double norm;
   ldptr();

   relax_(phi_v,phi_rhs_v,rho_v_mg,mask_mg,phys_bdy,x,y,z,&norm,&Nx,&Ny,&Nz,ghost_width);

   return norm;
}

//=============================================================================
// Computes the differential operator L acting on the current grid,
// in a region specified by the grid function "mask". Stores the result
// in "f_lop" for each MG variable "f"
//=============================================================================
void nps_L_op(void)
{
   ldptr();

   lop_(phi_lop_v,phi_v,rho_v_mg,mask_mg,x,y,z,&Nx,&Ny,&Nz);

   return;
}

//=============================================================================
// Calculations prior to saving info to disk.
//=============================================================================
void nps_pre_io_calc(void)
{
   return;
}

//=============================================================================
// Called after calculating the TRE for all variables
//=============================================================================
void nps_scale_tre(void)
{
   return;
}

//=============================================================================
// no post-regrid/tstep stuff
//=============================================================================
void nps_post_regrid(void)
{
}
   
void nps_post_tstep(int L)
{
   int valid; 
  valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
  while(valid)
    {
      ldptr();
      find_primitives_(d_n,Sx_n,Sy_n,Sz_n,E_n,phi_v_n,rho,vx,vy,vz,u,&gamma,&Nx,&Ny,&Nz,&rho_atm,&p_deplete);
      valid=PAMR_next_g();
    }

}

void nps_pre_tstep(int L)
{
  int valid;
  static int local_first=1;
  char name[256];
  static FILE *out1;
  char out1_name[256];
  int mpi_size;
  int Lf,Lc,write_lev;
  int has_origin,lhas_origin;
  real m, lm,ct;
  real m_ana;
  real lrho_c, rho_c;

  // This is where you need to resolve for the primitive variables.
  // But you must realize that the time levels have already been 
  // cycled.
  valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
  while(valid)
    {
      ldptr();
      find_primitives_(d_n,Sx_n,Sy_n,Sz_n,E_n,phi_v_n,rho,vx,vy,vz,u,&gamma,&Nx,&Ny,&Nz,&rho_atm,&p_deplete);
      valid=PAMR_next_g();
    }

  Lf=PAMR_get_max_lev(PAMR_AMRH);
  Lc=PAMR_get_min_lev(PAMR_AMRH);
  
  if (Lf==1) write_lev = Lf;
  if (Lf>1)  write_lev = Lc + 1; // not the shadow level

  // Next, write some diagnostics to a file.
  if (local_first && L==write_lev)
    {
      sprintf(out1_name,"%s.out1",AMRD_save_tag);
      if (my_rank==0) { if (!(out1=fopen(out1_name,"w"))) AMRD_stop("gh3d_post_tstep: error opening out1 file\n",""); }
      local_first=0;
    }

  if (L==write_lev) {
    valid=PAMR_init_s_iter(L,PAMR_AMRH,0);

    while(valid)
      {
	ldptr();
	total_mass_(&lm,rho,x,y,z,&Nx,&Ny,&Nz,ghost_width,&my_rank);
	get_rho_c_(&lhas_origin,&lrho_c,rho_v,x,y,z,&Nx,&Ny,&Nz,ghost_width);
	valid=PAMR_next_g();
      }
    
    rho_c = 0;
    m = 0;
    if (my_rank == 0) ct=PAMR_get_time(L);
    MPI_Allreduce(&lrho_c,&rho_c,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&lhas_origin,&has_origin,3,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (has_origin) rho_c = rho_c / has_origin;
    MPI_Allreduce(&lm,&m,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if (my_rank == 0) fprintf(out1,"%15.6E %15.6E %20.11E\n",ct,rho_c,m);
    fflush(out1);
  }

}

void nps_find_primitives()
{
  ldptr();
  find_primitives_(d_n,Sx_n,Sy_n,Sz_n,E_n,phi_v_n,rho,vx,vy,vz,u,&gamma,&Nx,&Ny,&Nz,&rho_atm,&p_deplete);
}

void nps_fill_ex_mask(double *mask, double *mask_c, int dim, int *shape, int *shape_c, double *bbox, double excised)
{
}

//=============================================================================
int main(int argc, char **argv)
{
  amrd_set_app_fcs_var_clear_hook(nps_fcs_var_clear);
  amrd_set_app_flux_correct_hook(nps_flux_correct);
  amrd_set_app_post_flux_correct_hook(nps_post_flux_correct);
  amrd_set_app_pre_tstep_hook(nps_pre_tstep);
  amrd(argc,argv,&nps_id,&nps_var_pre_init,
       &nps_var_post_init, &nps_AMRH_var_clear,
       &nps_free_data, &nps_t0_cnst_data,
       &nps_evo_residual, &nps_MG_residual,
       &nps_evolve, &nps_MG_relax, &nps_L_op, 
       &nps_pre_io_calc,&nps_scale_tre,
       &nps_post_regrid,&nps_post_tstep,
       &nps_fill_ex_mask,0);
}

// A hook function to apply the flux correction.

void nps_flux_correct(void)
{
  ldptr();

  apply_fc_(d_n,Sx_n,Sy_n,Sz_n,E_n,d_fcs,Sx_fcs,Sy_fcs,Sz_fcs,E_fcs,&Nx,&Ny,&Nz);
  
  return;
}  

// A hook function to take care of any unfinished business
// following the flux correction step.

void nps_post_flux_correct(void)
{
  ldptr();
  find_primitives_(d_n,Sx_n,Sy_n,Sz_n,E_n,phi_v_n,rho,vx,vy,vz,u,&gamma,&Nx,&Ny,&Nz,&rho_atm,&p_deplete);
  
}

//=============================================================================
// application functions/variables for Newtonian Boson Stars
//=============================================================================

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pamr_w.h"
#include "amrd_w.h"
#include "num_w.h"

//=============================================================================
// id parameters
//=============================================================================

real phi_r_amp_1,phi_r_r0_1,phi_r_delta_1,phi_r_x0_1[3],phi_r_ecc_1[3];
real phi_i_amp_1,phi_i_r0_1,phi_i_delta_1,phi_i_x0_1[3],phi_i_ecc_1[3];
real phi_r_amp_2,phi_r_r0_2,phi_r_delta_2,phi_r_x0_2[3],phi_r_ecc_2[3];
real phi_i_amp_2,phi_i_r0_2,phi_i_delta_2,phi_i_x0_2[3],phi_i_ecc_2[3];
real phi_r_amp_3,phi_r_r0_3,phi_r_delta_3,phi_r_x0_3[3],phi_r_ecc_3[3];
real phi_i_amp_3,phi_i_r0_3,phi_i_delta_3,phi_i_x0_3[3],phi_i_ecc_3[3];

//=============================================================================
// some convenient, "local" global variables
//=============================================================================
real *phi_r,*phi_i,*V,*V_lop,*V_res,*V_rhs;
real *phi_r_n,*phi_i_n,*V_n;
real *phi_r_np1,*phi_i_np1,*V_np1;
real *mask,*mask_mg,*w1,*se_res;

real *x,*y,*z;
int shape[3],ghost_width[6],Nx,Ny,Nz,phys_bdy[6],size;
real base_bbox[6],bbox[6],dx,dy,dz,dt;
int g_L;

int phi_r_gfn,phi_i_gfn,V_gfn,V_lop_gfn,V_res_gfn,V_rhs_gfn,mask_gfn,mask_mg_gfn;
int phi_r_n_gfn,phi_i_n_gfn,V_n_gfn;
int phi_r_np1_gfn,phi_i_np1_gfn,V_np1_gfn;
int w1_gfn,se_res_gfn;

//=============================================================================
// call after variables have been defined
//=============================================================================
void set_gfns(void)
{
    if ((phi_r_gfn=PAMR_get_gfn("phi_r",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((phi_r_np1_gfn=PAMR_get_gfn("phi_r",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((phi_r_n_gfn=PAMR_get_gfn("phi_r",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((phi_i_gfn=PAMR_get_gfn("phi_i",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((phi_i_np1_gfn=PAMR_get_gfn("phi_i",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((phi_i_n_gfn=PAMR_get_gfn("phi_i",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((V_gfn=PAMR_get_gfn("V",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((V_np1_gfn=PAMR_get_gfn("V",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((V_n_gfn=PAMR_get_gfn("V",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((V_res_gfn=PAMR_get_gfn("V_res",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((V_lop_gfn=PAMR_get_gfn("V_lop",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((V_rhs_gfn=PAMR_get_gfn("V_rhs",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mask_mg_gfn=PAMR_get_gfn("cmask",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mask_gfn=PAMR_get_gfn("cmask",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((w1_gfn=PAMR_get_gfn("w1",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((se_res_gfn=PAMR_get_gfn("se_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
}

//=============================================================================
// call with valid iter to set up globals:
//=============================================================================
void ldptr(void)
{
   int rank,dim,ngfs,lev,shape_c[3];
   real t,*x0[3],*gfs[100],dx0[3],*x0_c[3];
   static int first=1;

   if (first) 
   {
      first=0; 
      set_gfns();
      PAMR_get_global_bbox(base_bbox);
   }

   if (!(PAMR_get_g_attribs(&rank,&dim,shape,shape_c,bbox,ghost_width,&t,&ngfs,x0,x0_c,gfs))) 
      AMRD_stop("ldptr: PAMR_get_g_attribs failed\n","");

   phi_r=gfs[phi_r_gfn-1];
   phi_r_n=gfs[phi_r_n_gfn-1];
   phi_r_np1=gfs[phi_r_np1_gfn-1];

   phi_i=gfs[phi_i_gfn-1];
   phi_i_n=gfs[phi_i_n_gfn-1];
   phi_i_np1=gfs[phi_i_np1_gfn-1];

   V=gfs[V_gfn-1];
   V_n=gfs[V_n_gfn-1];
   V_np1=gfs[V_np1_gfn-1];

   mask=gfs[mask_gfn-1];
   mask_mg=gfs[mask_mg_gfn-1];
   V_res=gfs[V_res_gfn-1];
   V_rhs=gfs[V_rhs_gfn-1];
   V_lop=gfs[V_lop_gfn-1];

   w1=gfs[w1_gfn-1];
   se_res=gfs[se_res_gfn-1];

   x=x0[0]; dx=x[1]-x[0];
   y=x0[1]; dy=y[1]-y[0];
   z=x0[2]; dz=z[1]-z[0];
   PAMR_get_g_level(&lev);
   PAMR_get_dxdt(lev,dx0,&dt);

   if ((bbox[0]-base_bbox[0])<dx/2) phys_bdy[0]=1; else phys_bdy[0]=0;
   if ((base_bbox[1]-bbox[1])<dx/2) phys_bdy[1]=1; else phys_bdy[1]=0;
   if ((bbox[2]-base_bbox[2])<dy/2) phys_bdy[2]=1; else phys_bdy[2]=0;
   if ((base_bbox[3]-bbox[3])<dy/2) phys_bdy[3]=1; else phys_bdy[3]=0;
   if ((bbox[4]-base_bbox[4])<dz/2) phys_bdy[4]=1; else phys_bdy[4]=0;
   if ((base_bbox[5]-bbox[5])<dz/2) phys_bdy[5]=1; else phys_bdy[5]=0;

   Nx=shape[0];
   Ny=shape[1];
   Nz=shape[2];

   size=Nx*Ny*Nz;
}

//=============================================================================
// utility routines
//=============================================================================
void zero(real *f)
{
   int i;

   for (i=0; i<shape[0]*shape[1]*shape[2]; i++) f[i]=0;
}

//=============================================================================
// Routines required by amrd:
//=============================================================================

//=============================================================================
// Returns 0 to use default mechanism, or is expected to calculate
// the correct initial hierarchy and return 1:
//=============================================================================
int nbs_id(void)
{
   return 0;
}

//=============================================================================
// Sets custom parameters, variables, etc. Split up into two segments,
// one called before the pamr context is initialized and standard
// parameters are read, and the other afterwards
//=============================================================================
void nbs_var_pre_init(char *pfile)
{
   return;
}

void nbs_var_post_init(char *pfile)
{

   if (my_rank==0)
   {
      printf("===================================================================\n");
      printf("Reading NBS parameters:\n\n");
   }

   phi_r_amp_1=phi_r_r0_1=phi_r_x0_1[0]=phi_r_x0_1[1]=phi_r_x0_1[2]=phi_r_ecc_1[0]=phi_r_ecc_1[1]=phi_r_ecc_1[2]=0;
   phi_i_amp_1=phi_i_r0_1=phi_i_x0_1[0]=phi_i_x0_1[1]=phi_i_x0_1[2]=phi_i_ecc_1[0]=phi_i_ecc_1[1]=phi_i_ecc_1[2]=0;
   phi_i_delta_1=phi_r_delta_1=1;

   phi_r_amp_2=phi_r_r0_2=phi_r_x0_2[0]=phi_r_x0_2[1]=phi_r_x0_2[2]=phi_r_ecc_2[0]=phi_r_ecc_2[1]=phi_r_ecc_2[2]=0;
   phi_i_amp_2=phi_i_r0_2=phi_i_x0_2[0]=phi_i_x0_2[1]=phi_i_x0_2[2]=phi_i_ecc_2[0]=phi_i_ecc_2[1]=phi_i_ecc_2[2]=0;
   phi_i_delta_2=phi_r_delta_2=1;

   phi_r_amp_3=phi_r_r0_3=phi_r_x0_3[0]=phi_r_x0_3[1]=phi_r_x0_3[2]=phi_r_ecc_3[0]=phi_r_ecc_3[1]=phi_r_ecc_3[2]=0;
   phi_i_amp_3=phi_i_r0_3=phi_i_x0_3[0]=phi_i_x0_3[1]=phi_i_x0_3[2]=phi_i_ecc_3[0]=phi_i_ecc_3[1]=phi_i_ecc_3[2]=0;
   phi_i_delta_3=phi_r_delta_3=1;

   AMRD_real_param(pfile,"phi_r_amp_1",&phi_r_amp_1,1);
   AMRD_real_param(pfile,"phi_r_r0_1",&phi_r_r0_1,1);
   AMRD_real_param(pfile,"phi_r_delta_1",&phi_r_delta_1,1);
   AMRD_real_param(pfile,"phi_r_x0_1",phi_r_x0_1,3);
   AMRD_real_param(pfile,"phi_r_ecc_1",phi_r_ecc_1,3);

   AMRD_real_param(pfile,"phi_i_amp_1",&phi_i_amp_1,1);
   AMRD_real_param(pfile,"phi_i_r0_1",&phi_i_r0_1,1);
   AMRD_real_param(pfile,"phi_i_delta_1",&phi_i_delta_1,1);
   AMRD_real_param(pfile,"phi_i_x0_1",phi_i_x0_1,3);
   AMRD_real_param(pfile,"phi_i_ecc_1",phi_i_ecc_1,3);

   AMRD_real_param(pfile,"phi_r_amp_2",&phi_r_amp_2,1);
   AMRD_real_param(pfile,"phi_r_r0_2",&phi_r_r0_2,1);
   AMRD_real_param(pfile,"phi_r_delta_2",&phi_r_delta_2,1);
   AMRD_real_param(pfile,"phi_r_x0_2",phi_r_x0_2,3);
   AMRD_real_param(pfile,"phi_r_ecc_2",phi_r_ecc_2,3);

   AMRD_real_param(pfile,"phi_i_amp_2",&phi_i_amp_2,1);
   AMRD_real_param(pfile,"phi_i_r0_2",&phi_i_r0_2,1);
   AMRD_real_param(pfile,"phi_i_delta_2",&phi_i_delta_2,1);
   AMRD_real_param(pfile,"phi_i_x0_2",phi_i_x0_2,3);
   AMRD_real_param(pfile,"phi_i_ecc_2",phi_i_ecc_2,3);

   AMRD_real_param(pfile,"phi_r_amp_3",&phi_r_amp_3,1);
   AMRD_real_param(pfile,"phi_r_r0_3",&phi_r_r0_3,1);
   AMRD_real_param(pfile,"phi_r_delta_3",&phi_r_delta_3,1);
   AMRD_real_param(pfile,"phi_r_x0_3",phi_r_x0_3,3);
   AMRD_real_param(pfile,"phi_r_ecc_3",phi_r_ecc_3,3);

   AMRD_real_param(pfile,"phi_i_amp_3",&phi_i_amp_3,1);
   AMRD_real_param(pfile,"phi_i_r0_3",&phi_i_r0_3,1);
   AMRD_real_param(pfile,"phi_i_delta_3",&phi_i_delta_3,1);
   AMRD_real_param(pfile,"phi_i_x0_3",phi_i_x0_3,3);
   AMRD_real_param(pfile,"phi_i_ecc_3",phi_i_ecc_3,3);

   if (my_rank==0) printf("===================================================================\n");
   return;
}

//=============================================================================
// Sets all variables to their 'zero' values:
//=============================================================================
void nbs_AMRH_var_clear(void)
{
   ldptr();

   zero(phi_r_n); zero(phi_r_np1);
   zero(phi_i_n); zero(phi_i_np1);
   zero(V_n); zero(V_np1);
   
   return;
}

//=============================================================================
// Initial data for free fields: (at tn=2)
//=============================================================================
void nbs_free_data(void)
{
   int i;

   ldptr();

   gauss3d_(phi_r_n,&phi_r_amp_1,&phi_r_r0_1,&phi_r_delta_1,&phi_r_x0_1[0],&phi_r_x0_1[1],
            &phi_r_x0_1[2],&phi_r_ecc_1[0],&phi_r_ecc_1[1],&phi_r_ecc_1[2],x,y,z,&Nx,&Ny,&Nz);

   gauss3d_(w1,&phi_r_amp_2,&phi_r_r0_2,&phi_r_delta_2,&phi_r_x0_2[0],&phi_r_x0_2[1],
            &phi_r_x0_2[2],&phi_r_ecc_2[0],&phi_r_ecc_2[1],&phi_r_ecc_2[2],x,y,z,&Nx,&Ny,&Nz);
   for (i=0; i<size; i++) phi_r_n[i]+=w1[i];

   gauss3d_(w1,&phi_r_amp_3,&phi_r_r0_3,&phi_r_delta_3,&phi_r_x0_3[0],&phi_r_x0_3[1],
            &phi_r_x0_3[2],&phi_r_ecc_3[0],&phi_r_ecc_3[1],&phi_r_ecc_3[2],x,y,z,&Nx,&Ny,&Nz);
   for (i=0; i<size; i++) phi_r_n[i]+=w1[i];

   gauss3d_(phi_i_n,&phi_i_amp_1,&phi_i_r0_1,&phi_i_delta_1,&phi_i_x0_1[0],&phi_i_x0_1[1],
            &phi_i_x0_1[2],&phi_i_ecc_1[0],&phi_i_ecc_1[1],&phi_i_ecc_1[2],x,y,z,&Nx,&Ny,&Nz);

   gauss3d_(w1,&phi_i_amp_2,&phi_i_r0_2,&phi_i_delta_2,&phi_i_x0_2[0],&phi_i_x0_2[1],
            &phi_i_x0_2[2],&phi_i_ecc_2[0],&phi_i_ecc_2[1],&phi_i_ecc_2[2],x,y,z,&Nx,&Ny,&Nz);
   for (i=0; i<size; i++) phi_i_n[i]+=w1[i];

   gauss3d_(w1,&phi_i_amp_3,&phi_i_r0_3,&phi_i_delta_3,&phi_i_x0_3[0],&phi_i_x0_3[1],
            &phi_i_x0_3[2],&phi_i_ecc_3[0],&phi_i_ecc_3[1],&phi_i_ecc_3[2],x,y,z,&Nx,&Ny,&Nz);
   for (i=0; i<size; i++) phi_i_n[i]+=w1[i];

   return;
}  

//=============================================================================
// Initial constraint data --- called after each MG iteration:
//=============================================================================
void nbs_t0_cnst_data(void)
{
   return;
}

//=============================================================================
// Returns some norm of the residual for the evolution variables,
// called after an evolution iteration:
//=============================================================================
real nbs_evo_residual(void)
{
   real l2norm;

   ldptr();

   se_ires_(se_res,&l2norm,phi_r_n,phi_r_np1,phi_i_n,phi_i_np1,V_n,V_np1,AMRD_cmask, 
             x,y,z,&dt,&Nx,&Ny,&Nz);

   return l2norm;
}

//=============================================================================
// Returns some norm of the residual for the MG variables, *AND* 
// stores the point-wise residual in "f_res" for each MG variable "f" (for
// computing new RHS's)
//=============================================================================
real nbs_MG_residual(void)
{
   real norm;

   ldptr();

   residual_(V_res,V_rhs,V,phi_r,phi_i,mask_mg,x,y,z,&norm,&Nx,&Ny,&Nz);

   return norm;
}

//=============================================================================
// Performs 1 iteration of the evolution equations 
//=============================================================================
void nbs_evolve(int iter, int *ifc_mask)
{
   ldptr();

   phi_1step_cnc_(phi_r_n,phi_r_np1,phi_i_n,phi_i_np1,V_n,V_np1,
                  AMRD_cmask,x,y,z,&dt,&Nx,&Ny,&Nz);
   return;
}

//=============================================================================
// Performs 1 relaxation sweep of the MG equations, and returns an estimate
// of the norm of the residual.
//=============================================================================
real nbs_MG_relax(void)
{
   real norm;

   ldptr();

   relax_(V,V_rhs,phi_r,phi_i,mask_mg,phys_bdy,x,y,z,&norm,&Nx,&Ny,&Nz);

   return norm;
}

//=============================================================================
// Computes the differential operator L acting on the current grid,
// in a region specified by the grid function "mask". Stores the result
// in "f_lop" for each MG variable "f"
//=============================================================================
void nbs_L_op(void)
{
   ldptr();

   lop_(V_lop,V,phi_r,phi_i,mask_mg,x,y,z,&Nx,&Ny,&Nz);

   return;
}

//=============================================================================
// Calculations prior to saving info to disk.
//=============================================================================
void nbs_pre_io_calc(void)
{
   return;
}

//=============================================================================
// Called after calculating the TRE for all variables
//=============================================================================
void nbs_scale_tre(void)
{
   return;
}

//=============================================================================
// no post-regrid/tstep stuff
//=============================================================================
void nbs_post_regrid(void)
{
}
   
void nbs_post_tstep(int L)
{
}

void nbs_fill_ex_mask(real *mask, real *mask_c, int dim, int *shape, int *shape_c, real *bbox, real excised)
{
}

//=============================================================================
int main(int argc, char **argv)
{
   amrd(argc,argv,&nbs_id,&nbs_var_pre_init,
        &nbs_var_post_init, &nbs_AMRH_var_clear,
        &nbs_free_data, &nbs_t0_cnst_data,
        &nbs_evo_residual, &nbs_MG_residual,
        &nbs_evolve, &nbs_MG_relax, &nbs_L_op, 
        &nbs_pre_io_calc,&nbs_scale_tre,
        &nbs_post_regrid,&nbs_post_tstep,
        &nbs_fill_ex_mask,0);
}


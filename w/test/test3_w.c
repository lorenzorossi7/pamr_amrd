//=============================================================================
//
// test3.c ---- A couple of simple tests of the PAMR infrastructure
//              for CC support.
//
// Test types:
// 1  sync
// 2  c_to_v
// 3  v_to_c
// 4  interpolation
// 5  boundary interpolation
// 6  injection
// 7  regridding
// 8  multigrid hierarchy
// 9  excision
//
//=============================================================================

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include "pamr_w.h"
#include "../src/misc_w.h"
#include <math.h>


#define MAX_GRIDS 100
#define EXCISED 1
#define NOT_EXCISED 0
//------------------------------------------------------------------------------
// for excision ... excised a random sphere from the computational domain
//------------------------------------------------------------------------------
real ex_R,ex_x,ex_y,ex_z;
void set_ex_mask(real *mask, real *mask_c, int dim, int *shape, int *shape_c, real *bbox, real excised)
{
   int i,j,k,ie,je=0,ke=0;
   real x,y,z;
   real dx,dy,dz;

   ie=shape[0]-1; dx=(bbox[1]-bbox[0])/(shape[0]-1);
   if (dim>1) { je=shape[1]-1; dy=(bbox[3]-bbox[2])/(shape[1]-1); }
   if (dim>2) { ke=shape[2]-1; dz=(bbox[5]-bbox[4])/(shape[2]-1); }

   for (i=0; i<=ie; i++)
      for (j=0; j<=je; j++)
         for (k=0; k<=ke; k++)
         {
            x=bbox[0]+i*dx;
            if (dim>1) y=bbox[2]+j*dy; else y=ex_y;
            if (dim>2) z=bbox[4]+k*dz; else z=ex_z;

            if (sqrt((x-ex_x)*(x-ex_x)+(y-ex_y)*(y-ex_y)+(z-ex_z)*(z-ex_z))>ex_R) 
               mask[i+j*shape[0]+k*shape[0]*shape[1]]=NOT_EXCISED;
            else
               mask[i+j*shape[0]+k*shape[0]*shape[1]]=excised;
         }
}
  
  
  
  

// This will be the user-defined function for our tests.
double  user_func(double x0, double y0, double z0) {
  int profile_type = 0;                   // Test function profile.  1 for discontinuous, 0 for smooth.
  double xscale = 0.5;                    // Scale factors for test function.
  double yscale = 0.5;
  double zscale = 0.25;

  if (profile_type == 0) {
    return pow((x0/xscale),3)+pow((y0/yscale),3)+pow((z0/zscale),3) + 50.0;
    //return pow((y0/yscale),3) + 50.0;
  } else {
    // This is the shock profile
    if (y0 > 0.0) {
      return pow((x0/xscale),3)+pow((y0/yscale),3)+pow((z0/zscale),3) + 50.0;
    } else {
      return pow((x0/xscale),3)+pow((y0/yscale),3)+pow((z0/zscale),3) + 26.0;
    }
  }

}

struct availbox {
  double vol;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double zmin;
  double zmax;
  int unoccupied;
};

int main(int argc, char **argv)
{


   int rank,size,cnt=-1,shape_c[3],i,j,l,k,h,ind,h1,h2,h3,h4,h5,h6,nn,mn,sz;
   //   int shape[3]={65,33,33};
   //   int shape[3]={129,65,65};
   int shape[3]={65,65,65};
   //   int shape[3]={257,129,129};
   //char pre_tag[5] = "4_\0";
   char pre_tag[5] = "\0";
   //   real bbox[6]={-1,1,-2,2,-0.5,0.5};
   real bbox[6]={-1,1,-1,1,-0.5,0.5};
   real lambda,bbox_b[6];
   real bbox_list[6*MAX_GRIDS];
   real tmp_bbox_list[6*MAX_GRIDS];
   real merged_bbox_list[6*MAX_GRIDS];
   int lev_list[MAX_GRIDS];
   int levels=3;
   int rho_sp[3]={2,2,2},rho_tm[3]={2,2,4};
   int grho_sp[3],grho_tm[3];
   int ghost_width[3]={3,3,3},min_width[3]={9,7,13};
   int ghost_width6[6];
   int MG_min_cwidth[3]={5,4,7};
   int gghost_width[3],gmin_width[3];
   int gMG_min_cwidth[3],gMG_max_cwidth[3];
   real top_xy=1,top_xz=0.5,top_yz=1.5;
   real gtop_xy,gtop_xz,gtop_yz;
   real x1,x2,t,min_eff;
   int phys_bnd1[6]={PAMR_UNKNOWN,PAMR_UNKNOWN,PAMR_UNKNOWN,PAMR_UNKNOWN,PAMR_UNKNOWN,PAMR_UNKNOWN};
   int phys_bnd2[6]={PAMR_ODD,PAMR_ODD,PAMR_ODD,PAMR_ODD,PAMR_ODD,PAMR_ODD};
   int phys_bnd3[6]={PAMR_EVEN,PAMR_EVEN,PAMR_EVEN,PAMR_EVEN,PAMR_EVEN,PAMR_EVEN};
   int in_amrh,in_mgh,num_tl,amr_inject,amr_interp,amr_bdy_interp,amr_sync,mg_inject;
   int mg_interp,mg_sync,mg_noinj_to_amr,regrid_transfer,phys_bdy_type[6];
   int f1_n,f2_n,f3_n,f4_n,f5_n,f6_n,chr_n,chr_c_n;
   int ltrace=0,valid,ngfs;
   int ngl2,ngl3,seed,n,dim,gdm,nt,max_lev,coarsest;
   real dt,dx[3];
   real *x[3],*x_c[3],*gfs[500],x0,y0,z0;
   int periodic[3]={1,0,1},c_to_v,v_to_c,var_type;
   char *c_tag,*v_tag;
   char f1_name[256],f2_name[256],f6_name[256];
   char chr[256];
   char chr_c[256];
   double c_norm, v_norm;
   double c_norml, v_norml;
   int c_counter, v_counter;
   int c_counterl, v_counterl;

   struct availbox available[7];
   struct availbox avail3[8][7];
   double x_cut, y_cut, z_cut;
   int x_flag, y_flag, z_flag, n_exp;
   int i_lg, amghost, udiff, ldiff;
   double ufac, lfac, vol_lg;
   double xmin,xmax,ymin,ymax,zmin,zmax;
   int ilo, iup, jlo, jup, klo, kup;
   int wamrbdy_v, wamrbdy_c;
   int in_shadow, isx, isy, isz;
   int which_test;

   int SYNC_TEST=0,SAVE_SYNC_TEST=1;
   int C_TO_V_TEST=0,SAVE_C_TO_V_TEST=1;
   int V_TO_C_TEST=0,SAVE_V_TO_C_TEST=1;
   int INTERP_TEST=0,SAVE_INTERP_TEST=1;
   int BDY_INTERP_TEST=0,SAVE_BDY_INTERP_TEST=1;
   int INJECT_TEST=0,SAVE_INJECT_TEST=1;
   int REGRID_TEST=0,SAVE_REGRID_TEST=1;
   int MGH_TEST=0,SAVE_MGH_TEST=1;
   int EXCISE_TEST=0;
   int NESTED_1=1;
   int NESTED_2=0;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   if (rank==0) printf("done initializing MPI\n");

   if (argc<8 || argc>9)
   {
      printf("usage: %s which_test dim ngl2 ngl3 seed v_tag c_tag [min_eff=1]\n",argv[0]);
      goto cleanup;
   }
   which_test=atoi(argv[1]);
   dim=atoi(argv[2]);
   ngl2=atoi(argv[3]);
   ngl3=atoi(argv[4]);
   seed=atoi(argv[5]);
   srand(seed);
   v_tag=argv[6];
   c_tag=argv[7];
   if (argc==9) min_eff=atof(argv[8]); else min_eff=1;

   if (which_test==1) {
     SYNC_TEST=1;
   } else if (which_test==2) {
     C_TO_V_TEST=1;
   } else if (which_test==3) {
     V_TO_C_TEST=1;
   } else if (which_test==4) {
     INTERP_TEST=1;
   } else if (which_test==5) {
     BDY_INTERP_TEST=1;
   } else if (which_test==6) {
     INJECT_TEST=1;
   } else if (which_test==7) {
     REGRID_TEST=1;
   } else if (which_test==8) {
     MGH_TEST=1;
   } else {
     if (rank==0) printf("Test type not supported. Exiting. \n");
   }
   
   if (NESTED_1 || NESTED_2) {
     ngl2 = 1;
     ngl3 = 1;
   }
 
   MPI_Comm_size(MPI_COMM_WORLD,&size);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   if (rank==0) printf("MPI size =%i, (my) rank=%i\n",size,rank);

   PAMR_set_trace_lev(0);

   if (strlen(c_tag)==0 && strlen(v_tag)==0) cnt=PAMR_init_context(0,0,dim,shape,bbox,0,0,0);
   if (strlen(c_tag)!=0 && strlen(v_tag)==0) cnt=PAMR_init_context(0,0,dim,shape,bbox,0,v_tag,0);
   if (strlen(c_tag)==0 && strlen(v_tag)!=0) cnt=PAMR_init_context(0,0,dim,shape,bbox,0,0,c_tag);
   if (strlen(c_tag)!=0 && strlen(v_tag)!=0) cnt=PAMR_init_context(0,0,dim,shape,bbox,0,v_tag,c_tag);

   if (rank==0) printf("Init context returns %i\n",cnt);

   //   PAMR_set_periodic_bdy(periodic);

   // excision sphere parameters
   ex_R=0.25*((real)rand())/RAND_MAX+0.25;
   ex_x=1.0*((real)rand())/RAND_MAX-0.5;
   ex_y=1.0*((real)rand())/RAND_MAX-0.5;
   ex_z=1.0*((real)rand())/RAND_MAX-0.5;

   if (rank==0) printf("Excision sphere parameters: R=%lf, x=[%lf,%lf,%lf]\n",ex_R,ex_x,ex_y,ex_z);

   PAMR_set_interp_buffer(2);

   PAMR_set_lambda(0.5);
   PAMR_get_lambda(&lambda);
   IFL printf("lambda=%lf\n",lambda);

   PAMR_set_rho(rho_sp,rho_tm,3);
   PAMR_get_rho(grho_sp,grho_tm,3);
   IFL printf("rho_sp=[%i,%i,%i], rho_tm=[%i,%i,%i]\n",
          grho_sp[0],grho_sp[1],grho_sp[2],grho_tm[0],grho_tm[1],grho_tm[2]);

   PAMR_set_ghost_width(ghost_width);
   PAMR_get_ghost_width(gghost_width);
   IFL printf("ghost_width=[%i,%i,%i]\n",gghost_width[0],gghost_width[1],gghost_width[2]);

   PAMR_set_top_ratios(top_xy,top_xz,top_yz);
   PAMR_get_top_ratios(&gtop_xy,&gtop_xz,&gtop_yz);
   IFL printf("top_ratios=[%lf,%lf,%lf]\n",gtop_xy,gtop_xz,gtop_yz);

   PAMR_set_min_width(min_width);
   PAMR_get_min_width(gmin_width);
   IFL printf("min_width=[%i,%i,%i]\n",gmin_width[0],gmin_width[1],gmin_width[2]);

   PAMR_set_gdm(0);
   PAMR_get_gdm(&gdm);
   IFL printf("GDM=%i\n",gdm);

   PAMR_set_MG_coarse_width(MG_min_cwidth);
   PAMR_get_MG_coarse_width(gMG_min_cwidth);
   IFL printf("MG_min_width=[%i,%i,%i]\n",gMG_min_cwidth[0],gMG_min_cwidth[1],gMG_min_cwidth[2]);
   //   IFL printf("MG_max_width=[%i,%i,%i]\n",gMG_max_cwidth[0],gMG_max_cwidth[1],gMG_max_cwidth[2]);

   sprintf(f1_name,"%s%s","f1",v_tag);
   f1_n=PAMR_def_var_full(f1_name,1,1,4,PAMR_STRAIGHT_INJECT,PAMR_FOURTH_ORDER,PAMR_FOURTH_ORDER,PAMR_SYNC,
                         PAMR_STRAIGHT_INJECT,PAMR_FOURTH_ORDER,PAMR_SYNC,1,PAMR_FOURTH_ORDER,
                         PAMR_C_TO_V_SECOND_ORDER,PAMR_V_TO_C_SECOND_ORDER,phys_bnd1);

/*    f1_n=PAMR_def_var_full(f1_name,1,1,4,PAMR_STRAIGHT_INJECT,PAMR_SECOND_ORDER,PAMR_SECOND_ORDER,PAMR_SYNC, */
/*                          PAMR_STRAIGHT_INJECT,PAMR_SECOND_ORDER,PAMR_SYNC,1,PAMR_NO_TRANSFER, */
/*                          PAMR_C_TO_V_SECOND_ORDER,PAMR_V_TO_C_SECOND_ORDER,phys_bnd1); */
   //   only define as a conjugate var to f1 if having a ctag and vtag
   if (strlen(v_tag) || strlen(c_tag)) sprintf(f2_name,"%s%s","f1",c_tag); else sprintf(f1_name,"%s%s","f1",c_tag);
   f2_n=PAMR_def_var_full(f2_name,1,1,4,PAMR_NN_AVERAGE,PAMR_MC,PAMR_MC,PAMR_SYNC,
                         PAMR_NN_AVERAGE,PAMR_MC,PAMR_SYNC,1,PAMR_MC,
                         PAMR_C_TO_V_SECOND_ORDER,PAMR_V_TO_C_SECOND_ORDER,phys_bnd1);
/*    f2_n=PAMR_def_var_full(f2_name,1,1,4,PAMR_NN_AVERAGE,PAMR_FIRST_ORDER_CONS,PAMR_FIRST_ORDER_CONS,PAMR_SYNC, */
/*                          PAMR_NN_AVERAGE,PAMR_FIRST_ORDER_CONS,PAMR_SYNC,1,PAMR_FIRST_ORDER_CONS, */
/*                          PAMR_C_TO_V_SECOND_ORDER,PAMR_V_TO_C_SECOND_ORDER,phys_bnd1); */


   // MAKE SURE MASK IS LAST GFN
   sprintf(chr,"%s%s","chr",v_tag);
   chr_n=PAMR_def_var_full(chr,1,1,1,PAMR_NO_INJECT,PAMR_NO_INTERP,PAMR_NO_INTERP,PAMR_NO_SYNC,
                            PAMR_NO_INJECT,PAMR_NO_INTERP,PAMR_NO_SYNC,1,PAMR_NO_TRANSFER,
                            PAMR_C_TO_V_NO_TRANSFER,PAMR_V_TO_C_NO_TRANSFER,phys_bnd3);
   chr_c[0]=0;
   if (strlen(v_tag) || strlen(c_tag))
   {
      sprintf(chr_c,"%s%s","chr",c_tag);
      chr_c_n=PAMR_def_var_full(chr_c,1,1,1,PAMR_NO_INJECT,PAMR_NO_INTERP,PAMR_NO_INTERP,PAMR_NO_SYNC,
                                PAMR_NO_INJECT,PAMR_NO_INTERP,PAMR_NO_SYNC,1,PAMR_NO_TRANSFER,
                                PAMR_C_TO_V_NO_TRANSFER,PAMR_V_TO_C_NO_TRANSFER,phys_bnd3);
   }

   if (EXCISE_TEST)
   {
      if (rank==0) printf("testing excision, chr_n=%i\n",chr_n);
      PAMR_excision_on(chr,chr_c,set_ex_mask,EXCISED,1);
   }

/*    IFL printf("starting grid function numbers: f1,f2,f3,f4,f5,f6:%i,%i,%i,%i,%i,%i\n", */
/*           f1_n,f2_n,f3_n,f4_n,f5_n,f6_n); */

   // ******************** Setup initial grid hierarchy ******************* // 

   if (NESTED_1) {
     lev_list[0]=1;
     for (n=0; n<6; n++) bbox_list[n]=bbox[n];

     n = 1;
     x1 = 0.25;
     x2 = 0.75;
     PAMR_get_dxdt(1,&dx[0],&dt);

     bbox_list[2*dim*n]=max((int)(x1*(bbox[1]-bbox[0])/dx[0])*dx[0]+bbox[0],bbox[0]);
     bbox_list[2*dim*n+1]=min((int)(x2*(bbox[1]-bbox[0])/dx[0])*dx[0]+bbox[0],bbox[1]);
     if (dim > 1) {
       bbox_list[2*dim*n+2]=max((int)(x1*(bbox[3]-bbox[2])/dx[1])*dx[1]+bbox[2],bbox[2]);
       bbox_list[2*dim*n+3]=min((int)(x2*(bbox[3]-bbox[2])/dx[1])*dx[1]+bbox[2],bbox[3]);
       if (dim > 2) {
	 bbox_list[2*dim*n+4]=max((int)(x1*(bbox[5]-bbox[4])/dx[2])*dx[2]+bbox[4],bbox[4]);
	 bbox_list[2*dim*n+5]=min((int)(x2*(bbox[5]-bbox[4])/dx[2])*dx[2]+bbox[4],bbox[5]);
       }
     }

     lev_list[n] = 2;

     n = 2;
     x1 = 0.375;
     x2 = 0.625;
     PAMR_get_dxdt(2,&dx[0],&dt);

     bbox_list[2*dim*n]=max((int)(x1*(bbox[1]-bbox[0])/dx[0])*dx[0]+bbox[0],bbox[0]);
     bbox_list[2*dim*n+1]=min((int)(x2*(bbox[1]-bbox[0])/dx[0])*dx[0]+bbox[0],bbox[1]);
     if (dim > 1) {
       bbox_list[2*dim*n+2]=max((int)(x1*(bbox[3]-bbox[2])/dx[1])*dx[1]+bbox[2],bbox[2]);
       bbox_list[2*dim*n+3]=min((int)(x2*(bbox[3]-bbox[2])/dx[1])*dx[1]+bbox[2],bbox[3]);
       if (dim > 2) {
	 bbox_list[2*dim*n+4]=max((int)(x1*(bbox[5]-bbox[4])/dx[1])*dx[1]+bbox[4],bbox[4]);
	 bbox_list[2*dim*n+5]=min((int)(x2*(bbox[5]-bbox[4])/dx[1])*dx[1]+bbox[4],bbox[5]);
       }
     }

     lev_list[n] = 3;

   } else if (NESTED_2) {

     lev_list[0]=1;
     for (n=0; n<6; n++) bbox_list[n]=bbox[n];

     n = 1;
     PAMR_get_dxdt(1,&dx[0],&dt);
     x1 = 0.25;
     x2 = 0.75;
     bbox_list[2*dim*n]=max((int)(x1*(bbox[1]-bbox[0])/dx[0])*dx[0]+bbox[0],bbox[0]);
     bbox_list[2*dim*n+1]=min((int)(x2*(bbox[1]-bbox[0])/dx[0])*dx[0]+bbox[0],bbox[1]);
     x1 = 0.5;
     x2 = 1.0;
     if (dim > 1) {
       bbox_list[2*dim*n+2]=max((int)(x1*(bbox[3]-bbox[2])/dx[1])*dx[1]+bbox[2],bbox[2]);
       bbox_list[2*dim*n+3]=min((int)(x2*(bbox[3]-bbox[2])/dx[1])*dx[1]+bbox[2],bbox[3]);
       if (dim > 2) {
	 bbox_list[2*dim*n+4]=max((int)(x1*(bbox[5]-bbox[4])/dx[2])*dx[2]+bbox[4],bbox[4]);
	 bbox_list[2*dim*n+5]=min((int)(x2*(bbox[5]-bbox[4])/dx[2])*dx[2]+bbox[4],bbox[5]);
       }
     }
     lev_list[n] = 2;

     n = 2;
     PAMR_get_dxdt(2,&dx[0],&dt);
     x1 = 0.375;
     x2 = 0.625;
     bbox_list[2*dim*n]=max((int)(x1*(bbox[1]-bbox[0])/dx[0])*dx[0]+bbox[0],bbox[0]);
     bbox_list[2*dim*n+1]=min((int)(x2*(bbox[1]-bbox[0])/dx[0])*dx[0]+bbox[0],bbox[1]);
     x1 = 0.75;
     x2 = 1.0;
     if (dim > 1) {
       bbox_list[2*dim*n+2]=max((int)(x1*(bbox[3]-bbox[2])/dx[1])*dx[1]+bbox[2],bbox[2]);
       bbox_list[2*dim*n+3]=min((int)(x2*(bbox[3]-bbox[2])/dx[1])*dx[1]+bbox[2],bbox[3]);
       if (dim > 2) {
	 bbox_list[2*dim*n+4]=max((int)(x1*(bbox[5]-bbox[4])/dx[2])*dx[2]+bbox[4],bbox[4]);
	 bbox_list[2*dim*n+5]=min((int)(x2*(bbox[5]-bbox[4])/dx[2])*dx[2]+bbox[4],bbox[5]);
       }
     }
     lev_list[n] = 3;

   } else {

     // Compose a random set of grids
     // Prohibit overlapping grids on the same level.
     // Try to prevent any grids from being less than 4 cells wide.
     x_cut = 0.0;
     y_cut = 0.0;
     z_cut = 0.0;
     
     lev_list[0]=1;
     for (n=0; n<6; n++) bbox_list[n]=bbox[n];

     // Create level 2 grids.

     for (n=1; n<(1+ngl2); n++)
       {
	 PAMR_get_dxdt(1,&dx[0],&dt);
	 
	 if (n==1) {
	   // Set reference bbox to global.
	   xmin = bbox[0];
	   xmax = bbox[1];	
	   ymin = bbox[2];
	   ymax = bbox[3];
	   zmin = bbox[4];
	   zmax = bbox[5];
	 } else {
	   // set reference bbox to largest available subregion
	   i_lg = 0;
	   vol_lg = -1.0;
	   if (dim == 1) {
	     xmin = available[0].xmin;
	     xmax = available[0].xmax;
	     ymin = bbox[2];
	     ymax = bbox[3];
	     zmin = bbox[4];
	     zmax = bbox[5];
	     available[0].unoccupied = 0;
	   } else if (dim == 2) {
	     for (i=0;i<3;i++) {
	       if (available[i].vol*available[i].unoccupied > vol_lg) {
		 vol_lg = available[i].vol*available[i].unoccupied;
		 i_lg = i;
	       }
	     }
	     xmin = available[i_lg].xmin;
	     xmax = available[i_lg].xmax;
	     ymin = available[i_lg].ymin;
	     ymax = available[i_lg].ymax;
	     zmin = bbox[4];
	     zmax = bbox[5];
	     available[i_lg].unoccupied = 0;

	   } else if (dim == 3) {
	     for (i=0;i<7;i++) {
	       if (available[i].vol*available[i].unoccupied > vol_lg) {
		 vol_lg = available[i].vol*available[i].unoccupied;
		 i_lg = i;
	       }
	     }
	     xmin = available[i_lg].xmin;
	     xmax = available[i_lg].xmax;
	     ymin = available[i_lg].ymin;
	     ymax = available[i_lg].ymax;
	     zmin = available[i_lg].zmin;
	     zmax = available[i_lg].zmax;
	     available[i_lg].unoccupied = 0;
	   }
	 }
	 
	 x1=((real)rand())/RAND_MAX;
	 x2=((real)rand())/RAND_MAX;
	 if (x1>x2) {t=x1; x1=x2; x2=t;}
	 bbox_list[2*dim*n]=max((int)(x1*(xmax-xmin)/dx[0])*dx[0]+xmin,xmin+2*dx[0]);
	 bbox_list[2*dim*n+1]=min((int)(x2*(xmax-xmin)/dx[0])*dx[0]+xmin,xmax-2*dx[0]);

	 if (bbox_list[2*dim*n]>bbox_list[2*dim*n+1]) {
	   t=bbox_list[2*dim*n];
	   bbox_list[2*dim*n]   = bbox_list[2*dim*n+1];
	   bbox_list[2*dim*n+1] = t;
	 }

	 // Try to prevent boxes which are too small in the x direction
	 if ((bbox_list[2*dim*n+1] -  bbox_list[2*dim*n]) < 4*dx[0]) {
	   n_exp = 4 - (bbox_list[2*dim*n+1]-bbox_list[2*dim*n])/dx[0];
	   n_exp = (n_exp + 1)/2;
	   udiff = (xmax - 2*dx[0] - bbox_list[2*dim*n+1])/dx[0] + 0.01;
	   ldiff = (bbox_list[2*dim*n] - xmin - 2*dx[0])/dx[0] + 0.01;
	   lfac = 0.0;
	   ufac = 0.0;
	   if (udiff >= n_exp && ldiff >= n_exp) {
	     lfac = n_exp;
	     ufac = n_exp;
	   } else if (udiff < n_exp && ldiff >= n_exp) {
	     lfac = min(ldiff,(2*n_exp-udiff));
	     ufac = udiff;
	   } else if (udiff >= n_exp && ldiff < n_exp) {
	     lfac = ldiff;
	     ufac = min(udiff,(2*n_exp-ldiff));
	   } else if (udiff < n_exp && ldiff < n_exp) {
	     // expand as much as you can
	     lfac = max(ldiff,0);
	     ufac = max(udiff,0);
	   }
	   // Try to expand it
	   bbox_list[2*dim*n] -= lfac*dx[0];
	   bbox_list[2*dim*n+1] += ufac*dx[0];
	 }

	 // For good measure....
	 bbox_list[2*dim*n+2] = ymin;
	 bbox_list[2*dim*n+3] = ymax;
	 bbox_list[2*dim*n+4] = zmin;
	 bbox_list[2*dim*n+5] = zmax;
	 if (dim>1)
	   {
	     x1=((real)rand())/RAND_MAX;
	     x2=((real)rand())/RAND_MAX;
	     if (x1>x2) {t=x1; x1=x2; x2=t;}
	     bbox_list[2*dim*n+2]=max((int)(x1*(ymax-ymin)/dx[1])*dx[1]+ymin,ymin+2*dx[1]);
	     bbox_list[2*dim*n+3]=min((int)(x2*(ymax-ymin)/dx[1])*dx[1]+ymin,ymax-2*dx[1]);

	     if (bbox_list[2*dim*n+2]>bbox_list[2*dim*n+3]) {
	       t=bbox_list[2*dim*n+2];
	       bbox_list[2*dim*n+2] = bbox_list[2*dim*n+3];
	       bbox_list[2*dim*n+3] = t;
	     }

	     // Try to prevent boxes which are too small in the y direction
	     if ((bbox_list[2*dim*n+3] -  bbox_list[2*dim*n+2]) < 4*dx[1]) {
	       n_exp = 4 - (bbox_list[2*dim*n+3]-bbox_list[2*dim*n+2])/dx[1];
	       n_exp = (n_exp + 1)/2;
	       udiff = (ymax - 2*dx[1] - bbox_list[2*dim*n+3])/dx[1] + 0.01;
	       ldiff = (bbox_list[2*dim*n+2] - ymin - 2*dx[1])/dx[1] + 0.01;
	       lfac = 0.0;
	       ufac = 0.0;
	       if (udiff >= n_exp && ldiff >= n_exp) {
		 lfac = n_exp;
		 ufac = n_exp;
	       } else if (udiff < n_exp && ldiff >= n_exp) {
		 lfac = min(ldiff,(2*n_exp-udiff));
		 ufac = udiff;
	       } else if (udiff >= n_exp && ldiff < n_exp) {
		 lfac = ldiff;
		 ufac = min(udiff,(2*n_exp-ldiff));
	       } else if (udiff < n_exp && ldiff < n_exp) {
		 // expand as much as you can
		 lfac = max(ldiff,0);
		 ufac = max(udiff,0);
	       }
	       // Try to expand it
	       bbox_list[2*dim*n+2] -= lfac*dx[1];
	       bbox_list[2*dim*n+3] += ufac*dx[1];
	     }

	   }
	 if (dim>2)
	   {
	     x1=((real)rand())/RAND_MAX;
	     x2=((real)rand())/RAND_MAX;
	     if (x1>x2) {t=x1; x1=x2; x2=t;}
	     bbox_list[2*dim*n+4]=max((int)(x1*(zmax-zmin)/dx[2])*dx[2]+zmin,zmin+2*dx[2]);
	     bbox_list[2*dim*n+5]=min((int)(x2*(zmax-zmin)/dx[2])*dx[2]+zmin,zmax-2*dx[2]);

	     if (bbox_list[2*dim*n+4]>bbox_list[2*dim*n+5]) {
	       t=bbox_list[2*dim*n+4];
	       bbox_list[2*dim*n+4] = bbox_list[2*dim*n+5];
	       bbox_list[2*dim*n+5] = t;
	     }

	     // Try to prevent boxes which are too small in the z direction
	     if ((bbox_list[2*dim*n+5] -  bbox_list[2*dim*n+4]) < 4*dx[2]) {
	       n_exp = 4 - (bbox_list[2*dim*n+5]-bbox_list[2*dim*n+4])/dx[2];
	       n_exp = (n_exp + 1)/2;
	       udiff = (zmax - 2*dx[2] - bbox_list[2*dim*n+5])/dx[2] + 0.01;
	       ldiff = (bbox_list[2*dim*n+4] - zmin - 2*dx[2])/dx[2] + 0.01;
	       lfac = 0.0;
	       ufac = 0.0;
	       if (udiff >= n_exp && ldiff >= n_exp) {
		 lfac = n_exp;
		 ufac = n_exp;
	       } else if (udiff < n_exp && ldiff >= n_exp) {
		 lfac = min(ldiff,(2*n_exp-udiff));
		 ufac = udiff;
	       } else if (udiff >= n_exp && ldiff < n_exp) {
		 lfac = ldiff;
		 ufac = min(udiff,(2*n_exp-ldiff));
	       } else if (udiff < n_exp && ldiff < n_exp) {
		 // expand as much as you can
		 lfac = max(ldiff,0);
		 ufac = max(udiff,0);
	       }
	       // Try to expand it
	       bbox_list[2*dim*n+4] -= lfac*dx[2];
	       bbox_list[2*dim*n+5] += ufac*dx[2];
	     }
	   }
	 
	 lev_list[n]=2;
	 
	 if (rank==0) {
	   printf("n=%d, bbox = %f, %f, %f, %f, %f, %f\n",n,bbox_list[2*dim*n],bbox_list[2*dim*n+1],
		  bbox_list[2*dim*n+2],bbox_list[2*dim*n+3],bbox_list[2*dim*n+4],bbox_list[2*dim*n+5]);
	 }

	 if (n==1) {
	   // Figure out x_cut, y_cut, and z_cut.
	   if ((bbox_list[2*dim*n]-xmin)/(xmax-xmin) > (xmax-bbox_list[2*dim*n+1])/(xmax-xmin)) {
	     x_cut = bbox_list[2*dim*n];
	     x_flag = 0;
	   } else {
	     x_cut = bbox_list[2*dim*n+1];
	     x_flag = 1;
	   }
	   
	   if (dim > 1) {
	     if ((bbox_list[2*dim*n+2]-ymin)/(ymax-ymin) > (ymax-bbox_list[2*dim*n+3])/(ymax-ymin)) {
	       y_cut = bbox_list[2*dim*n+2];
	       y_flag = 0;
	     } else {
	       y_cut = bbox_list[2*dim*n+3];
	       y_flag = 1;
	     }
	   }

	   if (dim > 2) {
	     if ((bbox_list[2*dim*n+4]-zmin)/(zmax-zmin) > (zmax-bbox_list[2*dim*n+5])/(zmax-zmin)) {
	       z_cut = bbox_list[2*dim*n+4];
	       z_flag = 0;
	     } else {
	       z_cut = bbox_list[2*dim*n+5];
	       z_flag = 1;
	     }
	   }
	 
	   // Make list of available boxes.
	   if (dim == 1) {
	     // One dimension -> one available box
	     available[0].unoccupied = 1;
	     available[0].xmin = ( x_flag==0 ? x_cut : xmin  );
	     available[0].xmax = ( x_flag==0 ? xmax  : x_cut );
	     available[0].ymin = 0.0;
	     available[0].ymax = 0.0;
	     available[0].zmin = 0.0;
	     available[0].zmax = 0.0;
	     available[0].vol = available[0].xmax-available[0].xmin;
	   } else if (dim == 2) {
	     // Two dimensions -> three available boxes

	     for (i=0; i<3; i++) available[i].unoccupied = 1;
	     for (i=0; i<3; i++) available[i].zmin = 0.0;
	     for (i=0; i<3; i++) available[i].zmax = 0.0;

	     available[0].xmin = ( x_flag==1 ? x_cut : xmin  );
	     available[0].xmax = ( x_flag==1 ? xmax  : x_cut );
	     available[1].xmin = ( x_flag==1 ? x_cut : xmin  );
	     available[1].xmax = ( x_flag==1 ? xmax  : x_cut );
	     available[2].xmin = ( x_flag==1 ? xmin  : x_cut );
	     available[2].xmax = ( x_flag==1 ? x_cut : xmax  );

	     available[0].ymin = ymin;
	     available[0].ymax = y_cut;
	     available[1].ymin = y_cut;
	     available[1].ymax = ymax;
	     available[2].ymin = ( y_flag==1 ? y_cut : ymin  );
	     available[2].ymax = ( y_flag==1 ? ymax  : y_cut );

	     // Calculate volumes
	     for (i=0;i<3;i++) available[i].vol = (available[i].xmax-available[i].xmin) *
                                	         (available[i].ymax-available[i].ymin);

	   } else {
	     // Three dimensions -> seven available boxes
	     for (i=0; i<7; i++) available[i].unoccupied = 1;

	     // The group of 4
	     available[0].xmin = xmin;
	     available[0].xmax = x_cut;
	     available[0].ymin = ymin;
	     available[0].ymax = y_cut;

	     available[1].xmin = xmin;
	     available[1].xmax = x_cut;
	     available[1].ymin = y_cut;
	     available[1].ymax = ymax;

	     available[2].xmin = x_cut;
	     available[2].xmax = xmax;
	     available[2].ymin = ymin;
	     available[2].ymax = y_cut;

	     available[3].xmin = x_cut;
	     available[3].xmax = xmax;
	     available[3].ymin = y_cut;
	     available[3].ymax = ymax;

	     for (i=0; i<4; i++) {
	       available[i].zmin = ( z_flag==1 ? z_cut : zmin  );
	       available[i].zmax = ( z_flag==1 ? zmax  : z_cut );
	     } 

	     // The group of 3
	     available[4].xmin = ( x_flag==1 ? x_cut : xmin  );
	     available[4].xmax = ( x_flag==1 ? xmax  : x_cut );
	     available[5].xmin = ( x_flag==1 ? x_cut : xmin  );
	     available[5].xmax = ( x_flag==1 ? xmax  : x_cut );
	     available[6].xmin = ( x_flag==1 ? xmin  : x_cut );
	     available[6].xmax = ( x_flag==1 ? x_cut : xmax  );

	     available[4].ymin = ymin;
	     available[4].ymax = y_cut;
	     available[5].ymin = y_cut;
	     available[5].ymax = ymax;
	     available[6].ymin = ( y_flag==1 ? y_cut : ymin );
	     available[6].ymax = ( y_flag==1 ? ymax  : y_cut);

	     for (i=4; i<7; i++) {
	       available[i].zmin = ( z_flag==1 ? zmin  : z_cut );
	       available[i].zmax = ( z_flag==1 ? z_cut : zmax  );
	     }

	     // Calculate volumes
	     for (i=0;i<7;i++) available[i].vol = (available[i].xmax-available[i].xmin) *
                                        	 (available[i].ymax-available[i].ymin) *
                                	         (available[i].zmax-available[i].zmin);
	   }
	 } // end if n==1
       } // end loop on n for level 2 grids.

     // Create level 3 grids.
     for (n=ngl2+1; n<(1+ngl2+ngl3); n++)
       {
	 PAMR_get_dxdt(2,&dx[0],&dt);
	 // Change reference bbox to appropriate level 2 box.
	 i = (n-1)%ngl2 + 1;

	 if (n < (2*ngl2 + 1)) {
	   xmin = bbox_list[2*dim*i];
	   xmax = bbox_list[2*dim*i+1];
	   ymin = bbox_list[2*dim*i+2];
	   ymax = bbox_list[2*dim*i+3];
	   zmin = bbox_list[2*dim*i+4];
	   zmax = bbox_list[2*dim*i+5];
	 } else { // use largest remaining available region
	   i_lg = 0;
	   vol_lg = 0.0;
	   if (dim == 1) {
	     xmin = avail3[i][0].xmin;
	     xmax = avail3[i][0].xmax;
	     available[0].unoccupied = 0;
	   } else if (dim == 2) {
	     for (j=0;j<3;j++) {
	       if (avail3[i][j].vol*avail3[i][j].unoccupied > vol_lg) {
		 vol_lg = avail3[i][j].vol*avail3[i][j].unoccupied;
		 i_lg = j;
	       }
	     }
	     
	     xmin = avail3[i][i_lg].xmin;
	     xmax = avail3[i][i_lg].xmax;
	     ymin = avail3[i][i_lg].ymin;
	     ymax = avail3[i][i_lg].ymax;
	     avail3[i][i_lg].unoccupied = 0;

	   } else if (dim == 3) {
	     for (j=0;j<7;j++) {
	       if (avail3[i][j].vol*avail3[i][j].unoccupied > vol_lg) {
		 vol_lg = avail3[i][j].vol*avail3[i][j].unoccupied;
		 i_lg = i;
	       }
	     }
	     xmin = avail3[i][i_lg].xmin;
	     xmax = avail3[i][i_lg].xmax;
	     ymin = avail3[i][i_lg].ymin;
	     ymax = avail3[i][i_lg].ymax;
	     zmin = avail3[i][i_lg].zmin;
	     zmax = avail3[i][i_lg].zmax;
	     avail3[i][i_lg].unoccupied = 0;
	   }
	 }
	 
	 x1=((real)rand())/RAND_MAX;
	 x2=((real)rand())/RAND_MAX;
	 if (x1>x2) {t=x1; x1=x2; x2=t;}
	 bbox_list[2*dim*n]=max((int)(x1*(xmax-xmin)/dx[0])*dx[0]+xmin,xmin+2*dx[0]);
	 bbox_list[2*dim*n+1]=min((int)(x2*(xmax-xmin)/dx[0])*dx[0]+xmin,xmax-2*dx[0]);


	 if (bbox_list[2*dim*n]>bbox_list[2*dim*n+1]) {
	   t=bbox_list[2*dim*n];
	   bbox_list[2*dim*n]   = bbox_list[2*dim*n+1];
	   bbox_list[2*dim*n+1] = t;
	 }

	 if ((bbox_list[2*dim*n+1] -  bbox_list[2*dim*n]) < 4*dx[0]) {
	   n_exp = 4 - (bbox_list[2*dim*n+1]-bbox_list[2*dim*n])/dx[0];
	   n_exp = (n_exp + 1)/2;
	   udiff = (xmax - 2*dx[0] - bbox_list[2*dim*n+1])/dx[0] + 0.01;
	   ldiff = (bbox_list[2*dim*n] - xmin - 2*dx[0])/dx[0] + 0.01;
	   lfac = 0.0;
	   ufac = 0.0;
	   if (udiff >= n_exp && ldiff >= n_exp) {
	     lfac = n_exp;
	     ufac = n_exp;
	   } else if (udiff < n_exp && ldiff >= n_exp) {
	     lfac = min(ldiff,(2*n_exp-udiff));
	     ufac = udiff;
	   } else if (udiff >= n_exp && ldiff < n_exp) {
	     lfac = ldiff;
	     ufac = min(udiff,(2*n_exp-ldiff));
	   } else if (udiff < n_exp && ldiff < n_exp) {
	     // expand as much as you can
	     lfac = max(ldiff,0);
	     ufac = max(udiff,0);
	   }
	   // Try to expand it
	   bbox_list[2*dim*n] -= lfac*dx[0];
	   bbox_list[2*dim*n+1] += ufac*dx[0];
	 }

	 // For good measure....
	 bbox_list[2*dim*n+2] = ymin;
	 bbox_list[2*dim*n+3] = ymax;
	 bbox_list[2*dim*n+4] = zmin;
	 bbox_list[2*dim*n+5] = zmax;
	 if (dim>1)
	   {
	     x1=((real)rand())/RAND_MAX;
	     x2=((real)rand())/RAND_MAX;
	     if (x1>x2) {t=x1; x1=x2; x2=t;}
	     bbox_list[2*dim*n+2]=max((int)(x1*(ymax-ymin)/dx[1])*dx[1]+ymin,ymin+2*dx[1]);
	     bbox_list[2*dim*n+3]=min((int)(x2*(ymax-ymin)/dx[1])*dx[1]+ymin,ymax-2*dx[1]);

	     if (bbox_list[2*dim*n+2]>bbox_list[2*dim*n+3]) {
	       t=bbox_list[2*dim*n+2];
	       bbox_list[2*dim*n+2] = bbox_list[2*dim*n+3];
	       bbox_list[2*dim*n+3] = t;
	     }

	     if ((bbox_list[2*dim*n+3] -  bbox_list[2*dim*n+2]) < 4*dx[1]) {
	       n_exp = 4 - (bbox_list[2*dim*n+3]-bbox_list[2*dim*n+2])/dx[1];
	       n_exp = (n_exp + 1)/2;
	       udiff = (ymax - 2*dx[1] - bbox_list[2*dim*n+3])/dx[1] + 0.01;
	       ldiff = (bbox_list[2*dim*n+2] - ymin - 2*dx[1])/dx[1] + 0.01;
	       lfac = 0.0;
	       ufac = 0.0;
	       if (udiff >= n_exp && ldiff >= n_exp) {
		 lfac = n_exp;
		 ufac = n_exp;
	       } else if (udiff < n_exp && ldiff >= n_exp) {
		 lfac = min(ldiff,(2*n_exp-udiff));
		 ufac = udiff;
	       } else if (udiff >= n_exp && ldiff < n_exp) {
		 lfac = ldiff;
		 ufac = min(udiff,(2*n_exp-ldiff));
	       } else if (udiff < n_exp && ldiff < n_exp) {
		 // expand as much as you can
		 lfac = max(ldiff,0);
		 ufac = max(udiff,0);
	       }
	       // Try to expand it
	       bbox_list[2*dim*n+2] -= lfac*dx[1];
	       bbox_list[2*dim*n+3] += ufac*dx[1];
	     }

	   }
	 if (dim>2)
	   {
	     x1=((real)rand())/RAND_MAX;
	     x2=((real)rand())/RAND_MAX;
	     if (x1>x2) {t=x1; x1=x2; x2=t;}
	     bbox_list[2*dim*n+4]=max((int)(x1*(zmax-zmin)/dx[2])*dx[2]+zmin,zmin+2*dx[2]);
	     bbox_list[2*dim*n+5]=min((int)(x2*(zmax-zmin)/dx[2])*dx[2]+zmin,zmax-2*dx[2]);

	     if (bbox_list[2*dim*n+4]>bbox_list[2*dim*n+5]) {
	       t=bbox_list[2*dim*n+4];
	       bbox_list[2*dim*n+4] = bbox_list[2*dim*n+5];
	       bbox_list[2*dim*n+5] = t;
	     }

	     // Try to prevent boxes which are too small in the z direction
	     if ((bbox_list[2*dim*n+5] -  bbox_list[2*dim*n+4]) < 4*dx[2]) {
	       n_exp = 4 - (bbox_list[2*dim*n+5]-bbox_list[2*dim*n+4])/dx[2];
	       n_exp = (n_exp + 1)/2;
	       udiff = (zmax - 2*dx[2] - bbox_list[2*dim*n+5])/dx[2] + 0.01;
	       ldiff = (bbox_list[2*dim*n+4] - zmin - 2*dx[2])/dx[2] + 0.01;
	       lfac = 0.0;
	       ufac = 0.0;
	       if (udiff >= n_exp && ldiff >= n_exp) {
		 lfac = n_exp;
		 ufac = n_exp;
	       } else if (udiff < n_exp && ldiff >= n_exp) {
		 lfac = min(ldiff,(2*n_exp-udiff));
		 ufac = udiff;
	       } else if (udiff >= n_exp && ldiff < n_exp) {
		 lfac = ldiff;
		 ufac = min(udiff,(2*n_exp-ldiff));
	       } else if (udiff < n_exp && ldiff < n_exp) {
		 // expand as much as you can
		 lfac = max(ldiff,0);
		 ufac = max(udiff,0);
	       }
	       // Try to expand it
	       bbox_list[2*dim*n+4] -= lfac*dx[2];
	       bbox_list[2*dim*n+5] += ufac*dx[2];
	     }

	   }
	 
	 lev_list[n]=3; 
	 
	 if (rank==0) {
	   printf("n=%d, bbox = %f, %f, %f, %f, %f, %f\n",n,bbox_list[2*dim*n],bbox_list[2*dim*n+1],
		  bbox_list[2*dim*n+2],bbox_list[2*dim*n+3],bbox_list[2*dim*n+4],bbox_list[2*dim*n+5]);
	 }

	 // If this is the first visit to the given level 2 grid, then 
	 // calculate the available boxes.
	 if (n < (2*ngl2 + 1)) {

	   // Figure out x_cut, y_cut, and z_cut.
	   if ((bbox_list[2*dim*n]-xmin)/(xmax-xmin) > (xmax-bbox_list[2*dim*n+1])/(xmax-xmin)) {
	     x_cut = bbox_list[2*dim*n];
	     x_flag = 0;
	   } else {
	     x_cut = bbox_list[2*dim*n+1];
	     x_flag = 1;
	   }
	   
	   if (dim > 1) {
	     if ((bbox_list[2*dim*n+2]-ymin)/(ymax-ymin) > (ymax-bbox_list[2*dim*n+3])/(ymax-ymin)) {
	       y_cut = bbox_list[2*dim*n+2];
	       y_flag = 0;
	     } else {
	       y_cut = bbox_list[2*dim*n+3];
	       y_flag = 1;
	     }
	   }

	   if (dim > 2) {
	     if ((bbox_list[2*dim*n+4]-zmin)/(zmax-zmin) > (zmax-bbox_list[2*dim*n+5])/(zmax-zmin)) {
	       z_cut = bbox_list[2*dim*n+4];
	       z_flag = 0;
	     } else {
	       z_cut = bbox_list[2*dim*n+5];
	       z_flag = 1;
	     }
	   }
	 
	   // Make list of available boxes.
	   if (dim == 1) {
	     // One dimension -> one available box
	     avail3[i][0].unoccupied = 1;

	     avail3[i][0].xmin = ( x_flag == 0 ? x_cut : xmin  );
	     avail3[i][0].xmax = ( x_flag == 0 ? xmax  : x_cut );
	     avail3[i][0].ymin = 0.0;
	     avail3[i][0].ymax = 0.0;
	     avail3[i][0].zmin = 0.0;
	     avail3[i][0].zmax = 0.0;
	     avail3[i][0].vol = avail3[i][0].xmax-avail3[i][0].xmin;
	   } else if (dim == 2) {
	     // Two dimensions -> three available boxes

	     for (j=0; j<3; j++) avail3[i][j].unoccupied = 1;
	     for (j=0; j<3; j++) avail3[i][j].zmin = 0.0;
	     for (j=0; j<3; j++) avail3[i][j].zmax = 0.0;

	     avail3[i][0].xmin = ( x_flag==1 ? x_cut : xmin  );
	     avail3[i][0].xmax = ( x_flag==1 ? xmax  : x_cut );
	     avail3[i][1].xmin = ( x_flag==1 ? x_cut : xmin  );
	     avail3[i][1].xmax = ( x_flag==1 ? xmax  : x_cut );
	     avail3[i][2].xmin = ( x_flag==1 ? xmin  : x_cut );
	     avail3[i][2].xmax = ( x_flag==1 ? x_cut : xmax  );

	     avail3[i][0].ymin = ymin;
	     avail3[i][0].ymax = y_cut;
	     avail3[i][1].ymin = y_cut;
	     avail3[i][1].ymax = ymax;
	     avail3[i][2].ymin = ( y_flag==1 ? y_cut : ymin  );
	     avail3[i][2].ymax = ( y_flag==1 ? ymax  : y_cut );

	     // Calculate volumes
	     for (j=0;j<3;j++) avail3[i][j].vol = (avail3[i][j].xmax-avail3[i][j].xmin) *
                                	         (avail3[i][j].ymax-avail3[i][j].ymin);

	   } else {
	     // Three dimensions -> seven available boxes
	     for (j=0; j<7; j++) avail3[i][j].unoccupied = 1;

	     // The group of 4
	     avail3[i][0].xmin = xmin;
	     avail3[i][0].xmax = x_cut;
	     avail3[i][0].ymin = ymin;
	     avail3[i][0].ymax = y_cut;

	     avail3[i][1].xmin = xmin;
	     avail3[i][1].xmax = x_cut;
	     avail3[i][1].ymin = y_cut;
	     avail3[i][1].ymax = ymax;

	     avail3[i][2].xmin = x_cut;
	     avail3[i][2].xmax = xmax;
	     avail3[i][2].ymin = ymin;
	     avail3[i][2].ymax = y_cut;

	     avail3[i][3].xmin = x_cut;
	     avail3[i][3].xmax = xmax;
	     avail3[i][3].ymin = y_cut;
	     avail3[i][3].ymax = ymax;

	     for (j=0; j<4; j++) {
	       avail3[i][j].zmin = ( z_flag==1 ? z_cut : zmin  );
	       avail3[i][j].zmax = ( z_flag==1 ? zmax  : z_cut );
	     }

	     // The group of 3
	     avail3[i][4].xmin = ( x_flag==1 ? x_cut : xmin  );
	     avail3[i][4].xmax = ( x_flag==1 ? xmax  : x_cut );
	     avail3[i][5].xmin = ( x_flag==1 ? x_cut : xmin  );
	     avail3[i][5].xmax = ( x_flag==1 ? xmax  : x_cut );
	     avail3[i][6].xmin = ( x_flag==1 ? xmin  : x_cut );
	     avail3[i][6].xmax = ( x_flag==1 ? x_cut : xmax  );

	     avail3[i][4].ymin = ymin;
	     avail3[i][4].ymax = y_cut;
	     avail3[i][5].ymin = y_cut;
	     avail3[i][5].ymax = ymax;
	     avail3[i][6].ymin = ( y_flag==1 ? y_cut : ymin  );
	     avail3[i][6].ymax = ( y_flag==1 ? ymax  : y_cut );

	     for (j=4; j<7; j++) {
	       avail3[i][j].zmin = ( z_flag==1 ? zmin  : z_cut );
	       avail3[i][j].zmax = ( z_flag==1 ? z_cut : zmax  );
	     }

	     // Calculate volumes
	     for (j=0;j<7;j++) avail3[i][j].vol = (avail3[i][j].xmax-avail3[i][j].xmin) *
                                        	 (avail3[i][j].ymax-avail3[i][j].ymin) *
                                	         (avail3[i][j].zmax-avail3[i][j].zmin);
	   } // end if dim==1
	 } // end if (n < 2*ngl2 + 1)
        } // end loop on n for level 3 grids.
   } // end if NESTED_1

   PAMR_compose_hierarchy(1,3,1+ngl2+ngl3,lev_list,bbox_list,0.0);

   // ******************** Finished with initial grid hierarchy ******************* // 

   if (EXCISE_TEST) PAMR_save_gfn("chr",PAMR_AMRH,1,-1,-1.0,pre_tag,"_after_compose");

   // All tests will deal with these two grid functions.
   h1=PAMR_get_gfn(f1_name,PAMR_AMRH,1);  
   h2=PAMR_get_gfn(f2_name,PAMR_AMRH,1);  

   // for sync. tests ... set the gfs to something
   // interesting, but only *inside* ghost regions.

   if (SYNC_TEST)
   {
      if (rank==0) printf("SYNC_TEST start ... \n");

      for (l=1; l<4; l++)
      {
         valid=PAMR_init_s_iter(l,PAMR_AMRH,0);
         while(valid)
         {
            PAMR_get_g_attribs(&rank,&dim,shape,shape_c,bbox_b,ghost_width6,&t,&ngfs,x,x_c,gfs);
            sz=1;
	    // Do the following for vertex centered 
            for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape[i]=1;
    
            for (i=0; i<shape[0]; i++)
            {
               for (j=0; j<shape[1]; j++)
               {
                  for (k=0; k<shape[2]; k++)
                  {
		    
		    ind=i+j*shape[0]+k*shape[0]*shape[1];
		    x0=x[0][i];
		    if (dim>1) y0=x[1][j]; else y0=0;
		    if (dim>2) z0=x[2][k]; else z0=0;

		    // Determine whether this point is in a ghost region
		    amghost = 0;
		    if (ghost_width6[0] && (i<ghost_width6[0]))            amghost = 1;
		    if (ghost_width6[2] && (j<ghost_width6[2]) && (dim>1)) amghost = 1;
		    if (ghost_width6[4] && (k<ghost_width6[4]) && (dim>2)) amghost = 1;
		    if (ghost_width6[1] && (i>=(shape[0]-ghost_width6[1])))            amghost = 1;
		    if (ghost_width6[3] && (j>=(shape[1]-ghost_width6[3])) && (dim>1)) amghost = 1;
		    if (ghost_width6[5] && (k>=(shape[2]-ghost_width6[5])) && (dim>2)) amghost = 1;

		    if (l!=1 && amghost) (gfs[h1-1])[ind]=user_func(x0,y0,z0);
		    else (gfs[h1-1])[ind]=0;
                  }
               }
            }
	    
	    // Do the following for cell centered
            for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape_c[i]=1;
    
            for (i=0; i<shape_c[0]; i++)
            {
               for (j=0; j<shape_c[1]; j++)
               {
                  for (k=0; k<shape_c[2]; k++)
                  {

		    ind=i+j*shape_c[0]+k*shape_c[0]*shape_c[1];
		    x0=x_c[0][i];
		    if (dim>1) y0=x_c[1][j]; else y0=0;
		    if (dim>2) z0=x_c[2][k]; else z0=0;
		    
		    // Determine whether this point is in a ghost region
		    amghost = 0;
		    if (ghost_width6[0] && (i<ghost_width6[0]))            amghost = 1;
		    if (ghost_width6[2] && (j<ghost_width6[2]) && (dim>1)) amghost = 1;
		    if (ghost_width6[4] && (k<ghost_width6[4]) && (dim>2)) amghost = 1;
		    if (ghost_width6[1] && (i>=(shape_c[0]-ghost_width6[1])))            amghost = 1;
		    if (ghost_width6[3] && (j>=(shape_c[1]-ghost_width6[3])) && (dim>1)) amghost = 1;
		    if (ghost_width6[5] && (k>=(shape_c[2]-ghost_width6[5])) && (dim>2)) amghost = 1;
		    
		    if (l!=1 && amghost) (gfs[h2-1])[ind]=user_func(x0,y0,z0);
		    else (gfs[h2-1])[ind]=0;
                  }
               }
            }

            valid=PAMR_next_g();
         }
      }

      if (SAVE_SYNC_TEST) PAMR_save_gfn(f1_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_before_sync");
      if (SAVE_SYNC_TEST) PAMR_save_gfn(f2_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_before_sync");

      PAMR_sync(1,1,PAMR_AMRH,0);
      PAMR_sync(2,1,PAMR_AMRH,0);
      PAMR_sync(3,1,PAMR_AMRH,0);

      if (SAVE_SYNC_TEST) PAMR_save_gfn(f1_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_after_sync");
      if (SAVE_SYNC_TEST) PAMR_save_gfn(f2_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_after_sync");

      // Calculate errors
      v_norm = 0.0;
      v_counter = 0;
      v_norml = 0.0;
      v_counterl = 0;
      c_norm = 0.0;
      c_counter = 0;
      c_norml = 0.0;
      c_counterl = 0;

      for (l=1; l<4; l++)
      {
         valid=PAMR_init_s_iter(l,PAMR_AMRH,0);
         while(valid)
         {
            PAMR_get_g_attribs(&rank,&dim,shape,shape_c,bbox_b,ghost_width6,&t,&ngfs,x,x_c,gfs);
            sz=1;
            for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape[i]=1;
            for (i=0; i<shape[0]; i++)
            {
               for (j=0; j<shape[1]; j++)
               {
                  for (k=0; k<shape[2]; k++)
                  {
                     ind=i+j*shape[0]+k*shape[0]*shape[1];
		     v_norml += ((gfs[h1-1])[ind] * (gfs[h1-1])[ind]);
		     v_counterl++;
		  }
               }
            }
	    
            for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape_c[i]=1;
            for (i=0; i<shape_c[0]; i++)
            {
               for (j=0; j<shape_c[1]; j++)
               {
                  for (k=0; k<shape_c[2]; k++)
                  {
                     ind=i+j*shape_c[0]+k*shape_c[0]*shape_c[1];
		     c_norml += ((gfs[h2-1])[ind] * (gfs[h2-1])[ind]);
		     c_counterl++;
		  }
               }
            }

            valid=PAMR_next_g();
         }
      }

      MPI_Allreduce(&c_norml,&c_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&c_counterl,&c_counter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&v_norml,&v_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&v_counterl,&v_counter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
      if (rank == 0) {
	printf("l = %d, c_error_norm = %20.15e\n", l, sqrt(c_norm/c_counter));
	printf("l = %d, v_error_norm = %20.15e\n", l, sqrt(v_norm/v_counter));
      }

      if (rank==0) printf("... SYNC_TEST finished\n");
   }

   // for interpolation test ... set up with cubic polynomial in level 1
   // (4th order should be exact)

   if (INTERP_TEST)
   {
      if (rank==0) printf("INTERP_TEST start\n");

      for (l=1; l<4; l++)
      {
         valid=PAMR_init_s_iter(l,PAMR_AMRH,0);
         while(valid)
         {
            PAMR_get_g_attribs(&rank,&dim,shape,shape_c,bbox_b,ghost_width6,&t,&ngfs,x,x_c,gfs);
            sz=1;
	    // Do the following for vertex centered 
            for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape[i]=1;
    
            for (i=0; i<shape[0]; i++)
            {
               for (j=0; j<shape[1]; j++)
               {
                  for (k=0; k<shape[2]; k++)
                  {
		    ind=i+j*shape[0]+k*shape[0]*shape[1];
		    x0=x[0][i];
		    if (dim>1) y0=x[1][j]; else y0=0;
		    if (dim>2) z0=x[2][k]; else z0=0;
		    if (l==1) (gfs[h1-1])[ind]=user_func(x0,y0,z0);
		    else (gfs[h1-1])[ind]=0;
                  }
               }
            }
	    
	    // Do the following for cell centered
            for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape_c[i]=1;
    
            for (i=0; i<shape_c[0]; i++)
            {
               for (j=0; j<shape_c[1]; j++)
               {
                  for (k=0; k<shape_c[2]; k++)
                  {
                     ind=i+j*shape_c[0]+k*shape_c[0]*shape_c[1];
                     x0=x_c[0][i];
                     if (dim>1) y0=x_c[1][j]; else y0=0;
                     if (dim>2) z0=x_c[2][k]; else z0=0;
		     if (l==1) (gfs[h2-1])[ind]=user_func(x0,y0,z0);
		     else (gfs[h2-1])[ind]=0;
                  }
               }
            }

            valid=PAMR_next_g();
         }
      }

      if (SAVE_INTERP_TEST) PAMR_save_gfn(f1_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_before_interp");
      if (SAVE_INTERP_TEST) PAMR_save_gfn(f2_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_before_interp");

      PAMR_interp(1,1,PAMR_AMRH);
      PAMR_interp(2,1,PAMR_AMRH); 

      if (SAVE_INTERP_TEST) PAMR_save_gfn(f1_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_after_interp");
      if (SAVE_INTERP_TEST) PAMR_save_gfn(f2_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_after_interp");

      // Next calculate the errors.

      c_counter = 0;
      c_norm = 0.0;
      c_counterl = 0;
      c_norml = 0.0;
      v_counter = 0;
      v_norm = 0.0;
      v_counterl = 0;
      v_norml = 0.0;
    
      for (l=1; l<4; l++)
	{
	  valid=PAMR_init_s_iter(l,PAMR_AMRH,0);
	  while(valid)
	    {
	      PAMR_get_g_attribs(&rank,&dim,shape,shape_c,bbox_b,ghost_width6,&t,&ngfs,x,x_c,gfs);
	      sz=1;
	      // Do the following for vertex centered
	      ilo = 0;
	      jlo = 0;
	      klo = 0;
	      
	      if (ghost_width6[0]) ilo = ghost_width6[0] + 1;
	      iup = shape[0] - ghost_width6[1];
	      if (dim > 1) {
		if (ghost_width6[2]) jlo = ghost_width6[2] + 1;
		jup = shape[1] - ghost_width6[3];
	      } else {
		jup = shape[1];
	      }
	      if (dim > 2) {
		if (ghost_width6[4]) klo = ghost_width6[4] + 1;
		kup = shape[2] - ghost_width6[5];
	      } else {
		kup = shape[2];
	      }
	      
	      for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape[i]=1;
	      for (i=0; i<shape[0]; i++)
		{
		  for (j=0; j<shape[1]; j++)
		    {
		      for (k=0; k<shape[2]; k++)
			{
			  
			  ind=i+j*shape[0]+k*shape[0]*shape[1];
			  x0=x[0][i];
			  if (dim>1) y0=x[1][j]; else y0=0;
			  if (dim>2) z0=x[2][k]; else z0=0;
			  (gfs[h1-1])[ind] -= user_func(x0,y0,z0);
			  (gfs[h1-1])[ind] /= user_func(x0,y0,z0);  // Calculate relative error.
			  // Warning message
			  if ( fabs((gfs[h1-1])[ind]) > 0.1) {
			    printf("Warning, large errors at i=%d, j=%d, k=%d on proc %d\n",i,j,k,rank);
			  }

			  // Notice that this will give you problems if your function passes through zero anywhere.
			  if ( (i>=ilo && i<iup) && (j>=jlo && j<jup) && (k>=klo && k<kup)) {		      
			    v_norml += ((gfs[h1-1])[ind] * (gfs[h1-1])[ind]);
			    v_counterl++;
			  }
			}
		    }
		}
	      
	      ilo = ghost_width6[0];
	      iup = shape_c[0] - ghost_width6[1];
	      if (dim > 1) {
		jlo = ghost_width6[2];
		jup = shape_c[1] - ghost_width6[3];
	      } else {
		jlo = 0;
		jup = shape_c[1];
	      }
	      if (dim > 2) {
		klo = ghost_width6[4];
		kup = shape_c[2] - ghost_width6[5];
	      } else {
		klo = 0;
		kup = shape_c[2];
	      }
	      
	      // Do the following for cell centered
	      for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape_c[i]=1;
	      for (i=0; i<shape_c[0]; i++)
		{
		  for (j=0; j<shape_c[1]; j++)
		    {
		      for (k=0; k<shape_c[2]; k++)
			{
			  ind=i+j*shape_c[0]+k*shape_c[0]*shape_c[1];
			  x0=x_c[0][i];
			  if (dim>1) y0=x_c[1][j]; else y0=0;
			  if (dim>2) z0=x_c[2][k]; else z0=0;
			  (gfs[h2-1])[ind] -= user_func(x0,y0,z0); 
			  (gfs[h2-1])[ind] /= user_func(x0,y0,z0); 
			  if ( fabs((gfs[h2-1])[ind]) > 0.1) {
			    printf("Warning, large errors at i=%d, j=%d, k=%d on proc %d\n",i,j,k,rank);
			  }

			  if ((i>=ilo && i<iup) && (j>=jlo && j<jup) && (k>=klo && k<kup)) {
			    c_norml += ((gfs[h2-1])[ind] * (gfs[h2-1])[ind]);
			    c_counterl++;
			  }
			}
		    }
		}
	      valid=PAMR_next_g();
	    }
	}
      
      MPI_Allreduce(&v_norml,&v_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&v_counterl,&v_counter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&c_norml,&c_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&c_counterl,&c_counter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
      if (rank==0) {
	printf("l = %d, v_error_norm = %20.15e\n", l, sqrt(v_norm/v_counter));
	printf("l = %d, c_error_norm = %20.15e\n", l, sqrt(c_norm/c_counter));
      }

      if (SAVE_INTERP_TEST) PAMR_save_gfn(f1_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_error");
      if (SAVE_INTERP_TEST) PAMR_save_gfn(f2_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_error");

      if (rank==0) printf("... INTERP_TEST finished\n");
   }

   // boundary interpolation test

   if (BDY_INTERP_TEST)
   {
      if (rank==0) printf("BDY_INTERP_TEST start\n");

      wamrbdy_v = 2;
      wamrbdy_c = 1;

      int interior;
      int interiorx, interiory, interiorz;

      for (l=1; l<4; l++)
      {
         valid=PAMR_init_s_iter(l,PAMR_AMRH,0);
         while(valid)
         {
            PAMR_get_g_attribs(&rank,&dim,shape,shape_c,bbox_b,ghost_width6,&t,&ngfs,x,x_c,gfs);
            sz=1;
            for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape[i]=1;
            for (i=0; i<shape[0]; i++)
            {
               for (j=0; j<shape[1]; j++)
               {
                  for (k=0; k<shape[2]; k++)
                  {
		    interior = 0;
		    interiorx = 0;
		    interiory = 0;
		    interiorz = 0;
		    if ((i>=wamrbdy_v) && i<(shape[0]-wamrbdy_v)) interiorx = 1;
		    if ((j>=wamrbdy_v) && j<(shape[1]-wamrbdy_v)) interiory = 1;
		    if ((k>=wamrbdy_v) && k<(shape[2]-wamrbdy_v)) interiorz = 1;

		    // Now, if it's at a processor boundary, let's just call it interior.
		    if (ghost_width6[0] && (i<=ghost_width6[0])) interiorx = 1;
		    if (ghost_width6[1] && (i>=(shape[0]-ghost_width6[1]-1))) interiorx = 1;
		    if (ghost_width6[2] && (j<=ghost_width6[2])) interiory = 1;
		    if (ghost_width6[3] && (j>=(shape[1]-ghost_width6[3]-1))) interiory = 1;
		    if (ghost_width6[4] && (k<=ghost_width6[4])) interiorz = 1;
		    if (ghost_width6[5] && (k>=(shape[2]-ghost_width6[5]-1))) interiorz = 1;

		    interior = interiorx;
		    if (dim > 1) interior = interior * interiory;
		    if (dim > 2) interior = interior * interiorz;
		    
		    ind=i+j*shape[0]+k*shape[0]*shape[1];
		    x0=x[0][i];
		    if (dim>1) y0=x[1][j]; else y0=0;
		    if (dim>2) z0=x[2][k]; else z0=0;
		    
		    if (interior || l==1) {
		      (gfs[h1-1])[ind]=user_func(x0,y0,z0);
		    } else {
		      (gfs[h1-1])[ind]=-100.0; 
		    }
                  }
               }
            }

            for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape_c[i]=1;
            for (i=0; i<shape_c[0]; i++)
            {
               for (j=0; j<shape_c[1]; j++)
               {
                  for (k=0; k<shape_c[2]; k++)
                  {

		    interior = 0;
		    interiorx = 0;
		    interiory = 0;
		    interiorz = 0;
		    if ((i>=wamrbdy_c) && i<(shape_c[0]-wamrbdy_c)) interiorx = 1;
		    if ((j>=wamrbdy_c) && j<(shape_c[1]-wamrbdy_c)) interiory = 1;
		    if ((k>=wamrbdy_c) && k<(shape_c[2]-wamrbdy_c)) interiorz = 1;

		    // Now, if it's at a processor boundary, let's just call it interior.
		    if (ghost_width6[0] && (i<ghost_width6[0])) interiorx = 1;
		    if (ghost_width6[1] && (i>=(shape_c[0]-ghost_width6[1]))) interiorx = 1;
		    if (ghost_width6[2] && (j<ghost_width6[2])) interiory = 1;
		    if (ghost_width6[3] && (j>=(shape_c[1]-ghost_width6[3]))) interiory = 1;
		    if (ghost_width6[4] && (k<ghost_width6[4])) interiorz = 1;
		    if (ghost_width6[5] && (k>=(shape_c[2]-ghost_width6[5]))) interiorz = 1;

		    interior = interiorx;
		    if (dim > 1) interior = interior * interiory;
		    if (dim > 2) interior = interior * interiorz;

                     ind=i+j*shape_c[0]+k*shape_c[0]*shape_c[1];
                     x0=x_c[0][i];
                     if (dim>1) y0=x_c[1][j]; else y0=0;
                     if (dim>2) z0=x_c[2][k]; else z0=0;
		     if (interior || l==1) {
		       (gfs[h2-1])[ind]=user_func(x0,y0,z0);
		     } else {
		       (gfs[h2-1])[ind]=-100.0; 
		     }
                  }
               }
            }

            valid=PAMR_next_g();
         }
      }

      // Sync first to fix crap at processor boundaries.
      PAMR_sync(2,1,PAMR_AMRH,wamrbdy_c);
      PAMR_sync(3,1,PAMR_AMRH,wamrbdy_c);

      if (SAVE_BDY_INTERP_TEST) PAMR_save_gfn(f1_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_before_bdy_interp");
      if (SAVE_BDY_INTERP_TEST) PAMR_save_gfn(f2_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_before_bdy_interp");

      PAMR_AMR_bdy_interp(1,1,wamrbdy_v);
      PAMR_AMR_bdy_interp_c(1,1,wamrbdy_c);

      PAMR_AMR_bdy_interp(2,1,wamrbdy_v);
      PAMR_AMR_bdy_interp_c(2,1,wamrbdy_c);

      PAMR_sync(2,1,PAMR_AMRH,wamrbdy_c);
      PAMR_sync(3,1,PAMR_AMRH,wamrbdy_c);

      if (SAVE_BDY_INTERP_TEST) PAMR_save_gfn(f1_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_after_bdy_interp");
      if (SAVE_BDY_INTERP_TEST) PAMR_save_gfn(f2_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_after_bdy_interp");

      c_norm = 0.0;
      c_counter = 0;
      c_norml = 0.0;
      c_counterl = 0;
      v_norm = 0.0;
      v_counter = 0;
      v_norml = 0.0;
      v_counterl = 0;

      for (l=1; l<4; l++)
      {
         valid=PAMR_init_s_iter(l,PAMR_AMRH,0);
	 while(valid)
	   {
	     PAMR_get_g_attribs(&rank,&dim,shape,shape_c,bbox_b,ghost_width6,&t,&ngfs,x,x_c,gfs);
	     sz=1;
	  
	     for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape_c[i]=1;
	     for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape[i]=1;

	     ilo = ghost_width6[0];
	     iup = shape_c[0] - ghost_width6[1];
	     if (dim > 1) {
	       jlo = ghost_width6[2];
	       jup = shape_c[1] - ghost_width6[3];
	     } else {
	       jlo = 0;
	       jup = shape_c[1];
	     }
	     if (dim > 2) {
	       klo = ghost_width6[4];
	       kup = shape_c[2] - ghost_width6[5];
	     } else {
	       klo = 0;
	       kup = shape_c[2];
	     }
	     
	     for (i=0; i<shape_c[0]; i++)
	       {
		 for (j=0; j<shape_c[1]; j++)
		   {
		     for (k=0; k<shape_c[2]; k++)
		       {

			 interior = 0;
			 interiorx = 0;
			 interiory = 0;
			 interiorz = 0;
			 if ((i>=wamrbdy_c) && i<(shape_c[0]-wamrbdy_c)) interiorx = 1;
			 if ((j>=wamrbdy_c) && j<(shape_c[1]-wamrbdy_c)) interiory = 1;
			 if ((k>=wamrbdy_c) && k<(shape_c[2]-wamrbdy_c)) interiorz = 1;
			 
			 interior = interiorx;
			 if (dim > 1) interior = interior * interiory;
			 if (dim > 2) interior = interior * interiorz;
			 
			 ind=i+j*shape_c[0]+k*shape_c[0]*shape_c[1];
			 x0=x_c[0][i];
			 if (dim>1) y0=x_c[1][j]; else y0=0;
			 if (dim>2) z0=x_c[2][k]; else z0=0;
			 (gfs[h2-1])[ind] -= user_func(x0,y0,z0);
			 (gfs[h2-1])[ind] /= user_func(x0,y0,z0);
			  if ( fabs((gfs[h2-1])[ind]) > 0.1) {
			    printf("Warning, large errors at i=%d, j=%d, k=%d on proc %d\n",i,j,k,rank);
			  }

			 if ((i>=ilo && i<iup) && (j>=jlo && j<jup) && (k>=klo && k<kup)) {
			   c_norml += (gfs[h2-1])[ind]*(gfs[h2-1])[ind];
			   c_counterl++;
			 }
			 
		       }
		   }
	       }

	     ilo = ghost_width6[0];
	     iup = shape[0] - ghost_width6[1];
	     if (ghost_width6[1] > 0) --iup;
	     if (dim > 1) {
	       jlo = ghost_width6[2];
	       jup = shape[1] - ghost_width6[3];
	       if (ghost_width6[3] > 0) --jup;
	     } else {
	       jlo = 0;
	       jup = shape[1];
	     }
	     if (dim > 2) {
	       klo = ghost_width6[4];
	       kup = shape[2] - ghost_width6[5];
	       if (ghost_width6[5] > 0) --kup;
	     } else {
	       klo = 0;
	       kup = shape[2];
	     }
     
	     for (i=0; i<shape[0]; i++)
	       {
		 for (j=0; j<shape[1]; j++)
		   {
		     for (k=0; k<shape[2]; k++)
		       {

			 interior = 0;
			 interiorx = 0;
			 interiory = 0;
			 interiorz = 0;

			 if ((i>=wamrbdy_v) && i<(shape[0]-wamrbdy_v)) interiorx = 1;
			 if ((j>=wamrbdy_v) && j<(shape[1]-wamrbdy_v)) interiory = 1;
			 if ((k>=wamrbdy_v) && k<(shape[2]-wamrbdy_v)) interiorz = 1;

			 interior = interiorx;
			 if (dim > 1) interior = interior * interiory;
			 if (dim > 2) interior = interior * interiorz;

			 ind=i+j*shape[0]+k*shape[0]*shape[1];
			 x0=x[0][i];
			 if (dim>1) y0=x[1][j]; else y0=0;
			 if (dim>2) z0=x[2][k]; else z0=0;
			 
			 (gfs[h1-1])[ind] -= user_func(x0,y0,z0);
			 (gfs[h1-1])[ind] /= user_func(x0,y0,z0);
			  if ( fabs((gfs[h1-1])[ind]) > 0.1) {
			    printf("Warning, large errors at i=%d, j=%d, k=%d on proc %d\n",i,j,k,rank);
			  }
			 
			 if ((i>=ilo && i<iup) && (j>=jlo && j<jup) && (k>=klo && k<kup)) {
			   v_norml += (gfs[h1-1])[ind]*(gfs[h1-1])[ind];
			   v_counterl++;
			 }
		       }
		   }
	       }
	     valid=PAMR_next_g();
	   }
      }      

      MPI_Allreduce(&v_norml,&v_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&v_counterl,&v_counter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&c_norml,&c_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&c_counterl,&c_counter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
      if (rank == 0) {
	printf("v_error_norm = %20.15e\n", sqrt(v_norm/v_counter));
	printf("c_error_norm = %20.15e\n", sqrt(c_norm/c_counter));
      }

      if (SAVE_BDY_INTERP_TEST) PAMR_save_gfn(f1_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_error");
      if (SAVE_BDY_INTERP_TEST) PAMR_save_gfn(f2_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_error");
      
      if (rank==0) printf("... BDY_INTERP_TEST finished\n");
   }

   if (INJECT_TEST)
   {
      if (rank==0) printf("INJECT_TEST start\n");


      for (l=1; l<4; l++)
      {
         valid=PAMR_init_s_iter(l,PAMR_AMRH,0);
         while(valid)
         {
            PAMR_get_g_attribs(&rank,&dim,shape,shape_c,bbox_b,ghost_width6,&t,&ngfs,x,x_c,gfs);

            sz=1;
            for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape[i]=1;
            for (i=0; i<shape[0]; i++)
            {
               for (j=0; j<shape[1]; j++)
               {
                  for (k=0; k<shape[2]; k++)
                  {
                     ind=i+j*shape[0]+k*shape[0]*shape[1];
                     x0=x[0][i];
                     if (dim>1) y0=x[1][j]; else y0=0;
                     if (dim>2) z0=x[2][k]; else z0=0;

		     if (l<3) {
		       // Determine whether this point is in the shadow of a level 3 grid.
		       for (n=ngl2+1; n<(ngl2+ngl3+1); n++) {
			 in_shadow = 0;
			 isx = 0;
			 isy = 0;
			 isz = 0;

			 if ((x0 <= bbox_list[2*dim*n+1]) && (x0 >= bbox_list[2*dim*n  ])) isx = 1;
			 if ((y0 <= bbox_list[2*dim*n+3]) && (y0 >= bbox_list[2*dim*n+2])) isy = 1;
			 if ((z0 <= bbox_list[2*dim*n+5]) && (z0 >= bbox_list[2*dim*n+4])) isz = 1;
			 
			 in_shadow = isx;
			 if (dim > 1) in_shadow *= isy;
			 if (dim > 2) in_shadow *= isz;
			 if (in_shadow) break;
		       }
		     }

		     if (l==3 || !in_shadow) {
		       (gfs[h1-1])[ind]= user_func(x0,y0,z0);
		     } else {
		       (gfs[h1-1])[ind]=0.0;
		     }
                  }
               }
            }

            for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape_c[i]=1;
            for (i=0; i<shape_c[0]; i++)
            {
               for (j=0; j<shape_c[1]; j++)
               {
                  for (k=0; k<shape_c[2]; k++)
                  {
                     ind=i+j*shape_c[0]+k*shape_c[0]*shape_c[1];
                     x0=x_c[0][i];
                     if (dim>1) y0=x_c[1][j]; else y0=0;
                     if (dim>2) z0=x_c[2][k]; else z0=0;

		     if (l<3) {
		       // Determine whether this point is in the shadow of a level 3 grid.
		       for (n=ngl2+1; n<(ngl2+ngl3+1); n++) {
			 in_shadow = 0;
			 isx = 0;
			 isy = 0;
			 isz = 0;
			 if ((x0 <= bbox_list[2*dim*n+1]) && (x0 >= bbox_list[2*dim*n  ])) isx = 1;
			 if ((y0 <= bbox_list[2*dim*n+3]) && (y0 >= bbox_list[2*dim*n+2])) isy = 1;
			 if ((z0 <= bbox_list[2*dim*n+5]) && (z0 >= bbox_list[2*dim*n+4])) isz = 1;
			 
			 in_shadow = isx;
			 if (dim > 1) in_shadow *= isy;
			 if (dim > 2) in_shadow *= isz;
			 if (in_shadow) break;
		       }
		     }

		     if (l==3 || !in_shadow) {
		       (gfs[h2-1])[ind]=user_func(x0,y0,z0);
		     } else {
		       (gfs[h2-1])[ind]=0;
		     }
		  }
               }
            }

            valid=PAMR_next_g();
         }
      }

      if (SAVE_INJECT_TEST) PAMR_save_gfn(f1_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_before_inject");
      if (SAVE_INJECT_TEST) PAMR_save_gfn(f2_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_before_inject"); 

      PAMR_inject(3,1,PAMR_AMRH);
      PAMR_inject(2,1,PAMR_AMRH);

      if (SAVE_INJECT_TEST) PAMR_save_gfn(f1_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_after_inject");
      if (SAVE_INJECT_TEST) PAMR_save_gfn(f2_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_after_inject"); 

      v_norm = 0.0;
      v_counter = 0;
      v_norml = 0.0;
      v_counterl = 0;
      c_norm = 0.0;
      c_counter = 0;
      c_norml = 0.0;
      c_counterl = 0;

      for (l=1; l<4; l++)
      {
         valid=PAMR_init_s_iter(l,PAMR_AMRH,0);
         while(valid)
         {
            PAMR_get_g_attribs(&rank,&dim,shape,shape_c,bbox_b,ghost_width6,&t,&ngfs,x,x_c,gfs);
	    
            sz=1;

	    // Try excluding the ghost zones *and* vertices which 
	    // are *shared* by more than one processor.
	    ilo = 0;
	    jlo = 0;
	    klo = 0;
	    
	    if (ghost_width6[0]) ilo = ghost_width6[0] + 1;
	    iup = shape[0] - ghost_width6[1];
	    if (dim > 1) {
	      if (ghost_width6[2]) jlo = ghost_width6[2] + 1;
	      jup = shape[1] - ghost_width6[3];
	    } else {
	      jup = shape[1];
	    }
	    if (dim > 2) {
	      if (ghost_width6[4]) klo = ghost_width6[4] + 1;
	      kup = shape[2] - ghost_width6[5];
	    } else {
	      kup = shape[2];
	    }
	    
            for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape[i]=1;
            for (i=0; i<shape[0]; i++)
            {
               for (j=0; j<shape[1]; j++)
               {
                  for (k=0; k<shape[2]; k++)
                  {
                     ind=i+j*shape[0]+k*shape[0]*shape[1];
                     x0=x[0][i];
                     if (dim>1) y0=x[1][j]; else y0=0;
                     if (dim>2) z0=x[2][k]; else z0=0;
		     (gfs[h1-1])[ind] -= user_func(x0,y0,z0);
		     (gfs[h1-1])[ind] /= user_func(x0,y0,z0);
		     if ( fabs((gfs[h1-1])[ind]) > 0.1) {
		       printf("Warning, large errors at i=%d, j=%d, k=%d on proc %d\n",i,j,k,rank);
		     }

		     if ( (i>=ilo && i<iup) && (j>=jlo && j<jup) && (k>=klo && k<kup)) {		      		       
		       v_norml += ((gfs[h1-1])[ind] * (gfs[h1-1])[ind]);
		       v_counterl++;
		     }
		  }
               }
            }

	    // Exclude the ghost zones.
	    ilo = ghost_width6[0];
	    iup = shape_c[0] - ghost_width6[1];
	    if (dim > 1) {
	      jlo = ghost_width6[2];
	      jup = shape_c[1] - ghost_width6[3];
	    } else {
	      jlo = 0;
	      jup = shape_c[1];
	    }
	    if (dim > 2) {
	      klo = ghost_width6[4];
	      kup = shape_c[2] - ghost_width6[5];
	    } else {
	      klo = 0;
	      kup = shape_c[2];
	    }
	    
            for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape_c[i]=1;
            for (i=0; i<shape_c[0]; i++)
            {
               for (j=0; j<shape_c[1]; j++)
               {
                  for (k=0; k<shape_c[2]; k++)
                  {
                     ind=i+j*shape_c[0]+k*shape_c[0]*shape_c[1];
                     x0=x_c[0][i];
                     if (dim>1) y0=x_c[1][j]; else y0=0;
                     if (dim>2) z0=x_c[2][k]; else z0=0;
		     (gfs[h2-1])[ind] -= user_func(x0,y0,z0);
		     (gfs[h2-1])[ind] /= user_func(x0,y0,z0);
		     if ( fabs((gfs[h2-1])[ind]) > 0.1) {
		       printf("Warning, large errors at i=%d, j=%d, k=%d on proc %d\n",i,j,k,rank);
		     }

		     if ( (i>=ilo && i<iup) && (j>=jlo && j<jup) && (k>=klo && k<kup)) {		      		       
		       c_norml += ((gfs[h2-1])[ind] * (gfs[h2-1])[ind]);
		       c_counterl++;
		     }
		  }
               }
            }


            valid=PAMR_next_g();
         }
      }
      
      MPI_Allreduce(&c_norml,&c_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&c_counterl,&c_counter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&v_norml,&v_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&v_counterl,&v_counter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
      if (rank == 0) {
	printf("l = %d, c_error_norm = %20.15e\n", l, sqrt(c_norm/c_counter));
	printf("l = %d, v_error_norm = %20.15e\n", l, sqrt(v_norm/v_counter));
      }
      
      if (SAVE_INJECT_TEST) PAMR_save_gfn(f1_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_error");
      if (SAVE_INJECT_TEST) PAMR_save_gfn(f2_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_error"); 
      
      if (rank==0) printf("... INJECT_TEST finished\n");
   }

   if (C_TO_V_TEST)
    {
      if (rank==0) printf("C_TO_V_TEST start\n");

      for (l=1; l<4; l++)
      {
         valid=PAMR_init_s_iter(l,PAMR_AMRH,0);
         while(valid)
         {
            PAMR_get_g_attribs(&rank,&dim,shape,shape_c,bbox_b,ghost_width6,&t,&ngfs,x,x_c,gfs);
            sz=1;

	    // Do the following for vertex centered 
            for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape[i]=1;
    
            for (i=0; i<shape[0]; i++)
            {
               for (j=0; j<shape[1]; j++)
               {
                  for (k=0; k<shape[2]; k++)
                  {
		    ind=i+j*shape[0]+k*shape[0]*shape[1];
		    (gfs[h1-1])[ind]=0;
                  }
               }
            }
	    
	    // Do the following for cell centered
            for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape_c[i]=1;
    
            for (i=0; i<shape_c[0]; i++)
            {
               for (j=0; j<shape_c[1]; j++)
               {
                  for (k=0; k<shape_c[2]; k++)
                  {
		    ind=i+j*shape_c[0]+k*shape_c[0]*shape_c[1];
		    x0=x_c[0][i];
		    if (dim>1) y0=x_c[1][j]; else y0=0;
		    if (dim>2) z0=x_c[2][k]; else z0=0;
		    (gfs[h2-1])[ind]=user_func(x0,y0,z0);
		  }
	       }
            }
	    
            valid=PAMR_next_g();
	 }
      }
      
      if (SAVE_C_TO_V_TEST) PAMR_save_gfn(f1_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_before_c_to_v");
      if (SAVE_C_TO_V_TEST) PAMR_save_gfn(f2_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_before_c_to_v");

      PAMR_c_to_v(1,1,PAMR_AMRH,0);
      PAMR_c_to_v(2,1,PAMR_AMRH,0);
      PAMR_c_to_v(3,1,PAMR_AMRH,0);

      if (SAVE_C_TO_V_TEST) PAMR_save_gfn(f1_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_after_c_to_v");

      // Now calculate error.
      v_norm = 0.0;
      v_counter = 0;
      v_norml = 0.0;
      v_counterl = 0;

      for (l=1; l<4; l++)
      {

	valid=PAMR_init_s_iter(l,PAMR_AMRH,0);
	while(valid)
	  {
	    PAMR_get_g_attribs(&rank,&dim,shape,shape_c,bbox_b,ghost_width6,&t,&ngfs,x,x_c,gfs);
	    sz=1;
	  
	    // Try excluding the ghost zones *and* vertices which 
	    // are *shared* by more than one processor.
	    ilo = 0;
	    jlo = 0;
	    klo = 0;
	  
	    if (ghost_width6[0]) ilo = ghost_width6[0] + 1;
	    iup = shape[0] - ghost_width6[1];
	    if (dim > 1) {
	      if (ghost_width6[2]) jlo = ghost_width6[2] + 1;
	      jup = shape[1] - ghost_width6[3];
	    } else {
	      jup = shape[1];
	    }
	    if (dim > 2) {
	      if (ghost_width6[4]) klo = ghost_width6[4] + 1;
	      kup = shape[2] - ghost_width6[5];
	    } else {
	      kup = shape[2];
	    }
	    
	    for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape[i]=1;
	    
	    for (i=0; i<shape[0]; i++)
	      {
		for (j=0; j<shape[1]; j++)
		  {
		    for (k=0; k<shape[2]; k++)
		      {
			ind=i+j*shape[0]+k*shape[0]*shape[1];
			x0=x[0][i];
			if (dim>1) y0=x[1][j]; else y0=0;
			if (dim>2) z0=x[2][k]; else z0=0;
			(gfs[h1-1])[ind] -= user_func(x0,y0,z0);
			(gfs[h1-1])[ind] /= user_func(x0,y0,z0);
			if ( fabs((gfs[h1-1])[ind]) > 0.1) {
			  printf("Warning, large errors at i=%d, j=%d, k=%d on proc %d\n",i,j,k,rank);
			}

			if ( (i>=ilo && i<iup) && (j>=jlo && j<jup) && (k>=klo && k<kup)) {		      
			  v_norml += (gfs[h1-1])[ind]*(gfs[h1-1])[ind];
			  v_counterl++;
			}
		      }
		  }
	      }
	    
	    valid=PAMR_next_g();
	  }
      }

      MPI_Allreduce(&v_norml,&v_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&v_counterl,&v_counter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
      if (rank==0) printf("l = %d, v_error_norm = %20.15e\n", l, sqrt(v_norm/v_counter));

      if (SAVE_C_TO_V_TEST) PAMR_save_gfn(f1_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_error");
   }

   if (V_TO_C_TEST)
   {
      if (rank==0) printf("V_TO_C_TEST start\n");

      for (l=1; l<4; l++)
      {
         valid=PAMR_init_s_iter(l,PAMR_AMRH,0);
         while(valid)
         {
            PAMR_get_g_attribs(&rank,&dim,shape,shape_c,bbox_b,ghost_width6,&t,&ngfs,x,x_c,gfs);
            sz=1;
	    
	    // Do the following for vertex centered 
            for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape[i]=1;
    
            for (i=0; i<shape[0]; i++)
            {
               for (j=0; j<shape[1]; j++)
               {
                  for (k=0; k<shape[2]; k++)
                  {

		    ind=i+j*shape[0]+k*shape[0]*shape[1];
		    x0=x[0][i];
		    if (dim>1) y0=x[1][j]; else y0=0;
		    if (dim>2) z0=x[2][k]; else z0=0;
		    (gfs[h1-1])[ind]= user_func(x0,y0,z0);
		  }
               }
            }
	    
	    // Do the following for cell centered
            for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape_c[i]=1;
    
            for (i=0; i<shape_c[0]; i++)
            {
               for (j=0; j<shape_c[1]; j++)
               {
                  for (k=0; k<shape_c[2]; k++)
                  {
		    ind=i+j*shape_c[0]+k*shape_c[0]*shape_c[1];
		    (gfs[h2-1])[ind]=0.0;
                  }
               }
            }
	    
            valid=PAMR_next_g();
         }
      }

      if (SAVE_V_TO_C_TEST) PAMR_save_gfn(f1_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_before_v_to_c");
      if (SAVE_V_TO_C_TEST) PAMR_save_gfn(f2_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_before_v_to_c");

      PAMR_v_to_c(1,1,PAMR_AMRH,0);
      PAMR_v_to_c(2,1,PAMR_AMRH,0);
      PAMR_v_to_c(3,1,PAMR_AMRH,0);

      if (SAVE_V_TO_C_TEST) PAMR_save_gfn(f2_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_after_v_to_c");

      // Now calculate errors
      c_norm = 0.0;
      c_counter = 0;
      c_norml = 0.0;
      c_counterl = 0;

      for (l=1; l<4; l++)
      {
	valid=PAMR_init_s_iter(l,PAMR_AMRH,0);
	while(valid)
	  {
	    PAMR_get_g_attribs(&rank,&dim,shape,shape_c,bbox_b,ghost_width6,&t,&ngfs,x,x_c,gfs);
	    sz=1;
	    
	    for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape_c[i]=1;
	    
	    // Try excluding the ghost zones.
	    ilo = ghost_width6[0];
	    iup = shape_c[0] - ghost_width6[1];
	    if (dim > 1) {
	      jlo = ghost_width6[2];
	      jup = shape_c[1] - ghost_width6[3];
	    } else {
	      jlo = 0;
	      jup = shape_c[1];
	    }
	    if (dim > 2) {
	      klo = ghost_width6[4];
	      kup = shape_c[2] - ghost_width6[5];
	    } else {
	      klo = 0;
	      kup = shape_c[2];
	    }
	    
	    for (i=0; i<shape_c[0]; i++)
	      {
		for (j=0; j<shape_c[1]; j++)
		  {
		    for (k=0; k<shape_c[2]; k++)
		      {
			ind=i+j*shape_c[0]+k*shape_c[0]*shape_c[1];
			x0=x_c[0][i];
			if (dim>1) y0=x_c[1][j]; else y0=0;
			if (dim>2) z0=x_c[2][k]; else z0=0;
			(gfs[h2-1])[ind] -= user_func(x0,y0,z0); 
			(gfs[h2-1])[ind] /= user_func(x0,y0,z0); 
			if ( fabs((gfs[h2-1])[ind]) > 0.1) {
			  printf("Warning, large errors at i=%d, j=%d, k=%d on proc %d\n",i,j,k,rank);
			}

			if ( (i>=ilo && i<iup) && (j>=jlo && j<jup) && (k>=klo && k<kup)) {
			  c_norml += (gfs[h2-1])[ind]*(gfs[h2-1])[ind];
			  c_counterl++;
			}
		      
		      }
		  }
	      }
	  
	    valid=PAMR_next_g();
	  }
      }
      
      MPI_Allreduce(&c_norml,&c_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&c_counterl,&c_counter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
      if (rank == 0) printf("lev = %d, c_error_norm = %20.15e\n", l, sqrt(c_norm/c_counter));
      if (SAVE_V_TO_C_TEST) PAMR_save_gfn(f2_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_error");
   }


   if (MGH_TEST)
   {
     if (rank == 0) printf("MGH_TEST start\n");

      for (l=1; l<4; l++)
      {
         valid=PAMR_init_s_iter(l,PAMR_AMRH,0);
         while(valid)
         {
            PAMR_get_g_attribs(&rank,&dim,shape,shape_c,bbox_b,ghost_width6,&t,&ngfs,x,x_c,gfs);

            sz=1;
            for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape[i]=1;
            for (i=0; i<shape[0]; i++)
            {
               for (j=0; j<shape[1]; j++)
               {
                  for (k=0; k<shape[2]; k++)
                  {
                     ind=i+j*shape[0]+k*shape[0]*shape[1];
                     x0=x[0][i];
                     if (dim>1) y0=x[1][j]; else y0=0;
                     if (dim>2) z0=x[2][k]; else z0=0;
                     (gfs[h1-1])[ind]=(gfs[h4-1])[ind]=user_func(x0,y0,z0);
                  }
               }
            }


            for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape_c[i]=1;
            for (i=0; i<shape_c[0]; i++)
            {
               for (j=0; j<shape_c[1]; j++)
               {
                  for (k=0; k<shape_c[2]; k++)
                  {
                     ind=i+j*shape_c[0]+k*shape_c[0]*shape_c[1];
                     x0=x_c[0][i];
                     if (dim>1) y0=x_c[1][j]; else y0=0;
                     if (dim>2) z0=x_c[2][k]; else z0=0;
                     (gfs[h2-1])[ind]=user_func(x0,y0,z0);
                  }
               }
            }
            valid=PAMR_next_g();
         }
      }

      max_lev=PAMR_build_mgh(1,3,1);

      for (l=1; l<max_lev; l++)
      {
         valid=PAMR_init_s_iter(l,PAMR_MGH,0);
         while(valid)
         {
            PAMR_get_g_attribs(&rank,&dim,shape,shape_c,bbox_b,ghost_width6,&t,&ngfs,x,x_c,gfs);
            PAMR_get_g_coarsest(&coarsest);
            sz=1;

            for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape[i]=1;
            for (i=0; i<shape[0]; i++)
            {
               for (j=0; j<shape[1]; j++)
               {
                  for (k=0; k<shape[2]; k++)
                  {
                     ind=i+j*shape[0]+k*shape[0]*shape[1];
                     (gfs[h5-1])[ind]=coarsest;
                  }
               }
            }

            for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape_c[i]=1;
            for (i=0; i<shape_c[0]; i++)
            {
               for (j=0; j<shape_c[1]; j++)
               {
                  for (k=0; k<shape_c[2]; k++)
                  {
                     ind=i+j*shape_c[0]+k*shape_c[0]*shape_c[1];
                     (gfs[h6-1])[ind]=coarsest;
                  }
               }
            }

            valid=PAMR_next_g();
         }
      }

      if (SAVE_MGH_TEST) PAMR_save_gfn(f1_name,PAMR_MGH,1,-1,-1.0,pre_tag,"_after_build");
      if (SAVE_MGH_TEST) PAMR_save_gfn(f2_name,PAMR_MGH,1,-1,-1.0,pre_tag,"_after_build");
      if (SAVE_MGH_TEST) PAMR_save_gfn("f5",PAMR_MGH,1,-1,-1.0,pre_tag,"_after_build");
      if (SAVE_MGH_TEST) PAMR_save_gfn(f6_name,PAMR_MGH,1,-1,-1.0,pre_tag,"_after_build");
      if (SAVE_MGH_TEST && EXCISE_TEST) PAMR_save_gfn("chr",PAMR_MGH,1,-1,-1.0,pre_tag,"_after_build");

      PAMR_destroy_mgh();

      if (rank == 0) printf("... MGH_TEST finished\n");
   }
 
   if (REGRID_TEST) {
     // fill old hierarchy with a function, to see whether stuff transferred
     // correctly
     if (rank==0) printf("REGRID_TEST start\n");
     
     for (l=1; l<4; l++)
       {
         valid=PAMR_init_s_iter(l,PAMR_AMRH,0);
         while(valid)
	   {
	     PAMR_get_g_attribs(&rank,&dim,shape,shape_c,bbox_b,ghost_width6,&t,&ngfs,x,x_c,gfs);
	     sz=1;
	     for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape[i]=1;
	     for (i=0; i<shape[0]; i++)
	       {
		 for (j=0; j<shape[1]; j++)
		   {
		     for (k=0; k<shape[2]; k++)
		       {
			 ind=i+j*shape[0]+k*shape[0]*shape[1];
			 x0=x[0][i];
			 if (dim>1) y0=x[1][j]; else y0=0;
			 if (dim>2) z0=x[2][k]; else z0=0;
			 (gfs[h1-1])[ind]=user_func(x0,y0,z0);
		       }
		   }
               }
	     
	 
	     for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape_c[i]=1;
	     for (i=0; i<shape_c[0]; i++)
	       {
		 for (j=0; j<shape_c[1]; j++)
		   {
		     for (k=0; k<shape_c[2]; k++)
		       {
			 ind=i+j*shape_c[0]+k*shape_c[0]*shape_c[1];
			 x0=x_c[0][i];
			 if (dim>1) y0=x_c[1][j]; else y0=0;
			 if (dim>2) z0=x_c[2][k]; else z0=0;
			 (gfs[h2-1])[ind]=user_func(x0,y0,z0);
		       }
		   }
	       }

	     valid=PAMR_next_g();
	   }
       }
        
     if (SAVE_REGRID_TEST) PAMR_save_gfn(f1_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_before_regrid");
     if (SAVE_REGRID_TEST) PAMR_save_gfn(f2_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_before_regrid");
     
     // compose another set of random of grids on level 3 to test regridding
     // prevent overlapping grids
     // also, let each node contribute a unique set, then test merge

     n = 0;
     for (nn=ngl2+1; nn<(1+ngl2+ngl3); nn++)
       {
	 PAMR_get_dxdt(2,&dx[0],&dt);
	 // Change reference bbox to appropriate level 2 box.
	 i = (nn-1)%ngl2 + 1;

	 if (nn < (2*ngl2 + 1)) {
	   xmin = bbox_list[2*dim*i];
	   xmax = bbox_list[2*dim*i+1];
	   ymin = bbox_list[2*dim*i+2];
	   ymax = bbox_list[2*dim*i+3];
	   zmin = bbox_list[2*dim*i+4];
	   zmax = bbox_list[2*dim*i+5];
	 } else { // use largest remaining available region
	   i_lg = 0;
	   vol_lg = 0.0;
	   if (dim == 1) {
	     xmin = avail3[i][0].xmin;
	     xmax = avail3[i][0].xmax;
	     available[0].unoccupied = 0;
	   } else if (dim == 2) {
	     for (j=0;j<3;j++) {
	       if (avail3[i][j].vol*avail3[i][j].unoccupied > vol_lg) {
		 vol_lg = avail3[i][j].vol*avail3[i][j].unoccupied;
		 i_lg = j;
	       }
	     }
	     
	     xmin = avail3[i][i_lg].xmin;
	     xmax = avail3[i][i_lg].xmax;
	     ymin = avail3[i][i_lg].ymin;
	     ymax = avail3[i][i_lg].ymax;
	     avail3[i][i_lg].unoccupied = 0;

	   } else if (dim == 3) {
	     for (j=0;j<7;j++) {
	       if (avail3[i][j].vol*avail3[i][j].unoccupied > vol_lg) {
		 vol_lg = avail3[i][j].vol*avail3[i][j].unoccupied;
		 i_lg = j;
	       }
	     }
	     xmin = avail3[i][i_lg].xmin;
	     xmax = avail3[i][i_lg].xmax;
	     ymin = avail3[i][i_lg].ymin;
	     ymax = avail3[i][i_lg].ymax;
	     zmin = avail3[i][i_lg].zmin;
	     zmax = avail3[i][i_lg].zmax;
	     avail3[i][i_lg].unoccupied = 0;
	   }
	 }
	 
	 x1=((real)rand())/RAND_MAX;
	 x2=((real)rand())/RAND_MAX;
	 if (x1>x2) {t=x1; x1=x2; x2=t;}
	 tmp_bbox_list[2*dim*nn]=max((int)(x1*(xmax-xmin)/dx[0])*dx[0]+xmin,xmin+2*dx[0]);
	 tmp_bbox_list[2*dim*nn+1]=min((int)(x2*(xmax-xmin)/dx[0])*dx[0]+xmin,xmax-2*dx[0]);

	 if (tmp_bbox_list[2*dim*nn]>tmp_bbox_list[2*dim*nn+1]) {
	   t=tmp_bbox_list[2*dim*nn];
	   tmp_bbox_list[2*dim*nn]   = tmp_bbox_list[2*dim*nn+1];
	   tmp_bbox_list[2*dim*nn+1] = t;
	 }

	 if ((tmp_bbox_list[2*dim*nn+1] -  tmp_bbox_list[2*dim*nn]) < 4*dx[0]) {
	   n_exp = 4 - (tmp_bbox_list[2*dim*nn+1]-tmp_bbox_list[2*dim*nn])/dx[0];
	   n_exp = (n_exp + 1)/2;
	   udiff = (xmax - 2*dx[0] - tmp_bbox_list[2*dim*nn+1])/dx[0] + 0.01;
	   ldiff = (tmp_bbox_list[2*dim*nn] - xmin - 2*dx[0])/dx[0] + 0.01;
	   lfac = 0.0;
	   ufac = 0.0;
	   if (udiff >= n_exp && ldiff >= n_exp) {
	     lfac = n_exp;
	     ufac = n_exp;
	   } else if (udiff < n_exp && ldiff >= n_exp) {
	     lfac = min(ldiff,(2*n_exp-udiff));
	     ufac = udiff;
	   } else if (udiff >= n_exp && ldiff < n_exp) {
	     lfac = ldiff;
	     ufac = min(udiff,(2*n_exp-ldiff));
	   } else if (udiff < n_exp && ldiff < n_exp) {
	     // expand as much as you can
	     lfac = max(ldiff,0);
	     ufac = max(udiff,0);
	   }
	   // Try to expand it
	   tmp_bbox_list[2*dim*nn] -= lfac*dx[0];
	   tmp_bbox_list[2*dim*nn+1] += ufac*dx[0];
	 }

	 if ((nn%size)==rank) {
	   bbox_list[2*dim*n] = tmp_bbox_list[2*dim*nn];
	   bbox_list[2*dim*n+1] = tmp_bbox_list[2*dim*nn+1];
	   printf("node %i has grid %i\n",rank,nn);
	   printf("bbox=[%lf,%lf",tmp_bbox_list[2*dim*nn],tmp_bbox_list[2*dim*nn+1]);
	 }

	 // For good measure....
	 tmp_bbox_list[2*dim*nn+2] = ymin;
	 tmp_bbox_list[2*dim*nn+3] = ymax;
	 tmp_bbox_list[2*dim*nn+4] = zmin;
	 tmp_bbox_list[2*dim*nn+5] = zmax;
	 if (dim>1)
	   {
	     x1=((real)rand())/RAND_MAX;
	     x2=((real)rand())/RAND_MAX;
	     if (x1>x2) {t=x1; x1=x2; x2=t;}
	     tmp_bbox_list[2*dim*nn+2]=max((int)(x1*(ymax-ymin)/dx[1])*dx[1]+ymin,ymin+2*dx[1]);
	     tmp_bbox_list[2*dim*nn+3]=min((int)(x2*(ymax-ymin)/dx[1])*dx[1]+ymin,ymax-2*dx[1]);

	     if (tmp_bbox_list[2*dim*nn+2]>tmp_bbox_list[2*dim*nn+3]) {
	       t=tmp_bbox_list[2*dim*nn+2];
	       tmp_bbox_list[2*dim*nn+2] = tmp_bbox_list[2*dim*nn+3];
	       tmp_bbox_list[2*dim*nn+3] = t;
	     }

	     if ((tmp_bbox_list[2*dim*nn+3] -  tmp_bbox_list[2*dim*nn+2]) < 4*dx[1]) {
	       n_exp = 4 - (tmp_bbox_list[2*dim*nn+3]-tmp_bbox_list[2*dim*nn+2])/dx[1];
	       n_exp = (n_exp + 1)/2;
	       udiff = (ymax - 2*dx[1] - tmp_bbox_list[2*dim*nn+3])/dx[1] + 0.01;
	       ldiff = (tmp_bbox_list[2*dim*nn+2] - ymin - 2*dx[1])/dx[1] + 0.01;
	       lfac = 0.0;
	       ufac = 0.0;
	       if (udiff >= n_exp && ldiff >= n_exp) {
		 lfac = n_exp;
		 ufac = n_exp;
	       } else if (udiff < n_exp && ldiff >= n_exp) {
		 lfac = min(ldiff,(2*n_exp-udiff));
		 ufac = udiff;
	       } else if (udiff >= n_exp && ldiff < n_exp) {
		 lfac = ldiff;
		 ufac = min(udiff,(2*n_exp-ldiff));
	       } else if (udiff < n_exp && ldiff < n_exp) {
		 // expand as much as you can
		 lfac = max(ldiff,0);
		 ufac = max(udiff,0);
	       }
	       // Try to expand it
	       tmp_bbox_list[2*dim*nn+2] -= lfac*dx[1];
	       tmp_bbox_list[2*dim*nn+3] += ufac*dx[1];
	     }

	     if ((nn%size)==rank) {
	       bbox_list[2*dim*n+2] = tmp_bbox_list[2*dim*nn+2];
	       bbox_list[2*dim*n+3] = tmp_bbox_list[2*dim*nn+3];
       
	       printf("][%lf,%lf",tmp_bbox_list[2*dim*nn+2],tmp_bbox_list[2*dim*nn+3]);
	     }

	   }
	 if (dim>2)
	   {
	     x1=((real)rand())/RAND_MAX;
	     x2=((real)rand())/RAND_MAX;
	     if (x1>x2) {t=x1; x1=x2; x2=t;}
	     tmp_bbox_list[2*dim*nn+4]=max((int)(x1*(zmax-zmin)/dx[2])*dx[2]+zmin,zmin+2*dx[2]);
	     tmp_bbox_list[2*dim*nn+5]=min((int)(x2*(zmax-zmin)/dx[2])*dx[2]+zmin,zmax-2*dx[2]);

	     if (tmp_bbox_list[2*dim*nn+4]>tmp_bbox_list[2*dim*nn+5]) {
	       t=tmp_bbox_list[2*dim*nn+4];
	       tmp_bbox_list[2*dim*nn+4] = tmp_bbox_list[2*dim*nn+5];
	       tmp_bbox_list[2*dim*nn+5] = t;
	     }

	     // Try to prevent boxes which are too small in the z direction
	     if ((tmp_bbox_list[2*dim*nn+5] -  tmp_bbox_list[2*dim*nn+4]) < 4*dx[2]) {
	       n_exp = 4 - (tmp_bbox_list[2*dim*nn+5]-tmp_bbox_list[2*dim*nn+4])/dx[2];
	       n_exp = (n_exp + 1)/2;
	       udiff = (zmax - 2*dx[2] - tmp_bbox_list[2*dim*nn+5])/dx[2] + 0.01;
	       ldiff = (tmp_bbox_list[2*dim*nn+4] - zmin - 2*dx[2])/dx[2] + 0.01;
	       lfac = 0.0;
	       ufac = 0.0;
	       if (udiff >= n_exp && ldiff >= n_exp) {
		 lfac = n_exp;
		 ufac = n_exp;
	       } else if (udiff < n_exp && ldiff >= n_exp) {
		 lfac = min(ldiff,(2*n_exp-udiff));
		 ufac = udiff;
	       } else if (udiff >= n_exp && ldiff < n_exp) {
		 lfac = ldiff;
		 ufac = min(udiff,(2*n_exp-ldiff));
	       } else if (udiff < n_exp && ldiff < n_exp) {
		 // expand as much as you can
		 lfac = max(ldiff,0);
		 ufac = max(udiff,0);
	       }
	       // Try to expand it
	       tmp_bbox_list[2*dim*nn+4] -= lfac*dx[2];
	       tmp_bbox_list[2*dim*nn+5] += ufac*dx[2];
	     }

	     if ((nn%size)==rank) {
	       bbox_list[2*dim*n+4] = tmp_bbox_list[2*dim*nn+4];
	       bbox_list[2*dim*n+5] = tmp_bbox_list[2*dim*nn+5];
	       
	       printf("][%lf,%lf",tmp_bbox_list[2*dim*nn+4],tmp_bbox_list[2*dim*nn+5]);
	     }
	   }

	 if ((nn%size)==rank) printf("]\n");
	 if ((nn%size)==rank) ++n;
	 lev_list[n]=3; 
	 
	 // If this is the first visit to the given level 2 grid, then 
	 // calculate the available boxes.
	 if (nn < (2*ngl2 + 1)) {

	   // Figure out x_cut, y_cut, and z_cut.
	   if ((tmp_bbox_list[2*dim*nn]-xmin)/(xmax-xmin) > (xmax-tmp_bbox_list[2*dim*nn+1])/(xmax-xmin)) {
	     x_cut = tmp_bbox_list[2*dim*nn];
	     x_flag = 0;
	   } else {
	     x_cut = tmp_bbox_list[2*dim*nn+1];
	     x_flag = 1;
	   }
	   
	   if (dim > 1) {
	     if ((tmp_bbox_list[2*dim*nn+2]-ymin)/(ymax-ymin) > (ymax-tmp_bbox_list[2*dim*nn+3])/(ymax-ymin)) {
	       y_cut = tmp_bbox_list[2*dim*nn+2];
	       y_flag = 0;
	     } else {
	       y_cut = tmp_bbox_list[2*dim*nn+3];
	       y_flag = 1;
	     }
	   }

	   if (dim > 2) {
	     if ((tmp_bbox_list[2*dim*nn+4]-zmin)/(zmax-zmin) > (zmax-tmp_bbox_list[2*dim*nn+5])/(zmax-zmin)) {
	       z_cut = tmp_bbox_list[2*dim*nn+4];
	       z_flag = 0;
	     } else {
	       z_cut = tmp_bbox_list[2*dim*nn+5];
	       z_flag = 1;
	     }
	   }
	 
	   // Make list of available boxes.
	   if (dim == 1) {
	     // One dimension -> one available box
	     avail3[i][0].unoccupied = 1;

	     avail3[i][0].xmin = ( x_flag == 0 ? x_cut : xmin  );
	     avail3[i][0].xmax = ( x_flag == 0 ? xmax  : x_cut );
	     avail3[i][0].ymin = 0.0;
	     avail3[i][0].ymax = 0.0;
	     avail3[i][0].zmin = 0.0;
	     avail3[i][0].zmax = 0.0;
	     avail3[i][0].vol = avail3[i][0].xmax-avail3[i][0].xmin;
	   } else if (dim == 2) {
	     // Two dimensions -> three available boxes

	     for (j=0; j<3; j++) avail3[i][j].unoccupied = 1;
	     for (j=0; j<3; j++) avail3[i][j].zmin = 0.0;
	     for (j=0; j<3; j++) avail3[i][j].zmax = 0.0;

	     avail3[i][0].xmin = ( x_flag==1 ? x_cut : xmin  );
	     avail3[i][0].xmax = ( x_flag==1 ? xmax  : x_cut );
	     avail3[i][1].xmin = ( x_flag==1 ? x_cut : xmin  );
	     avail3[i][1].xmax = ( x_flag==1 ? xmax  : x_cut );
	     avail3[i][2].xmin = ( x_flag==1 ? xmin  : x_cut );
	     avail3[i][2].xmax = ( x_flag==1 ? x_cut : xmax  );

	     avail3[i][0].ymin = ymin;
	     avail3[i][0].ymax = y_cut;
	     avail3[i][1].ymin = y_cut;
	     avail3[i][1].ymax = ymax;
	     avail3[i][2].ymin = ( y_flag==1 ? y_cut : ymin  );
	     avail3[i][2].ymax = ( y_flag==1 ? ymax  : y_cut );

	     // Calculate volumes
	     for (j=0;j<3;j++) avail3[i][j].vol = (avail3[i][j].xmax-avail3[i][j].xmin) *
                                	         (avail3[i][j].ymax-avail3[i][j].ymin);

	   } else {
	     // Three dimensions -> seven available boxes
	     for (j=0; j<7; j++) avail3[i][j].unoccupied = 1;

	     // The group of 4
	     avail3[i][0].xmin = xmin;
	     avail3[i][0].xmax = x_cut;
	     avail3[i][0].ymin = ymin;
	     avail3[i][0].ymax = y_cut;

	     avail3[i][1].xmin = xmin;
	     avail3[i][1].xmax = x_cut;
	     avail3[i][1].ymin = y_cut;
	     avail3[i][1].ymax = ymax;

	     avail3[i][2].xmin = x_cut;
	     avail3[i][2].xmax = xmax;
	     avail3[i][2].ymin = ymin;
	     avail3[i][2].ymax = y_cut;

	     avail3[i][3].xmin = x_cut;
	     avail3[i][3].xmax = xmax;
	     avail3[i][3].ymin = y_cut;
	     avail3[i][3].ymax = ymax;

	     for (j=0; j<4; j++) {
	       avail3[i][j].zmin = ( z_flag==1 ? z_cut : zmin  );
	       avail3[i][j].zmax = ( z_flag==1 ? zmax  : z_cut );
	     }

	     // The group of 3
	     avail3[i][4].xmin = ( x_flag==1 ? x_cut : xmin  );
	     avail3[i][4].xmax = ( x_flag==1 ? xmax  : x_cut );
	     avail3[i][5].xmin = ( x_flag==1 ? x_cut : xmin  );
	     avail3[i][5].xmax = ( x_flag==1 ? xmax  : x_cut );
	     avail3[i][6].xmin = ( x_flag==1 ? xmin  : x_cut );
	     avail3[i][6].xmax = ( x_flag==1 ? x_cut : xmax  );

	     avail3[i][4].ymin = ymin;
	     avail3[i][4].ymax = y_cut;
	     avail3[i][5].ymin = y_cut;
	     avail3[i][5].ymax = ymax;
	     avail3[i][6].ymin = ( y_flag==1 ? y_cut : ymin  );
	     avail3[i][6].ymax = ( y_flag==1 ? ymax  : y_cut );

	     for (j=4; j<7; j++) {
	       avail3[i][j].zmin = ( z_flag==1 ? zmin  : z_cut );
	       avail3[i][j].zmax = ( z_flag==1 ? z_cut : zmax  );
	     }

	     // Calculate volumes
	     for (j=0;j<7;j++) avail3[i][j].vol = (avail3[i][j].xmax-avail3[i][j].xmin) *
                                        	 (avail3[i][j].ymax-avail3[i][j].ymin) *
                                	         (avail3[i][j].zmax-avail3[i][j].zmin);
	   } // end if dim==1
	 } // end if (n < 2*ngl2 + 1)
       } // end loop on n for level 3 grids.
     
     mn=MAX_GRIDS;
     PAMR_merge_bboxes(bbox_list,n,merged_bbox_list,&mn,min_eff);
     for (k=0; k<mn; k++) lev_list[k]=3;
     printf("Node %i, level 3: Number of grids/merged grids = %i/%i\n",rank,n,mn);
     PAMR_compose_hierarchy(3,3,mn,lev_list,merged_bbox_list,0.0);

     if (SAVE_REGRID_TEST) PAMR_save_gfn(f1_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_after_regrid");
     if (SAVE_REGRID_TEST) PAMR_save_gfn(f2_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_after_regrid");
     
     // Next calculate errors.
     c_norm = 0.0;
     c_counter = 0;
     c_norml = 0.0;
     c_counterl = 0;
     v_norm = 0.0;
     v_counter = 0;
     v_norml = 0.0;
     v_counterl = 0;
     
     for (l=1; l<4; l++) {
       
       valid=PAMR_init_s_iter(l,PAMR_AMRH,0);
       while(valid)
	 {
	   PAMR_get_g_attribs(&rank,&dim,shape,shape_c,bbox_b,ghost_width6,&t,&ngfs,x,x_c,gfs);
	   sz=1;
	   
	   for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape_c[i]=1;
	   for (i=0; i<PAMR_MAX_DIM; i++) if (dim<(i+1)) shape[i]=1;
	   
	   ilo = 0;
	   jlo = 0;
	   klo = 0;
	   
	   if (ghost_width6[0]) ilo = ghost_width6[0] + 1;
	   iup = shape[0] - ghost_width6[1];
	   if (dim > 1) {
	     if (ghost_width6[2]) jlo = ghost_width6[2] + 1;
	     jup = shape[1] - ghost_width6[3];
	   } else {
	     jup = shape[1];
	   }
	   if (dim > 2) {
	     if (ghost_width6[4]) klo = ghost_width6[4] + 1;
	     kup = shape[2] - ghost_width6[5];
	   } else {
	     kup = shape[2];
	   }
	   
	   // Do the following for vertex centered 
	   for (i=0; i<shape[0]; i++)
	     {
	       for (j=0; j<shape[1]; j++)
		 {
		   for (k=0; k<shape[2]; k++)
		      {
			
			ind=i+j*shape[0]+k*shape[0]*shape[1];
			x0=x[0][i];
			if (dim>1) y0=x[1][j]; else y0=0;
			if (dim>2) z0=x[2][k]; else z0=0;
			(gfs[h1-1])[ind] -= user_func(x0,y0,z0);
			(gfs[h1-1])[ind] /= user_func(x0,y0,z0);
			if ( fabs((gfs[h1-1])[ind]) > 0.1) {
			  printf("Warning, large errors at i=%d, j=%d, k=%d on proc %d\n",i,j,k,rank);
			}

			if ( (i>=ilo && i<iup) && (j>=jlo && j<jup) && (k>=klo && k<kup)) {		      
			  v_norml += (gfs[h1-1])[ind]*(gfs[h1-1])[ind];
			  v_counterl++;
			}
		      }
		 }
	     }
	   
	   ilo = 0;
	   jlo = 0;
	   klo = 0;
	   ilo = ghost_width6[0];
	   iup = shape_c[0] - ghost_width6[1];	  
	   if (dim > 1) {
	     jlo = ghost_width6[2];
	     jup = shape_c[1] - ghost_width6[3];
	   } else {
	     jlo = 0;
	     jup = shape_c[1];
	   }
	   if (dim > 2) {
	     klo = ghost_width6[4];
	     kup = shape_c[2] - ghost_width6[5];
	   } else {
	     klo = 0;
	     kup = shape_c[2];
	   }
	   
	   for (i=0; i<shape_c[0]; i++)
	     {
	       for (j=0; j<shape_c[1]; j++)
		 {
		   for (k=0; k<shape_c[2]; k++)
		     {
		       ind=i+j*shape_c[0]+k*shape_c[0]*shape_c[1];
		       x0=x_c[0][i];
		       if (dim>1) y0=x_c[1][j]; else y0=0;
		       if (dim>2) z0=x_c[2][k]; else z0=0;
		       (gfs[h2-1])[ind] -= user_func(x0,y0,z0);
		       (gfs[h2-1])[ind] /= user_func(x0,y0,z0);
		       if ( fabs((gfs[h2-1])[ind]) > 0.1) {
			 printf("Warning, large errors at i=%d, j=%d, k=%d on proc %d\n",i,j,k,rank);
		       }

		       if ( (i>=ilo && i<iup) && (j>=jlo && j<jup) && (k>=klo && k<kup)) {
			 c_norml += (gfs[h2-1])[ind]*(gfs[h2-1])[ind];
			 c_counterl++;
		       }
		     }
		 }
	     }

	   valid=PAMR_next_g();
	 }
     }
     
     MPI_Allreduce(&c_norml,&c_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
     MPI_Allreduce(&c_counterl,&c_counter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
     MPI_Allreduce(&v_norml,&v_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
     MPI_Allreduce(&v_counterl,&v_counter,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
     if (rank == 0) {
       printf("l = %d, c_error_norm = %20.15e\n", l, sqrt(c_norm/c_counter));
       printf("l = %d, v_error_norm = %20.15e\n", l, sqrt(v_norm/v_counter));
     }
     
     if (SAVE_REGRID_TEST) PAMR_save_gfn(f1_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_error");
     if (SAVE_REGRID_TEST) PAMR_save_gfn(f2_name,PAMR_AMRH,1,-1,-1.0,pre_tag,"_error");
     
     if (rank==0) printf("... REGRID_TEST finished\n");
   }

   PAMR_tick(1); PAMR_tick(2); PAMR_tick(3);

cleanup:
   if (cnt>=0) PAMR_free_context(cnt);
   if (rank==0) printf("Calling MPI_Finalize...");
   MPI_Finalize();
   if (rank==0) printf("done\nexiting\n");
} 

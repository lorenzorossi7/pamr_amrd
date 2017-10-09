//=============================================================================
//
// test1.c ---- A set of tests for pamr
//
//=============================================================================

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include "pamr.h"
#include "../src/misc.h"
#include <math.h>

#define MAX_GRIDS 100
#define EXCISED 1
#define NOT_EXCISED 0
//------------------------------------------------------------------------------
// for excision ... excised a random sphere from the computational domain
//------------------------------------------------------------------------------
real ex_R,ex_x,ex_y,ex_z;
void set_ex_mask(real *mask, int dim, int *shape, real *bbox, real excised)
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

int main(int argc, char **argv)
{
   int rank,size,cnt=-1,shape[3]={65,33,33},i,j,l,k,h,ind,h1,h3,h4,h5,nn,mn,sz;
   real bbox[6]={-1,1,-2,2,-0.5,0.5},lambda,bbox_b[6];
   real bbox_list[6*MAX_GRIDS];
   real merged_bbox_list[6*MAX_GRIDS];
   int lev_list[MAX_GRIDS];
   int levels=3;
   int rho_sp[3]={2,4,2},rho_tm[3]={2,2,4};
   int grho_sp[3],grho_tm[3];
   int ghost_width[3]={4,1,3},min_width[3]={9,7,13};
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
   int f1_n,f2_n,f3_n,f4_n,f5_n,f6_n,chr_n;
   int ltrace=0,valid,ngfs;
   int ngl2,ngl3,seed,n,dim,gdm,nt,max_lev,coarsest;
   real dt,dx[3];
   real *x[3],*gfs[500],x0,y0,z0;
   int periodic[3]={1,1,1};

   int SYNC_TEST=1,SAVE_SYNC_TEST=1;
   int INTERP_TEST=0,SAVE_INTERP_TEST=1;
   int BDY_INTERP_TEST=0,SAVE_BDY_INTERP_TEST=1;
   int INJECT_TEST=0,SAVE_INJECT_TEST=1;
   int REGRID_TEST=0,SAVE_REGRID_TEST=1;
   int MGH_TEST=0,SAVE_MGH_TEST=1;
   int EXCISE_TEST=0;

   printf("Calling MPI_Init...\n");
   MPI_Init(&argc,&argv);
   printf("done\n");

   if (argc<5 || argc>6)
   {
      printf("usage: %s dim ngl2 ngl3 seed [min_eff=1]\n",argv[0]);
      goto cleanup;
   }
   dim=atoi(argv[1]);
   ngl2=atoi(argv[2]);
   ngl3=atoi(argv[3]);
   seed=atoi(argv[4]);
   if (argc==6) min_eff=atof(argv[5]); else min_eff=1;
   
   printf("%s dim=%i, ngl2=%i, ngl3=%i, seed=%i\n",argv[0],dim,ngl2,ngl3,seed);
   srand(seed);

   MPI_Comm_size(MPI_COMM_WORLD,&size);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   printf("MPI size =%i, (my) rank=%i\n",size,rank);

   PAMR_set_trace_lev(4);

   cnt=PAMR_init_context(0,0,dim,shape,bbox,0,0);
   printf("Init context returns %i\n",cnt);

   PAMR_set_periodic_bdy(periodic);

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
   IFL printf("MG_max_width=[%i,%i,%i]\n",gMG_max_cwidth[0],gMG_max_cwidth[1],gMG_max_cwidth[2]);

   f1_n=PAMR_def_var_full("f1",1,1,4,PAMR_HW_RESTR,PAMR_FOURTH_ORDER,PAMR_FOURTH_ORDER,PAMR_SYNC,
                         PAMR_HW_RESTR,PAMR_FOURTH_ORDER,PAMR_SYNC,1,PAMR_NO_TRANSFER,phys_bnd1);
   f2_n=PAMR_def_var_full("f2",1,1,4,PAMR_NO_INJECT,PAMR_FOURTH_ORDER,PAMR_FOURTH_ORDER,PAMR_SYNC,
                         PAMR_HW_RESTR,PAMR_SECOND_ORDER,PAMR_SYNC,1,PAMR_SECOND_ORDER,phys_bnd1);
   f3_n=PAMR_def_var_full("f3",1,0,2,PAMR_STRAIGHT_INJECT,PAMR_NO_INTERP,PAMR_NO_INTERP,PAMR_NO_SYNC,
                         PAMR_STRAIGHT_INJECT,PAMR_NO_INTERP,PAMR_NO_SYNC,0,PAMR_FOURTH_ORDER,phys_bnd2);
   f4_n=PAMR_def_var_full("f4",1,0,2,PAMR_FW_RESTR,PAMR_SECOND_ORDER,PAMR_SECOND_ORDER,PAMR_NO_SYNC,
                         PAMR_FW_RESTR,PAMR_SECOND_ORDER,PAMR_NO_SYNC,0,PAMR_NO_TRANSFER,phys_bnd2);
   f5_n=PAMR_def_var_full("f5",0,1,0,PAMR_NO_INJECT,PAMR_SECOND_ORDER,PAMR_SECOND_ORDER,PAMR_SYNC,
                         PAMR_SECOND_ORDER,PAMR_SECOND_ORDER,PAMR_SYNC,1,PAMR_NO_TRANSFER,phys_bnd3);
   f6_n=PAMR_def_var_brief("f6");

   // MAKE SURE MASK IS LAST GFN
   chr_n=PAMR_def_var_full("chr",1,1,1,PAMR_NO_INJECT,PAMR_NO_INTERP,PAMR_NO_INTERP,PAMR_NO_SYNC,
                            PAMR_NO_INJECT,PAMR_NO_INTERP,PAMR_NO_SYNC,1,PAMR_NO_TRANSFER,phys_bnd3);

   if (EXCISE_TEST)
   {
      if (rank==0) printf("testing excision, chr_n=%i\n",chr_n);
      PAMR_excision_on("chr",set_ex_mask,EXCISED,1);
   }

   IFL printf("starting grid function numbers: f1,f2,f3,f4,f5,f6:%i,%i,%i,%i,%i,%i\n",
          f1_n,f2_n,f3_n,f4_n,f5_n,f6_n);

   IFL printf("\nattributes for f2:\n");
   f2_n=PAMR_get_var_attribs("f2",&in_amrh,&in_mgh,&num_tl,&amr_inject,
                        &amr_interp,&amr_bdy_interp,&amr_sync,&mg_inject,&mg_interp,&mg_sync,
                        &mg_noinj_to_amr,&regrid_transfer,phys_bdy_type);
   IFL printf("f2_n=%i,   in_amrh=%i,   in_mgh=%i,   num_tl=%i,   amr_inject=%i\n"
          "amr_interp=%i,   amr_sync=%i,    mg_inject=%i,   mg_interp=%i\n"
          "mg_sync=%i,   mg_noinj_to_amr=%i,   regrid_transfer=%i\n"
          "phys_bdy_type=[%i,%i,%i,%i,%i,%i]\n",f2_n,in_amrh,in_mgh,num_tl,amr_inject,
          amr_interp,amr_sync,mg_inject,mg_interp,mg_sync,mg_noinj_to_amr,regrid_transfer,
          phys_bdy_type[0],phys_bdy_type[1],phys_bdy_type[2],phys_bdy_type[3],phys_bdy_type[4],
          phys_bdy_type[5]);

   IFL printf("\nattributes for f4:\n");
   f4_n=PAMR_get_var_attribs("f4",&in_amrh,&in_mgh,&num_tl,&amr_inject,
                        &amr_interp,&amr_bdy_interp,&amr_sync,&mg_inject,&mg_interp,&mg_sync,
                        &mg_noinj_to_amr,&regrid_transfer,phys_bdy_type);
   IFL printf("f4_n=%i,   in_amrh=%i,   in_mgh=%i,   num_tl=%i,   amr_inject=%i\n"
          "amr_interp=%i,   amr_sync=%i,    mg_inject=%i,   mg_interp=%i\n"
          "mg_sync=%i,   mg_noinj_to_amr=%i,   regrid_transfer=%i\n"
          "phys_bdy_type=[%i,%i,%i,%i,%i,%i]\n",f4_n,in_amrh,in_mgh,num_tl,amr_inject,
          amr_interp,amr_sync,mg_inject,mg_interp,mg_sync,mg_noinj_to_amr,regrid_transfer,
          phys_bdy_type[0],phys_bdy_type[1],phys_bdy_type[2],phys_bdy_type[3],phys_bdy_type[4],
          phys_bdy_type[5]);

   IFL printf("\nattributes for f6:\n");
   f6_n=PAMR_get_var_attribs("f6",&in_amrh,&in_mgh,&num_tl,&amr_inject,
                        &amr_interp,&amr_bdy_interp,&amr_sync,&mg_inject,&mg_interp,&mg_sync,
                        &mg_noinj_to_amr,&regrid_transfer,phys_bdy_type);
   IFL printf("f6_n=%i,   in_amrh=%i,   in_mgh=%i,   num_tl=%i,   amr_inject=%i\n"
          "amr_interp=%i,   amr_sync=%i,    mg_inject=%i,   mg_interp=%i\n"
          "mg_sync=%i,   mg_noinj_to_amr=%i,   regrid_transfer=%i\n"
          "phys_bdy_type=[%i,%i,%i,%i,%i,%i]\n",f6_n,in_amrh,in_mgh,num_tl,amr_inject,
          amr_interp,amr_sync,mg_inject,mg_interp,mg_sync,mg_noinj_to_amr,regrid_transfer,
          phys_bdy_type[0],phys_bdy_type[1],phys_bdy_type[2],phys_bdy_type[3],phys_bdy_type[4],
          phys_bdy_type[5]);

   // compose a random set of grids
   lev_list[0]=1;
   for (n=0; n<6; n++) bbox_list[n]=bbox[n];
   for (n=1; n<(1+ngl2+ngl3); n++)
   {
      if (n>=(ngl2+1)) PAMR_get_dxdt(2,&dx[0],&dt);
      else PAMR_get_dxdt(1,&dx[0],&dt);

      x1=((real)rand())/RAND_MAX;
      x2=((real)rand())/RAND_MAX;
      if (x1>x2) {t=x1; x1=x2; x2=t;}
      bbox_list[2*dim*n]=max((int)(x1*(bbox[1]-bbox[0])/dx[0])*dx[0]+bbox[0]-dx[0],bbox[0]);
      bbox_list[2*dim*n+1]=min((int)(x2*(bbox[1]-bbox[0])/dx[0])*dx[0]+bbox[0]+dx[0],bbox[1]);
      if (dim>1)
      {
         x1=((real)rand())/RAND_MAX;
         x2=((real)rand())/RAND_MAX;
         if (x1>x2) {t=x1; x1=x2; x2=t;}
         bbox_list[2*dim*n+2]=max((int)(x1*(bbox[3]-bbox[2])/dx[1])*dx[1]+bbox[2]-dx[1],bbox[2]);
         bbox_list[2*dim*n+3]=min((int)(x2*(bbox[3]-bbox[2])/dx[1])*dx[1]+bbox[2]+dx[1],bbox[3]);
      }
      if (dim>2)
      {
         x1=((real)rand())/RAND_MAX;
         x2=((real)rand())/RAND_MAX;
         if (x1>x2) {t=x1; x1=x2; x2=t;}
         bbox_list[2*dim*n+4]=max((int)(x1*(bbox[5]-bbox[4])/dx[2])*dx[2]+bbox[4]-dx[2],bbox[4]);
         bbox_list[2*dim*n+5]=min((int)(x2*(bbox[5]-bbox[4])/dx[2])*dx[2]+bbox[4]+dx[2],bbox[5]);
      }
      if (n>=(ngl2+1)) lev_list[n]=3; else lev_list[n]=2;
   }

   PAMR_compose_hierarchy(1,3,1+ngl2+ngl3,lev_list,bbox_list,0.0);

   if (EXCISE_TEST) PAMR_save_gfn("chr",PAMR_AMRH,1,-1,-1.0,"","_after_compose");

   // for sync. tests ... set up grid function data so that 
   // each grid has a unique value prior to sync.
   
   if (SYNC_TEST)
   {
      printf("node %i SYNC_TEST start ... \n",rank);
      for (l=1; l<4; l++)
      {
         nt=0;
         valid=PAMR_init_s_iter(l,PAMR_AMRH,0); while(valid) {valid=PAMR_next_g(); nt++;}
         valid=PAMR_init_s_iter(l,PAMR_AMRH,0);
         n=0;
         while(valid)
         {
            PAMR_get_g_attribs(&rank,&dim,shape,bbox_b,ghost_width,&t,&ngfs,x,gfs);
            sz=1;
            for (i=0; i<dim; i++) sz*=shape[i];
            for (i=0; i<(ngfs-2); i++) // ngfs-2 to not overwrite chr!!!
            {
               if (gfs[i]) for (j=0; j<sz; j++) (gfs[i])[j]=rank+l/4.0+n/4.0/(nt+1);
            }
            valid=PAMR_next_g();
            n++;
         }
      }

      if (SAVE_SYNC_TEST) PAMR_save_gfn("f1",PAMR_AMRH,1,-1,-1.0,"","_before_sync");

      PAMR_sync(1,1,PAMR_AMRH,0);
      PAMR_sync(2,1,PAMR_AMRH,0);
      PAMR_sync(3,1,PAMR_AMRH,0);

      if (SAVE_SYNC_TEST) PAMR_save_gfn("f1",PAMR_AMRH,1,-1,-1.0,"","_after_sync");

      printf("... node %i SYNC_TEST finished\n",rank);
   }

   // for interpolation test ... set up with cubic polynomial in level 1
   // (4th order should be exact)

   if (INTERP_TEST)
   {
      printf("node %i INTERP_TEST start\n",rank);

      for (l=1; l<4; l++)
      {
         valid=PAMR_init_s_iter(l,PAMR_AMRH,0);
         while(valid)
         {
            PAMR_get_g_attribs(&rank,&dim,shape,bbox_b,ghost_width,&t,&ngfs,x,gfs);
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
                     for (h=0; h<(ngfs-2); h++)
                     {
                        if (l==1 && gfs[h] && h!=(chr_n-1)) (gfs[h])[ind]=x0*x0*x0+y0*y0*y0+z0*z0*z0; 
                        else if (gfs[h]) (gfs[h])[ind]=0;
                     }
                  }
               }
            }
            valid=PAMR_next_g();
         }
      }

      PAMR_interp(1,1,PAMR_AMRH);
      PAMR_interp(2,1,PAMR_AMRH);

      if (SAVE_INTERP_TEST) PAMR_save_gfn("f1",PAMR_AMRH,1,-1,-1.0,"","_after_interp");
      if (SAVE_INTERP_TEST) PAMR_save_gfn("f4",PAMR_AMRH,1,-1,-1.0,"","_after_interp");

      printf("... node %i INTERP_TEST finished\n",rank);
   }

   // boundary interpolation test

   if (BDY_INTERP_TEST)
   {
      printf("node %i BDY_INTERP_TEST start\n",rank);

      for (l=1; l<4; l++)
      {
         valid=PAMR_init_s_iter(l,PAMR_AMRH,0);
         while(valid)
         {
            PAMR_get_g_attribs(&rank,&dim,shape,bbox_b,ghost_width,&t,&ngfs,x,gfs);
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
                     for (h=0; h<(ngfs-2); h++)
                     {
                        if (l==1 && gfs[h]) (gfs[h])[ind]=x0*x0*x0+y0*y0*y0+z0*z0*z0; 
                        else if (gfs[h]) (gfs[h])[ind]=0;
                     }
                  }
               }
            }
            valid=PAMR_next_g();
         }
      }

      PAMR_AMR_bdy_interp(1,1,1);
      PAMR_AMR_bdy_interp(2,1,2);

      if (SAVE_BDY_INTERP_TEST) PAMR_save_gfn("f1",PAMR_AMRH,1,-1,-1.0,"","_after_bdy_interp");
      if (SAVE_BDY_INTERP_TEST) PAMR_save_gfn("f4",PAMR_AMRH,1,-1,-1.0,"","_after_bdy_interp");

      printf("... node %i INTERP_TEST finished\n",rank);
   }

   if (INJECT_TEST)
   {
      printf("node %i INJECT_TEST start\n",rank);

      h1=PAMR_get_gfn("f1",PAMR_AMRH,1);  // f1 is HW
      h3=PAMR_get_gfn("f3",PAMR_AMRH,1);  // f3 is straight
      h4=PAMR_get_gfn("f4",PAMR_AMRH,1);  // f4 is FW

      for (l=1; l<4; l++)
      {
         valid=PAMR_init_s_iter(l,PAMR_AMRH,0);
         while(valid)
         {
            PAMR_get_g_attribs(&rank,&dim,shape,bbox_b,ghost_width,&t,&ngfs,x,gfs);
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
                     // to test the different stencils:
//                     (gfs[h1-1])[ind]=(gfs[h3-1])[ind]=(gfs[h4-1])[ind]=pow(-1,i)*pow(-1,j)*pow(-1,k);
//                     if (i%2==1 && j%2==1 && k%2==1) (gfs[h1-1])[ind]=(gfs[h3-1])[ind]=(gfs[h4-1])[ind]=0;
                     // for a simpler, visual test
                     if (l==3) (gfs[h1-1])[ind]=(gfs[h3-1])[ind]=(gfs[h4-1])[ind]=x0*x0*x0*x0+y0*y0*y0*y0+z0*z0*z0*z0; 
                     else (gfs[h1-1])[ind]=(gfs[h3-1])[ind]=(gfs[h4-1])[ind]=0;
                     //TEST!!
                     //if (l==2) (gfs[h1-1])[ind]=(gfs[h3-1])[ind]=(gfs[h4-1])[ind]=x0*x0*x0*x0+y0*y0*y0*y0+z0*z0*z0*z0; 
                     //else (gfs[h1-1])[ind]=(gfs[h3-1])[ind]=(gfs[h4-1])[ind]=0;
                  }
               }
            }
            valid=PAMR_next_g();
         }
      }

      if (SAVE_INJECT_TEST) PAMR_save_gfn("f1",PAMR_AMRH,1,-1,-1.0,"","_before_inject");
      if (SAVE_INJECT_TEST) PAMR_save_gfn("f3",PAMR_AMRH,1,-1,-1.0,"","_before_inject");
      if (SAVE_INJECT_TEST) PAMR_save_gfn("f4",PAMR_AMRH,1,-1,-1.0,"","_before_inject");

      PAMR_inject(3,1,PAMR_AMRH);
      PAMR_inject(2,1,PAMR_AMRH);

      if (SAVE_INJECT_TEST) PAMR_save_gfn("f1",PAMR_AMRH,1,-1,-1.0,"","_after_inject");
      if (SAVE_INJECT_TEST) PAMR_save_gfn("f3",PAMR_AMRH,1,-1,-1.0,"","_after_inject");
      if (SAVE_INJECT_TEST) PAMR_save_gfn("f4",PAMR_AMRH,1,-1,-1.0,"","_after_inject");

      printf("... node %i INJECT_TEST finished\n",rank);
   }

   if (MGH_TEST)
   {
      printf("node %i MGH_TEST start\n",rank);

      h1=PAMR_get_gfn("f1",PAMR_AMRH,1);  // f1 is in AMRH and MGH
      h4=PAMR_get_gfn("f4",PAMR_AMRH,1);  // f4 is only in AMRH
      h5=PAMR_get_gfn("f5",PAMR_AMRH,1);  // f5 is only in MGH
      printf("AMRH h5=%i\n",h5);
      h5=PAMR_get_gfn("f5",PAMR_MGH,1);
      printf("MGH h5=%i\n",h5);

      for (l=1; l<4; l++)
      {
         valid=PAMR_init_s_iter(l,PAMR_AMRH,0);
         while(valid)
         {
            PAMR_get_g_attribs(&rank,&dim,shape,bbox_b,ghost_width,&t,&ngfs,x,gfs);
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
                     (gfs[h1-1])[ind]=(gfs[h4-1])[ind]=x0*x0*x0+y0*y0*y0+z0*z0*z0; 
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
            PAMR_get_g_attribs(&rank,&dim,shape,bbox_b,ghost_width,&t,&ngfs,x,gfs);
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
            valid=PAMR_next_g();
         }
      }


      if (SAVE_MGH_TEST) PAMR_save_gfn("f1",PAMR_MGH,1,-1,-1.0,"","_after_build");
      if (SAVE_MGH_TEST) PAMR_save_gfn("f5",PAMR_MGH,1,-1,-1.0,"","_after_build");
      if (SAVE_MGH_TEST && EXCISE_TEST) PAMR_save_gfn("chr",PAMR_MGH,1,-1,-1.0,"","_after_build");

      PAMR_destroy_mgh();

      printf("... node %i MGH_TEST finished\n",rank);
   }
 
   if (REGRID_TEST)
   {
 
      // fill old hierarchy with a function, to see wether stuff transferred
      // correctly
      printf("node %i REGRID_TEST start\n",rank);

      for (l=1; l<4; l++)
      {
         valid=PAMR_init_s_iter(l,PAMR_AMRH,0);
         while(valid)
         {
            PAMR_get_g_attribs(&rank,&dim,shape,bbox_b,ghost_width,&t,&ngfs,x,gfs);
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
                     for (h=0; h<(ngfs-2); h++)
                     {
                        if (gfs[h]) (gfs[h])[ind]=x0*x0*x0+y0*y0*y0+z0*z0*z0; 
                        else if (gfs[h]) (gfs[h])[ind]=0;
                     }
                  }
               }
            }
            valid=PAMR_next_g();
         }
      }

      if (SAVE_REGRID_TEST) PAMR_save_gfn("f1",PAMR_AMRH,1,-1,-1.0,"","_before_regrid");
      if (SAVE_REGRID_TEST) PAMR_save_gfn("f2",PAMR_AMRH,1,-1,-1.0,"","_before_regrid");
      if (SAVE_REGRID_TEST) PAMR_save_gfn("f3",PAMR_AMRH,1,-1,-1.0,"","_before_regrid");

      // compose another set of random of grids to test regridding
      // also, let each node contribute a unique set, then test merge
      n=0;
      for (nn=1; nn<(1+ngl2+ngl3); nn++)
      {
         if ((nn%size)==rank)
         {
            if (nn>=(ngl2+1)) PAMR_get_dxdt(2,&dx[0],&dt);
            else PAMR_get_dxdt(1,&dx[0],&dt);

            x1=((real)rand())/RAND_MAX;
            x2=((real)rand())/RAND_MAX;
            if (x1>x2) {t=x1; x1=x2; x2=t;}
            bbox_list[2*dim*n]=max((int)(x1*(bbox[1]-bbox[0])/dx[0])*dx[0]+bbox[0]-dx[0],bbox[0]);
            bbox_list[2*dim*n+1]=min((int)(x2*(bbox[1]-bbox[0])/dx[0])*dx[0]+bbox[0]+dx[0],bbox[1]);
            printf("node %i has grid %i\n",rank,nn);
            printf("bbox=[%lf,%lf",bbox_list[2*dim*n],bbox_list[2*dim*n+1]);
            if (dim>1)
            {
               x1=((real)rand())/RAND_MAX;
               x2=((real)rand())/RAND_MAX;
               if (x1>x2) {t=x1; x1=x2; x2=t;}
               bbox_list[2*dim*n+2]=max((int)(x1*(bbox[3]-bbox[2])/dx[1])*dx[1]+bbox[2]-dx[1],bbox[2]);
               bbox_list[2*dim*n+3]=min((int)(x2*(bbox[3]-bbox[2])/dx[1])*dx[1]+bbox[2]+dx[1],bbox[3]);
               printf("][%lf,%lf",bbox_list[2*dim*n+2],bbox_list[2*dim*n+3]);
            }
            if (dim>2)
            {
               x1=((real)rand())/RAND_MAX;
               x2=((real)rand())/RAND_MAX;
               if (x1>x2) {t=x1; x1=x2; x2=t;}
               bbox_list[2*dim*n+4]=max((int)(x1*(bbox[5]-bbox[4])/dx[2])*dx[2]+bbox[4]-dx[2],bbox[4]);
               bbox_list[2*dim*n+5]=min((int)(x2*(bbox[5]-bbox[4])/dx[2])*dx[2]+bbox[4]+dx[2],bbox[5]);
               printf("][%lf,%lf",bbox_list[2*dim*n+4],bbox_list[2*dim*n+5]);
            }
            printf("]\n");
            n++;
         }
         else { rand(); rand(); if (dim>1) rand(); rand(); if (dim>2) rand(); rand(); }
         if (nn==(ngl2))
         {
            mn=MAX_GRIDS;
            PAMR_merge_bboxes(bbox_list,n,merged_bbox_list,&mn,min_eff);
            for (k=0; k<mn; k++) lev_list[k]=2;
            printf("Node %i, level 2: Number of grids/merged grids = %i/%i\n",rank,n,mn);
            PAMR_compose_hierarchy(2,2,mn,lev_list,merged_bbox_list,0.0);
            n=0;
         }
      }

      mn=MAX_GRIDS;
      printf("%i\n",n);
      PAMR_merge_bboxes(bbox_list,n,merged_bbox_list,&mn,min_eff);
      for (k=0; k<mn; k++) lev_list[k]=3;
      printf("Node %i, level 3: Number of grids/merged grids = %i/%i\n",rank,n,mn);
      PAMR_compose_hierarchy(3,3,mn,lev_list,merged_bbox_list,0.0);

      if (SAVE_REGRID_TEST) PAMR_save_gfn("f1",PAMR_AMRH,1,-1,-1.0,"","_after_regrid");
      if (SAVE_REGRID_TEST) PAMR_save_gfn("f2",PAMR_AMRH,1,-1,-1.0,"","_after_regrid");
      if (SAVE_REGRID_TEST) PAMR_save_gfn("f3",PAMR_AMRH,1,-1,-1.0,"","_after_regrid");

      printf("... node %i REGRID_TEST finished\n",rank);
   }

   PAMR_tick(1); PAMR_tick(2); PAMR_tick(3);

cleanup:
   if (cnt>=0) PAMR_free_context(cnt);
   printf("Calling MPI_Finalize...");
   MPI_Finalize();
   printf("done\nexiting\n");
} 

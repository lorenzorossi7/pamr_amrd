//=============================================================================
// util.c --- misc. utility routines
//=============================================================================

#include <stdio.h>
#include <math.h>
#include <bbhutil.h>
#include <stdlib.h>
#include <string.h>
#include "util_w.h"
#include "globals_w.h"
#include "regrid_script_w.h"
#include "evolve_w.h"
#include "m_util_r8_w.h"

//=============================================================================
// clean up stuff before we AMRD_stop(), or end normally
//=============================================================================
void clean_up(void)
{
   if (rgs_stream) fclose(rgs_stream);
}

//=============================================================================
// error handler
//=============================================================================
void AMRD_stop(char *err1, char *err2)
{
   clean_up();
   printf("\n%s%s\n\nStopping.\n",err1,err2);
   if (pamr_context) PAMR_free_context(pamr_context);
   MPI_Finalize();
   exit(1);
}

//=============================================================================
// calls application function f for all grids at level L in hierarchy hier
//=============================================================================
void call_app(void (*f)(void), int L, int hier)
{
   int valid;

   if (!f) return;

   valid=PAMR_init_s_iter(L,hier,0);

   while(valid)
   {
      f();
      valid=PAMR_next_g();
   }
}

//=============================================================================
// applies (multiplies for action=0, divides for action=1) weighted average 
// factor wavg to amr_inject_wavg vars for all grids at level L in hierarchy hier
//=============================================================================
void apply_wavg(int action, int L, int hier)
{
   
   if(!AMRD_num_inject_wavg_vars) return;
   
   int i,j;
   int num_wavg_fn;
   int *wavg_gfn;
   int N;
   int valid;	
   real *gfs[PAMR_MAX_GFNS];
   int shape_c[AMRD_dim];
   
   real *f, *wavg_factor;

   valid =PAMR_init_s_iter(L,hier,0);

   /*
   Note the assumption that all the variables with wavg injection are cell-centered
   */
   while(valid)
   {  	
	N = 1;
	PAMR_get_g_shape_c(shape_c);
	for(i=0; i<AMRD_dim; i++) N*=shape_c[i];
 	PAMR_get_g_gfs(gfs);
	
        if(hier==PAMR_MGH){
		wavg_factor = gfs[AMRD_wavg_mg_c_gfn-1];
                wavg_gfn = AMRD_MG_inject_wavg_gfn;
                num_wavg_fn = AMRD_num_MG_inject_wavg_vars;
	} else if(hier==PAMR_AMRH){
		wavg_factor = gfs[AMRD_wavg_c_gfn-1];
                wavg_gfn = AMRD_inject_wavg_gfn;
                num_wavg_fn = AMRD_num_inject_wavg_vars;
	} else {
		AMRD_stop("Invalid hierarchy in apply_wavg() !\n","");
	}
	
	for(i=0; i<num_wavg_fn; i++){
		f = gfs[wavg_gfn[i]-1];
		for(j=0; j<N; j++){
			if(action==0) f[j]*=wavg_factor[j];
			else if(action==1) f[j]/=wavg_factor[j];
		}
	}
   	valid=PAMR_next_g();
   }
}

//=============================================================================
// calls application function f(arg) for all grids at level L in hierarchy hier
//=============================================================================
void call_app_1int(void (*f)(int arg1), int arg, int L, int hier)
{
   int valid;

   if (!f) return;

   valid=PAMR_init_s_iter(L,hier,0);

   while(valid)
   {
      f(arg);

      valid=PAMR_next_g();
   }
}


//=============================================================================
// Versions for calling with the flux correction mask.
//=============================================================================
void call_app_1int_ifc_mask(void (*f)(int arg1, int *ifc_mask), int arg, int L, int hier)
{
   int valid;
   int *ifc_mask;

   if (!f) return;

   valid=PAMR_init_s_iter(L,hier,0);

   while(valid)
   {
     ev_ldptr();
     ifc_mask = 0;
     if (AMRD_fc_mask_gfn) ifc_mask = get_ifc_mask(AMRD_fc_mask,AMRD_g_size_c);
     f(arg,ifc_mask);
     valid=PAMR_next_g();
   }
}

void call_app_ifc_mask(void (*f)(int *ifc_mask), int L, int hier)
{
   int valid;
   int *ifc_mask;

   if (!f) return;

   valid=PAMR_init_s_iter(L,hier,0);

   while(valid)
   {
     ev_ldptr();
     ifc_mask = 0;
     if (AMRD_fc_mask_gfn) ifc_mask = get_ifc_mask(AMRD_fc_mask,AMRD_g_size_c);
     f(ifc_mask);
     valid=PAMR_next_g();
   }
}


//=============================================================================
// as call_app(), but f returns a real ... this function returns 
// some MPI reduction of the local results
//=============================================================================
real call_rr_app(real (*f)(void), int L, int hier, int comm_only, MPI_Op op)
{
   int valid,first=1,comm;
   real res=0,lres=0;

   if (!f) return 0.0e0;

   valid=PAMR_init_s_iter(L,hier,0);

   while(valid)
   {
      PAMR_get_g_comm(&comm);
      if (comm || !(comm_only))
      {
         res=f();
         if (first) lres=res; 
         if (op==MPI_MAX)
         {
            lres=max(res,lres);
         }
         else
         {
            AMRD_stop("call_rr_app: unsupported reduction\n",0);
         }
         first=0;
      }
      valid=PAMR_next_g();
   }

   MPI_Allreduce(&lres,&res,1,MPI_DOUBLE,op,MPI_COMM_WORLD);

   return res;
}

//=============================================================================
// parameter input ... due to some as-of-yet unknown problems on 
// arcturus (SGI cluster), we offer several different options of 
// reading parameters
//=============================================================================
#define USE_SGET_PARAM 0
#define BCAST_PARAMS 1
#define MAX_SFILE_SIZE 100000
char *read_pfile(char *pfile)
{
   static int first=1;
   static char *sfile=0;
   size_t r_num;
   FILE *f;

   if (!first) return sfile;

   first=0;

   if (!(sfile=(char *)malloc(sizeof(char)*MAX_SFILE_SIZE)))
      AMRD_stop("read_pfile: allocating memory for parameter file\n","");

   if (!(f=fopen(pfile,"r"))) AMRD_stop("read_pfile: error opening parameter file\n","");

   r_num=fread(sfile,sizeof(char),MAX_SFILE_SIZE-1,f);

   if (r_num==MAX_SFILE_SIZE-1) AMRD_stop("read_pfile: error ... internal buffer too small\n","");
   sfile[r_num]=0;

   return sfile;
}

void AMRD_int_param(char *pfile, char *name, int *var, int size)
{
   int i,found,eof=0;
   char *sfile,*p;

   if (USE_SGET_PARAM && (my_rank==0 || !BCAST_PARAMS))
   {
      sfile=read_pfile(pfile);
 
      found=0;
      while(!found && *sfile)
      {
         p=sfile;
         while (*p && *p!='\n') p++; 
         if (!(*p)) eof=1; else *p=0;
         found=sget_int_param(sfile,name,var,size);
         if (found==-1) AMRD_stop("Error reading parameter file, variable = ",name);
         if (found==-2) AMRD_stop("Error parsing requested number of values, variable = ",name);
         if (!eof) { *p='\n'; sfile=p+1; }
      }
   }
   else if (my_rank==0 || !BCAST_PARAMS)
   {
      found=get_int_param(pfile,name,var,size);
      if (found==-1) AMRD_stop("Error reading parameter file, variable = ",name);
      if (found==-2) AMRD_stop("Error parsing requested number of values, variable = ",name);
   }

   if (BCAST_PARAMS) MPI_Bcast(var,size,MPI_INT,0,MPI_COMM_WORLD);

   if (AMRD_echo_params && my_rank==0) 
   {
      printf("%s = [%i",name,*var);
      for (i=1; i<size; i++) printf(",%i",var[i]);
      printf("]\n");
   }
}

void AMRD_real_param(char *pfile, char *name, real *var, int size)
{
   int i,found,eof=0;
   char *sfile,*p;

   if (USE_SGET_PARAM && (my_rank==0 || !BCAST_PARAMS))
   {
      sfile=read_pfile(pfile);

      found=0;
      while(!found && *sfile)
      {
         p=sfile;
         while (*p && *p!='\n') p++; 
         if (!(*p)) eof=1; else *p=0;
         found=sget_real_param(sfile,name,var,size);
         if (found==-1) AMRD_stop("Error reading parameter file, variable = ",name);
         if (found==-2) AMRD_stop("Error parsing requested number of values, variable = ",name);
         if (!eof) { *p='\n'; sfile=p+1; }
      }
   }
   else if (my_rank==0 || !BCAST_PARAMS)
   {
      found=get_real_param(pfile,name,var,size);
      if (found==-1) AMRD_stop("Error reading parameter file, variable = ",name);
      if (found==-2) AMRD_stop("Error parsing requested number of values, variable = ",name);
   }

   if (BCAST_PARAMS) MPI_Bcast(var,size,MPI_DOUBLE,0,MPI_COMM_WORLD);

   if (AMRD_echo_params && my_rank==0)
   {
      printf("%s = [%lf",name,*var);
      for (i=1; i<size; i++) printf(",%lf",var[i]);
      printf("]\n");
   }
}

//-----------------------------------------------------------------------------
// !NOTE: with str_param, we assume this routine is only called *once* per name,
// and so clear var[0..size-1]. If this is called more than
// once with the same name, users should "free" the corresponding memory
// before the second call
//-----------------------------------------------------------------------------
void AMRD_str_param(char *pfile, char *name, char **var, int size)
{
   int i,found,eof=0,lsize[AMRD_MAX_VARS];
   char *sfile,*p;

   for (i=0; i<size; i++) var[i]=0;

   if (USE_SGET_PARAM && (my_rank==0 || !BCAST_PARAMS))
   {
      sfile=read_pfile(pfile);

      found=0;
      while(!found && *sfile)
      {
         p=sfile;
         while (*p && *p!='\n') p++; 
         if (!(*p)) eof=1; else *p=0;
         found=sget_str_param(sfile,name,var,size);
         if (found==-1) AMRD_stop("Error reading parameter file, variable = ",name);
         if (found==-2) AMRD_stop("Error parsing requested number of values, variable = ",name);
         if (!eof) { *p='\n'; sfile=p+1; }
      }
   }
   else if (my_rank==0 || !BCAST_PARAMS)
   {
      found=get_str_param(pfile,name,var,size);
      if (found==-1) AMRD_stop("Error reading parameter file, variable = ",name);
      if (found==-2) AMRD_stop("Error parsing requested number of values, variable = ",name);
   }

   if (BCAST_PARAMS) 
   {
      if (my_rank==0) for (i=0; i<size; i++) if (var[i]) lsize[i]=strlen(var[i])+1; else lsize[i]=0;
      MPI_Bcast(lsize,size,MPI_INT,0,MPI_COMM_WORLD);
      for (i=0; i<size; i++)
      {
         if (my_rank!=0)
         {
            var[i]=0;
            if (lsize[i]) if (!(var[i]=(char *)malloc(lsize[i]))) AMRD_stop("AMRD_str_param ... out of memory","");
         }
         if (lsize[i]) MPI_Bcast(var[i],lsize[i],MPI_CHAR,0,MPI_COMM_WORLD);
      }
   }

   if (AMRD_echo_params && my_rank==0)
   {
      if (size==1) printf("%s = [ \"%s\" ]\n",name,var[0]);
      else
      {
         printf("%s = \n[ \"%s\" ",name,var[0]);
         for (i=1; i<size; i++) printf("\n  \"%s\" ",var[i]);
         printf("]\n");
      }
   }
}

void AMRD_ivec_param(char *pfile, char *name, int *var, int size)
{
   int i,found,eof=0,ivsize;
   char *sfile,*p;

   if (USE_SGET_PARAM && (my_rank==0 || !BCAST_PARAMS))
   {
      sfile=read_pfile(pfile);

      found=0;
      while(!found && *sfile)
      {
         p=sfile;
         while (*p && *p!='\n') p++; 
         if (!(*p)) eof=1; else *p=0;
         found=sget_ivec_param(sfile,name,var,size);
         if (found==-1) AMRD_stop("Error reading parameter file, variable = ",name);
         if (found==-2) AMRD_stop("Error parsing requested number of values, variable = ",name);
         if (!eof) { *p='\n'; sfile=p+1; }
      }
   }
   else if (my_rank==0 || !BCAST_PARAMS)
   {
      found=get_ivec_param(pfile,name,var,size);
      if (found==-1) AMRD_stop("Error reading parameter file, variable = ",name);
      if (found==-2) AMRD_stop("Error parsing requested number of values, variable = ",name);
   }

   if (BCAST_PARAMS) 
   {
      // WARNING ... the following uses 'private' info from bbhutil, namely that
      // ivec[0] is the size, and each element takes 4 ints of memory!!
      ivsize=var[0]*4; 
      MPI_Bcast(&ivsize,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(var,ivsize,MPI_INT,0,MPI_COMM_WORLD);
   }

   if (AMRD_echo_params && my_rank==0)
   {
      printf("%s = [ <ivec> ]\n",name);
   }
}

//=============================================================================
// the following versions of the parameter routines are for variable
// sized lists, as supported by the new version of RNPL. USE_SGET_PARAM
// and BCAST_PARAMS are *NOT* supported here (flags are ignored if set)
//=============================================================================
void AMRD_str_param_v(char *pfile, char *name, char ***var, int *size)
{
   int i,found;

   found=get_str_param_v(pfile,name,var,size);
   if (found==-1) AMRD_stop("Error reading parameter file, variable = ",name);

   if (AMRD_echo_params && my_rank==0)
   {
      if (found)
      {
         if (*size==1) printf("%s = [ \"%s\" ]\n",name,(*var)[0]);
         else
         {
            printf("%s = \n[ \"%s\" ",name,(*var)[0]);
            for (i=1; i<*size; i++) printf("\n  \"%s\" ",(*var)[i]);
            printf("] (%i elements in list)\n",*size);
         }
      }
      else printf("%s not found in parameter file %s)\n",name,pfile);
   }
}

void AMRD_real_param_v(char *pfile, char *name, real **var, int *size)
{
   int i,found;

   found=get_real_param_v(pfile,name,var,size);
   if (found==-1) AMRD_stop("Error reading parameter file, variable = ",name);

   if (AMRD_echo_params && my_rank==0)
   {
      if (found)
      {
         if (*size==1) printf("%s = [ \"%e\" ]\n",name,(*var)[0]);
         else
         {
            printf("%s = \n[ \"%e\" ",name,(*var)[0]);
            for (i=1; i<*size; i++) printf("\n  \"%e\" ",(*var)[i]);
            printf("] (%i elements in list)\n",*size);
         }
      }
      else printf("%s not found in parameter file %s)\n",name,pfile);
   }
}

void AMRD_int_param_v(char *pfile, char *name, int **var, int *size)
{
   int i,found;

   found=get_int_param_v(pfile,name,var,size);
   if (found==-1) AMRD_stop("Error reading parameter file, variable = ",name);

   if (AMRD_echo_params && my_rank==0)
   {
      if (found)
      {
         if (*size==1) printf("%s = [ \"%i\" ]\n",name,(*var)[0]);
         else
         {
            printf("%s = \n[ \"%i\" ",name,(*var)[0]);
            for (i=1; i<*size; i++) printf("\n  \"%i\" ",(*var)[i]);
            printf("] (%i elements in list)\n",*size);
         }
      }
      else printf("%s not found in parameter file %s)\n",name,pfile);
   }
}

//=============================================================================
// calculates the *approximate* (because the overlap regions will be
// double-counted) L2 norm of a grid function over the level
//=============================================================================
real approx_l2norm(int gfn, int L, int hier)
{
   int valid,ngfs,i,l_size;
   real *gfs[PAMR_MAX_GFNS],t;
   int l_rank,l_dim,l_shape[PAMR_MAX_DIM],l_shape_c[PAMR_MAX_DIM],l_ghost_width[2*PAMR_MAX_DIM];
   real *l_x[PAMR_MAX_DIM],*l_x_c[PAMR_MAX_DIM],l_bbox[2*PAMR_MAX_DIM],lnorm[2],gnorm[2];

   valid=PAMR_init_s_iter(L,hier,0);

   lnorm[0]=0;
   lnorm[1]=0;
   while(valid)
   {
      PAMR_get_g_attribs(&l_rank,&l_dim,l_shape,l_shape_c,l_bbox,l_ghost_width,&t,&ngfs,l_x,l_x_c,gfs);
      l_size=1; for (i=0; i<l_dim; i++) l_size*=l_shape[i];
      for (i=0; i<l_size; i++) lnorm[0]+=(gfs[gfn-1])[i]*(gfs[gfn-1])[i];
      lnorm[1]+=l_size;
      valid=PAMR_next_g();
   }

   MPI_Allreduce(&lnorm,&gnorm,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

   return (real)sqrt((double)(gnorm[0]/gnorm[1]));
}

//=============================================================================
// saves a list of bounding boxes to file "tag_L#_my_rank.sdf" ... value of grids
// set to my_rank.
//
// (if lev=0, all levels set to glev)
//=============================================================================
void debug_save_bbox_list(real *bbox, int *lev, int glev, int num, real t, char *tag)
{
   int l,i,n,ltrace=1;

   real data[100];
   int shape[PAMR_MAX_DIM],lvc;
   char buffer[256];

   for (n=0; n<100; n++) data[n]=my_rank;
   for (n=0; n<PAMR_MAX_DIM; n++) shape[n]=2;

   for (n=0; n<num; n++)
   {
      if (lev) lvc=lev[n]; else lvc=glev;
      sprintf(buffer,"%s_L%i_%i",tag,lvc,my_rank);
      gft_out_bbox(buffer,t,shape,AMRD_dim,&bbox[2*AMRD_dim*n],data);
      IFL 
      { 
         if (AMRD_dim==1)
            printf("   debug_save_bbox_list: L=%i, rank=%i, bbox=[%lf,%lf]\n",lvc,my_rank,
                   bbox[2*AMRD_dim*n],bbox[2*AMRD_dim*n+1]);
         else if (AMRD_dim==2)
            printf("   debug_save_bbox_list: L=%i, rank=%i, bbox=[%lf,%lf][%lf,%lf]\n",lvc,my_rank,
                   bbox[2*AMRD_dim*n],bbox[2*AMRD_dim*n+1],bbox[2*AMRD_dim*n+2],bbox[2*AMRD_dim*n+3]);
         else if (AMRD_dim==3)
            printf("   debug_save_bbox_list: L=%i, rank=%i, bbox=[%lf,%lf][%lf,%lf][%lf,%lf]\n",lvc,my_rank,
                   bbox[2*AMRD_dim*n],bbox[2*AMRD_dim*n+1],bbox[2*AMRD_dim*n+2],bbox[2*AMRD_dim*n+3],
                   bbox[2*AMRD_dim*n+4],bbox[2*AMRD_dim*n+5]);
      }
   }
}

//=============================================================================
// prints some statistics for the AMR grid hierarchy
//=============================================================================
void gh_stats(void)
{
   int L,Lf,valid,ngfs,agfs,i,mem,ngs,tmem;
   real *gfs[PAMR_MAX_GFNS];

   if (my_rank!=0) return;

   printf("%s\n",line_break);
   printf("AMR hierarchy statistics ... node %i\n\n",my_rank);

   Lf=PAMR_get_max_lev(PAMR_AMRH);
   L=1;
   valid=0;
   while(!valid && L<=Lf) valid=PAMR_init_s_iter(L++,PAMR_AMRH,0);
   if (!valid) return;
   PAMR_get_g_gfs(gfs);
   PAMR_get_g_ngfs(&ngfs);

   printf("Number of levels: %i\n",Lf);

   agfs=0; for (i=0; i<ngfs; i++) if (gfs[i]) agfs++;

   printf("%i (of %i) grid functions allocated\n\n",agfs,ngfs);

   tmem=0;
   for (L=1; L<=Lf; L++)
   {
      printf("Level %i:\n",L);
      valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
      mem=0;
      ngs=0;
      while(valid)
      {
         ngs++;
         ev_ldptr();
	 // Note this will only be approximate if some of the gfs are
	 // cell centered.
         mem+=AMRD_g_size*agfs*sizeof(real);
         valid=PAMR_next_g();
      }
      tmem+=mem;
      printf("   %i grid(s) using %f MB\n",ngs,((float)mem)/1000000);
   }

   printf("\n Total memory on node %i: %f MB\n",my_rank,((float)tmem)/1000000);

   printf("\n%s\n",line_break);
   return;
}

//=============================================================================
// returns the total memory in megabytes, for this node
//=============================================================================
int total_mem(void)
{
   int L,Lf,valid,agfs,i,mem,tmem,ngfs;
   real *gfs[PAMR_MAX_GFNS];

   if (my_rank!=0) return 0;

   Lf=PAMR_get_max_lev(PAMR_AMRH);
   L=1;
   valid=0;
   while(!valid && L<=Lf) valid=PAMR_init_s_iter(L++,PAMR_AMRH,0);
   if (!valid) return 0;
   PAMR_get_g_gfs(gfs);
   PAMR_get_g_ngfs(&ngfs);

   agfs=0; for (i=0; i<ngfs; i++) if (gfs[i]) agfs++;

   tmem=0;
   for (L=1; L<=Lf; L++)
   {
      valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
      mem=0;
      while(valid)
      {
         ev_ldptr();
	 // Note this will only be approximate if some of the grid
	 // functions are cell centered.
         mem+=AMRD_g_size*agfs*sizeof(real);
         valid=PAMR_next_g();
      }
      tmem+=mem;
   }

   return tmem/1000000;
}

int in_array(int x, int *a, int size)
{
   int i=size;

   while(i--) if (a[i]==x) return 1;

   return 0;
}

//=============================================================================
// a utility function that will call dmrepop3d1_ for all AMRH
// variables at all time levels, together with apropriate synchronization.
// --- repeats the process n-times for repopulation of up to n zones
// on the finest level
//
// NOTE: we only sync after each stage, as PAMR's interpolation routines
// use existing mask and hence will not interpolate to or from newly 
// repopulated regions.
//
// at this stage we only repopulate hyperbolic_vars, elliptic_vars, and
// work variables in the AMRD_work_repop list.
//
// We also dissipate rg_diss_vars if rg_diss_eps>0
//
// SPECIAL FLAG: if def_order <0, then use order=-def_order for *finest* level
//               only, order=1 for rest
//=============================================================================
void AMRD_repopulate(int n, int def_order)
{
   int Lf,L,i,j,valid,first=1,order,k;
   int ltrace=0,even=PAMR_EVEN,odd=PAMR_ODD,pbt[2*PAMR_MAX_DIM],stride;
   real mask_off=AMRD_ex-1;
   real *c_eps;

   Lf=PAMR_get_max_lev(PAMR_AMRH);

   PAMR_set_tf_bits(-1,PAMR_AMRH,PAMR_TF_SYNC);
   if (AMRD_w1_gfn) PAMR_set_tf_bit(AMRD_w1_gfn,PAMR_TF_SYNC);
   if (AMRD_w1_c_gfn) PAMR_set_tf_bit(AMRD_w1_c_gfn,PAMR_TF_SYNC);

   if (ltrace) 
   {
      if (my_rank==0) printf("AMRD_repopulate: n=%i, def_order=%i\n",n,def_order);
      // evo_dump(1,Lf,"before_repop",0,0);
   }

   while(n--)
   {
      for (L=1; L<=Lf; L++)
      {
         valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
         if (ltrace && my_rank==0) printf("AMRD_repopulate: L=%i\n",L);
         while(valid)
         {
            ev_ldptr();
            if (first) {
	      if (AMRD_w1_gfn) for (j=0; j<AMRD_g_size; j++) AMRD_w1[j]=AMRD_chr[j];
	      if (AMRD_w1_c_gfn) for (j=0; j<AMRD_g_size_c; j++) AMRD_w1_c[j]=AMRD_chr_c[j];
	    }
            if (ltrace && my_rank==0) printf("AMRD_repopulate: doing hyperbolics\n");

            for (i=0; i<AMRD_num_hyperbolic_vars; i++)
            {
               if (def_order<0 && L==Lf) order=-def_order;
               else if (def_order<0) order=1;
               else order=def_order;
               if (in_array(AMRD_AMR_f_gfn[0][i],AMRD_ex_repop1_sgfn,AMRD_num_ex_repop1_vars)) order=1;
               if (in_array(AMRD_AMR_f_gfn[0][i],AMRD_ex_repop2_sgfn,AMRD_num_ex_repop2_vars)) order=2;
               if (in_array(AMRD_AMR_f_gfn[0][i],AMRD_ex_repop3_sgfn,AMRD_num_ex_repop3_vars)) order=3;
               if (in_array(AMRD_AMR_f_gfn[0][i],AMRD_ex_repop4_sgfn,AMRD_num_ex_repop4_vars)) order=4;
	       if (PAMR_var_type(AMRD_hyperbolic_vars[i])==PAMR_VERTEX_CENTERED) {
		 for (j=0; j<AMRD_num_evo_tl; j++)
		   dmrepop3d1_(AMRD_AMR_f[j][i],AMRD_w1,&AMRD_ex,&order,&AMRD_g_Nx,&AMRD_g_Ny,&AMRD_g_Nz);
	       } else {
		 order=1;
		 for (j=0; j<AMRD_num_evo_tl; j++)
		   dmrepop3d1_c_(AMRD_AMR_f[j][i],AMRD_w1_c,&AMRD_ex,&order,&AMRD_g_Nx_c,&AMRD_g_Ny_c,&AMRD_g_Nz_c);
	       }
            }

            if (ltrace && my_rank==0) printf("AMRD_repopulate: doing elliptics\n");
            for (i=0; i<AMRD_num_elliptic_vars; i++)
            {
               if (def_order<0 && L==Lf) order=-def_order;
               else if (def_order<0) order=1;
               else order=def_order;
               if (in_array(AMRD_AMR_mgf_gfn[0][i],AMRD_ex_repop1_sgfn,AMRD_num_ex_repop1_vars)) order=1;
               if (in_array(AMRD_AMR_mgf_gfn[0][i],AMRD_ex_repop2_sgfn,AMRD_num_ex_repop2_vars)) order=2;
               if (in_array(AMRD_AMR_mgf_gfn[0][i],AMRD_ex_repop3_sgfn,AMRD_num_ex_repop3_vars)) order=3;
               if (in_array(AMRD_AMR_mgf_gfn[0][i],AMRD_ex_repop4_sgfn,AMRD_num_ex_repop4_vars)) order=4;
	       if (PAMR_var_type(AMRD_elliptic_vars[i])==PAMR_VERTEX_CENTERED) {
		 for (j=0; j<AMRD_num_evo_tl; j++)
		   dmrepop3d1_(AMRD_AMR_mgf[j][i],AMRD_w1,&AMRD_ex,&order,&AMRD_g_Nx,&AMRD_g_Ny,&AMRD_g_Nz);
	       } else {
		 order=1;
		 for (j=0; j<AMRD_num_evo_tl; j++)
		   dmrepop3d1_c_(AMRD_AMR_mgf[j][i],AMRD_w1_c,&AMRD_ex,&order,&AMRD_g_Nx_c,&AMRD_g_Ny_c,&AMRD_g_Nz_c);
	       }
            }
            if (ltrace && my_rank==0) printf("AMRD_repopulate: doing work vars\n");
            //-----------------------------------------------------------------
            // NOTE: The following was added to repopulate a list of chosen 
	    // work variables.  This is likely to be necessary for hydrodynamic
	    // primitive variables, for example.
            //-----------------------------------------------------------------
            for (i=0; i<AMRD_num_work_repop_vars; i++)
            {
               if (def_order<0 && L==Lf) order=-def_order;
               else if (def_order<0) order=1;
               else order=def_order;
               if (in_array(AMRD_work_repop_gfn[i],AMRD_ex_repop1_sgfn,AMRD_num_ex_repop1_vars)) order=1;
               if (in_array(AMRD_work_repop_gfn[i],AMRD_ex_repop2_sgfn,AMRD_num_ex_repop2_vars)) order=2;
               if (in_array(AMRD_work_repop_gfn[i],AMRD_ex_repop3_sgfn,AMRD_num_ex_repop3_vars)) order=3;
               if (in_array(AMRD_work_repop_gfn[i],AMRD_ex_repop4_sgfn,AMRD_num_ex_repop4_vars)) order=4;
	       if (PAMR_var_type(AMRD_work_repop_vars[i])==PAMR_VERTEX_CENTERED) {
		 dmrepop3d1_(AMRD_work_repop[i],AMRD_w1,&AMRD_ex,&order,&AMRD_g_Nx,&AMRD_g_Ny,&AMRD_g_Nz);
	       } else {
		 order=1;
		 dmrepop3d1_c_(AMRD_work_repop[i],AMRD_w1_c,&AMRD_ex,&order,&AMRD_g_Nx_c,&AMRD_g_Ny_c,&AMRD_g_Nz_c);
	       }
            }
            if (ltrace && my_rank==0) printf("AMRD_repopulate: doing w1's\n");
            if (AMRD_w1_gfn) dmrepop3dc1_(AMRD_w1,&AMRD_ex,&AMRD_g_Nx,&AMRD_g_Ny,&AMRD_g_Nz);
	    if (AMRD_w1_c_gfn) dmrepop3dc1_(AMRD_w1_c,&AMRD_ex,&AMRD_g_Nx_c,&AMRD_g_Ny_c,&AMRD_g_Nz_c);
            //-----------------------------------------------------------------
            // NOTE: using the cmask as a temporary variable
            //-----------------------------------------------------------------
            if (ltrace && my_rank==0) printf("AMRD_repopulate: AMRD_rg_eps_diss=%lf, AMRD_num_rg_diss_vars=%i\n",AMRD_rg_eps_diss,AMRD_num_rg_diss_vars);
            if (AMRD_rg_eps_diss!=0 && AMRD_num_rg_diss_vars>0)
            {
               for (j=0; j<AMRD_tnum_rg_diss_vars;j++)
               {
                  for (k=0; k<2*AMRD_dim; k++) if (AMRD_g_phys_bdy[k]) pbt[k]=AMRD_rg_diss_pbt[j][k]; else pbt[k]=PAMR_UNKNOWN;
                  c_eps=&AMRD_rg_eps_diss; if (AMRD_do_ex<0) c_eps=AMRD_w1;
                  for (stride=AMRD_diss_stride; stride>0; stride--)
                     dmdiss3d_(AMRD_rg_diss[j],AMRD_cmask,c_eps,&AMRD_repop_diss_bdy,
                            pbt,&even,&odd,
                            AMRD_w1,&mask_off,&AMRD_g_Nx,&AMRD_g_Ny,&AMRD_g_Nz,
                            AMRD_w1,&AMRD_ex,&AMRD_do_ex,&AMRD_diss_use_6th_order,&stride);
               }
            }
            valid=PAMR_next_g();
         }
         PAMR_sync(L,-1,PAMR_AMRH,0);
      }
      first=0;
   }

   PAMR_thaw_tf_bits();

   // if (ltrace) evo_dump(1,Lf,"after_repop",0,0);
   if (ltrace && my_rank==0) printf("AMRD_repopulate: done\n");

}


// Note that "name" is the output from this function.
void AMRD_append_tag(char *base_name, char *tag, char *v_tag, char *c_tag)
{
  int i, lbase, lctag, lvtag, ltag;
  int check_vtag, check_ctag;
  int vtag_found, ctag_found;
  char *tmpchar;

  check_vtag = 1;
  check_ctag = 1;
  vtag_found = 0;
  ctag_found = 0;

  lbase=0;
  lctag=0;
  lvtag=0;

  if (base_name) lbase=strlen(base_name);
  if (c_tag) lctag=strlen(c_tag);
  if (v_tag) lvtag=strlen(v_tag);
  ltag=strlen(tag);

  if (lctag == 0 || lctag >= lbase) check_ctag = 0;
  if (lvtag == 0 || lvtag >= lbase) check_vtag = 0;

  char tmpn[lbase+lctag+lvtag+ltag];

  // Check for c_tag
  if (check_ctag) {
    char tmpc[lctag+1];
    
    for (i=0; i<lctag; i++) {
      tmpc[i] = base_name[lbase-lctag+i];
    }
    tmpc[lctag]= '\0';

    if (!(strcmp(tmpc,c_tag))) ctag_found=1;

    // remove ctag from temporary name
    for (i=0; i<(lbase-lctag); i++) tmpn[i] = base_name[i];
    tmpn[lbase-lctag] = '\0';
  }

  // Check for v_tag
  if (check_vtag) {
    char tmpv[lvtag+1];

    for (i=0; i<lvtag; i++) {
      tmpv[i] = base_name[lbase-lvtag+i];
    }
    tmpv[lvtag] = '\0';

    if (!(strcmp(tmpv,v_tag))) vtag_found=1;

    // remove vtag from temporary name
    for (i=0; i<(lbase-lvtag); i++) tmpn[i] = base_name[i];
    tmpn[lbase-lvtag]='\0';
  }

  if (ctag_found || vtag_found) {
    // Append the desired tag first, then append the v or c tag.
    tmpchar = strcat(tmpn,tag);
    if (vtag_found) {
      strcpy(base_name,strcat(tmpchar,v_tag));
    } else {
      strcpy(base_name,strcat(tmpchar,c_tag));
    }
  } else {
    // Just append as normal
    base_name = strcat(base_name,tag);
  }

  return;

}

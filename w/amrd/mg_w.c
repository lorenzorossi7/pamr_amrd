//=============================================================================
// mg.c --- multigrid routines
//=============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "globals_w.h"
#include "mg_w.h"
#include "mpi.h"
#include "util_w.h"
#include "math.h"
#include "evolve_w.h"

//=============================================================================
// solves the set of elliptic equations from levels L..Lf at time level tl.
//
// if save_tre, then l2norm of tre is saved to mg_tre (which is f_tre[0])
//=============================================================================
void sync_elliptic_vars_t0_data(int L);
void inject_elliptic_vars_t0_data(int Lf);
void interp_elliptic_vars_t0_data(int Lc);

void solve_elliptics(int L, int tl, int save_tre)
{
   int l,Lf,Lf_mg;
   real res=AMRD_MG_tol+1;
   real t1,tLf;

   Lf=PAMR_get_max_lev(PAMR_AMRH);

   if ( (AMRD_num_elliptic_vars==0 && !AMRD_is_t0) || 
        (AMRD_num_elliptic_vars==0 && AMRD_num_elliptic_vars_t0==0 && AMRD_is_t0) ) return;

   if (!(PAMR_build_mgh(L,Lf,tl))) AMRD_stop("solve_elliptics: PAMR_build_mgh() failed\n",0);
   
   // If the TRE gfns are cell centered, then it doesn't make sense to try 
   // to calculate them in this MG setting.  And at present, if the first 
   // TRE gnf is cell centered, all of them must be.
   if (PAMR_gfn_var_type(AMRD_mg_tre_gfn)==PAMR_CELL_CENTERED) save_tre=0;

   Lf_mg=PAMR_get_max_lev(PAMR_MGH);

   if (AMRD_MG_trace && my_rank==0) 
   {
      printf("\n\n");
      for (l=0; l<(L-1); l++) printf(" ");
      printf("MG L = %i - %i, MG Lf=%i, tl=%i ...",L,Lf,Lf_mg,tl);
   }
   if (AMRD_MG_trace>1 && my_rank==0) printf("\n");

   Lf=PAMR_get_max_lev(PAMR_MGH);

   if (app_elliptic_vars_t0_init && AMRD_is_t0)
   {
      for (l=1; l<=Lf; l++)
      {
         call_app(app_elliptic_vars_t0_init,l,PAMR_MGH);
         sync_elliptic_vars_t0_data(l);
         if (l<Lf) interp_elliptic_vars_t0_data(l); 
      }
      for (l=Lf; l>1; l--) inject_elliptic_vars_t0_data(l);
   }

   AMRD_curr_MG_iter=0;
   while ((AMRD_curr_MG_iter < AMRD_MG_min_iter || res > AMRD_MG_tol) && AMRD_curr_MG_iter < AMRD_MG_max_iter)
   {
      res=vcycle(AMRD_curr_MG_iter,save_tre);
      if (AMRD_num_MG_cnst_data_vars>0 && AMRD_is_t0)
      {
         for (l=1; l<=Lf; l++) 
         {
            call_app(app_t0_cnst_data,l,PAMR_MGH);
            sync_cnst_data(l);
            if (l<Lf) interp_cnst_data(l); 
         }
         for (l=Lf; l>1; l--) inject_cnst_data(l);
      }
      AMRD_curr_MG_iter++;
      if (AMRD_MG_trace>1 && my_rank==0) printf("iter %i\t\tres(Lf)=%12.8e\n",AMRD_curr_MG_iter,res);
   }

   if (AMRD_curr_MG_iter==AMRD_MG_max_iter && my_rank==0 && res > AMRD_MG_tol) 
      printf("\nWARNING ... failed to solve elliptics to within %12.8e in %i iterations; res=%12.8e\n",AMRD_MG_tol,AMRD_MG_max_iter,res);
   else if (AMRD_MG_trace && my_rank==0) printf("  iters=%i\n",AMRD_curr_MG_iter);

   PAMR_destroy_mgh();
}

//=============================================================================
// 1 vcycle ... returns the residual the finest level
//=============================================================================
real vcycle(int giter, int save_tre)
{
   int L,Lf,iter,i,cAMR_bdy_width;
   real res[PAMR_MAX_LEVS],ctol;
   static int cwarn_count=10;

   Lf=PAMR_get_max_lev(PAMR_MGH);
   clear_MG_vars(Lf);

   //--------------------------------------------------------------------------
   // down the V
   //--------------------------------------------------------------------------
   for (L=Lf; L>0; L--)
   {
      if (L<Lf)
      {
         //--------------------------------------------------------------------
         // compute the residual on L+1, and synchronize (necessary for
         // convergence of MG)
         //--------------------------------------------------------------------
         res[L+1]=call_rr_app(app_MG_residual,L+1,PAMR_MGH,0,MPI_MAX);
         set_res_sync();
         PAMR_sync(L+1,0,PAMR_MGH,0); PAMR_thaw_tf_bits();
         // mg_dump(L+1,"after_calc_residual_",giter,1);

         //--------------------------------------------------------------------
         // inject grid functions (including residual) from L+1 to L
         //--------------------------------------------------------------------
         clear_MG_vars(L);
         PAMR_inject(L+1,0,PAMR_MGH);
         //--------------------------------------------------------------------
         // need to sync here, so that copy of f saved in MG_rhs is consistent
         // for subsequent cgc calculations ... essential for MG 
         PAMR_sync(L,0,PAMR_MGH,0);

         //--------------------------------------------------------------------
         // compute new RHS in regions where interior child grids are.
         // note that compute_MG_rhs() saves f(L) in f_lop(L)
         //--------------------------------------------------------------------
         set_cmask_child(L,PAMR_MGH);
         call_app(app_L_op,L,PAMR_MGH);  
         mg_dump(L,"before_compute_rhs_",giter,0);
         compute_MG_rhs(L,save_tre);
      }

      //-----------------------------------------------------------------------
      // set cmask to ON at interior and physical boundary points, and sweep.
      //-----------------------------------------------------------------------
      set_cmask_bdy_w1(L,PAMR_MGH);
      mg_dump(L,"before_preswp_",giter,0);

      for (iter=0; iter < AMRD_MG_pre_swp; iter++)
      {
         res[L]=call_rr_app(app_MG_relax,L,PAMR_MGH,0,MPI_MAX);
         if (AMRD_MG_trace>3  && my_rank==0) { for (i=L-1; i<Lf; i++) printf("   "); printf("(L%i(%i):%12.8e)\n",L,iter,res[L]); }
         if (iter==0) cAMR_bdy_width=0; else cAMR_bdy_width=1;
         PAMR_sync(L,0,PAMR_MGH,cAMR_bdy_width);
      }
      if (AMRD_MG_trace>2 && my_rank==0) { for (i=L-1; i<Lf; i++) printf("   "); printf("L%i(%i):%12.8e\n",L,iter,res[L]); }

      mg_dump(L,"after_preswp_",giter,0);

   }

   //--------------------------------------------------------------------------
   // up the V
   //--------------------------------------------------------------------------
   for (L=1; L<=Lf; L++)
   {
      ctol=res[L]*AMRD_MG_crtol;

      //---------------------------------------------------------------------------
      // apply CGC. compute_cgc uses the prior f saved in f_lop by compute_MG_rhs()
      //---------------------------------------------------------------------------
      if (L>1) compute_cgc(L);

      // mg_dump(L,"after_cgc_",giter,1);

      for (iter=0; iter<AMRD_MG_pst_swp; iter++)
      {
         res[L]=call_rr_app(app_MG_relax,L,PAMR_MGH,0,MPI_MAX);
         if (AMRD_MG_trace>3  && my_rank==0) { for (i=L-1; i<Lf; i++) printf("   "); printf("(L%i(%i):%12.8e)\n",L,iter,res[L]); }
         if (iter==0) cAMR_bdy_width=0; else cAMR_bdy_width=1;
         PAMR_sync(L,0,PAMR_MGH,cAMR_bdy_width);
      }

      //-----------------------------------------------------------------------
      // if coarsest grids exist here, sweep over them
      //-----------------------------------------------------------------------
      if (is_coarsest(L))
      {
         res[L]=ctol+1;
         set_comm_coarsest(L);
         while(iter < AMRD_MG_max_citer && res[L] > ctol)
         {
            res[L]=call_rr_app(app_MG_relax,L,PAMR_MGH,1,MPI_MAX);
            PAMR_sync(L,0,PAMR_MGH,1);
            iter++;
         }
         if (iter==AMRD_MG_max_citer && my_rank==0 && cwarn_count)
         {
            printf("WARNING ... failed to solve elliptics on coarsest grid to within"
                   " %12.8e in %i iterations; res=%12.8e\n",ctol,AMRD_MG_max_citer,res[L]);
            cwarn_count--;
            if (cwarn_count) printf("(only issuing %i more of these warnings)\n",cwarn_count);
         }
         reset_comm(L);
      }
      //-----------------------------------------------------------------------
      // return the accurate, final level residual residual
      //-----------------------------------------------------------------------
      if (L==Lf) res[Lf]=call_rr_app(app_MG_residual,L,PAMR_MGH,0,MPI_MAX);
      if (AMRD_MG_trace>2  && my_rank==0) { for (i=L-1; i<Lf; i++) printf("   "); printf("L%i(%i):%12.8e\n",L,iter,res[L]); }

      mg_dump(L,"after_pstswp_",giter,1);
   }

   return res[Lf];
}

//=============================================================================
// internal routines
//=============================================================================

//=============================================================================
// fills the data pointers for current grid (io.c initialized gfn's)
//=============================================================================
void mg_ldptr()
{
   int ngfs,i;
   real *gfs[PAMR_MAX_GFNS],t;

   AMRD_g_shape[0]=AMRD_g_shape[1]=AMRD_g_shape[2]=1;
   AMRD_g_shape_c[0]=AMRD_g_shape_c[1]=AMRD_g_shape_c[2]=1;
   if (!(PAMR_get_g_attribs(&AMRD_g_rank,&AMRD_g_dim,AMRD_g_shape,AMRD_g_shape_c,AMRD_g_bbox,AMRD_g_ghost_width,&t,&ngfs,AMRD_g_x,AMRD_g_x_c,gfs)))
      AMRD_stop("mg_ldptr: PAMR_get_g_attribs failed\n","");

   AMRD_g_size=1; 
   for (i=0; i<AMRD_g_dim; i++) { AMRD_g_size*=AMRD_g_shape[i]; AMRD_g_dx[i]=(AMRD_g_x[i])[2]-(AMRD_g_x[i])[1]; }

   AMRD_cmask_mg=gfs[AMRD_cmask_mg_gfn-1];
   for (i=0; i<(AMRD_num_elliptic_vars+AMRD_num_elliptic_vars_t0); i++)
   {
      AMRD_MG_res[i]=gfs[AMRD_MG_res_gfn[i]-1];
      AMRD_MG_lop[i]=gfs[AMRD_MG_lop_gfn[i]-1];
      AMRD_MG_rhs[i]=gfs[AMRD_MG_rhs_gfn[i]-1];
      if (i<AMRD_num_elliptic_vars) AMRD_MG_brs[i]=gfs[AMRD_MG_brs_gfn[i]-1]; // this var is actually in the AMR hierarchy
      AMRD_MG_f[i]=gfs[AMRD_MG_f_gfn[i]-1];
   }

   for (i=0; i<AMRD_num_MG_cnst_data_vars; i++) AMRD_MG_cnst_data[i]=gfs[AMRD_MG_cnst_data_gfn[i]-1];

   if (AMRD_mg_tre_gfn) AMRD_mg_tre=gfs[AMRD_mg_tre_gfn-1];

   if (AMRD_do_ex) 
   {
      if (AMRD_chr_mg_gfn) AMRD_chr_mg=gfs[AMRD_chr_mg_gfn-1]; else AMRD_chr_mg=AMRD_cmask_mg;
      if (AMRD_chr_mg_c_gfn) AMRD_chr_mg_c=gfs[AMRD_chr_mg_c_gfn-1]; 
   }
}

//=============================================================================
// sets the MG variables to 0
//=============================================================================
void clear_MG_vars(int L)
{
   int valid,i,j,num; 

   valid=PAMR_init_s_iter(L,PAMR_MGH,0);

   num=AMRD_num_elliptic_vars; if (AMRD_is_t0) num+=AMRD_num_elliptic_vars_t0;

   while(valid)
   {
      mg_ldptr();
      for (i=0; i<num; i++)
      {
         for (j=0; j<AMRD_g_size; j++) 
         {
            (AMRD_MG_res[i])[j]=(AMRD_MG_rhs[i])[j]=(AMRD_MG_lop[i])[j]=0;
         }
      }
      valid=PAMR_next_g();
   }
}

//=============================================================================
// computes the new rhs vectors as f_rhs(L) = f_res(L+1)*MG_w0 + f_lop(L)
// AND afterward saves a copy of f, prior to pre-sweeps, in f_lop.
//
// If (save_tre), then we save the infinity norm of the MG_lop in mg_tre,
// and interpolate it to the parent level.
//=============================================================================
void compute_MG_rhs(int L, int save_tre)
{
   int valid,i,j,num; 

   valid=PAMR_init_s_iter(L,PAMR_MGH,0);

   if (save_tre && !(AMRD_mg_tre_gfn)) AMRD_stop("compute_MG_rhs: error ... save_tre=1, but no tre variables\n",0);

   num=AMRD_num_elliptic_vars; if (AMRD_is_t0) num+=AMRD_num_elliptic_vars_t0;

   while(valid)
   {
      mg_ldptr();

      for (i=0; i<num; i++)
      {
         for (j=0; j<AMRD_g_size; j++) 
         {
            if (save_tre) 
            {
               if (i==0) AMRD_mg_tre[j]=0;
               //-------------------------------------------------------------------
               // multiply by sqrt(AMRD_num_f_tre_vars), as calc_tre will divide by this
               //-------------------------------------------------------------------
               AMRD_mg_tre[j]=max(AMRD_mg_tre[j],fabs((AMRD_MG_lop[i])[j])*sqrt(AMRD_num_f_tre_vars));
            }
            if ((AMRD_cmask_mg)[j]==AMRD_CMASK_ON) 
               (AMRD_MG_rhs[i])[j]=-(AMRD_MG_res[i])[j]*AMRD_MG_w0_r+(AMRD_MG_lop[i])[j];
            else 
               (AMRD_MG_rhs[i])[j]=0;
            (AMRD_MG_lop[i])[j]=(AMRD_MG_f[i])[j];
         }
      }
      valid=PAMR_next_g();
   }

   if (save_tre)
   {
      set_gfn_sync(AMRD_mg_tre_gfn);
      PAMR_sync(L,0,PAMR_MGH,0); 
      set_gfn_in(AMRD_mg_tre_gfn,PAMR_SECOND_ORDER);
      PAMR_interp(L,0,PAMR_MGH); 
      PAMR_thaw_tf_bits();
   }
}

//=============================================================================
// computes the coarse-grid correction to MG variables f at level L (>1)
//
// NOTE: this function re-uses f_lop var's, AND assumes that f_lop has
//       been set to the value of f prior to relaxation, so that we don't 
//       need to re-restrict f for computing the CGC
//=============================================================================
void compute_cgc(int L)
{
   int valid,i,j,num; 

   num=AMRD_num_elliptic_vars; if (AMRD_is_t0) num+=AMRD_num_elliptic_vars_t0;
   //--------------------------------------------------------------------------
   // only f_res's should be defined with an active mg_interp ... so
   // copy correction to f_res on L-1, then interpolate.
   //--------------------------------------------------------------------------

   valid=PAMR_init_s_iter(L-1,PAMR_MGH,0);
   while(valid)
   {
      mg_ldptr();
      for (i=0; i<num; i++) 
      {
         for (j=0; j<AMRD_g_size; j++) (AMRD_MG_res[i])[j]=((AMRD_MG_f[i])[j]-(AMRD_MG_lop[i])[j])/AMRD_MG_w0_i;
      }
      valid=PAMR_next_g();
   }

   PAMR_interp(L-1,0,PAMR_MGH); 

   valid=PAMR_init_s_iter(L,PAMR_MGH,0);
   while(valid)
   {
      mg_ldptr();
      for (i=0; i<num; i++) 
      {
         for (j=0; j<AMRD_g_size; j++) {
	   (AMRD_MG_f[i])[j]+=(AMRD_MG_res[i])[j];
	 }
      }
      valid=PAMR_next_g();
   }


   if (AMRD_MG_reinterp_bdy) { 
      /*
      Interpolate the values on the AMR boundaries for points 
      that don't exist on the parent grid from those that do. 
      */
      int rho_sp[PAMR_MAX_LEVS],rho_tm[PAMR_MAX_LEVS]; //we need rho's to use interp_AMR_bdy
      PAMR_get_rho(rho_sp,rho_tm,PAMR_MAX_LEVS);  //write rho into the data holders

      valid=PAMR_init_s_iter(L,PAMR_MGH,0);
      while (valid) {
         ev_ldptr();
         mg_ldptr();
         for (i=0; i<num; i++) {
            //Note that we are using res as a work function here
            interp_AMR_bdy(AMRD_MG_f[i], AMRD_MG_res[i], rho_sp[L-2]);  //see evolve_w.c
         }
         valid=PAMR_next_g();
      }
      PAMR_sync(L,0,PAMR_MGH,0);
   }
}

//=============================================================================
// sets the mg_cmask (if MGH) or cmask (if AMRH) of local grids to 1 in the 
// region (-1 in from interior boundaries) where (any) child grids are, 
// 0 elsewhere
//=============================================================================
void set_cmask_child(int L, int hier)
{
   int valid,is,js,ks,ie,je,ke,i,j,k;
   real cbbox[2*PAMR_MAX_DIM];
   real *mask;

   valid=PAMR_init_s_iter(L,hier,0);

   while(valid)
   {
      if (hier==PAMR_MGH) 
      {
         mg_ldptr(); 
         mask=AMRD_cmask_mg;
      }
      else 
      {
         ev_ldptr();
         mask=AMRD_cmask;
      }
      for (j=0; j<AMRD_g_size; j++) mask[j]=AMRD_CMASK_OFF;
      PAMR_push_iter();
      valid=PAMR_init_s_iter(L+1,hier,1);
      while(valid)
      {
         PAMR_get_g_bbox(cbbox);
         js=je=ks=ke=0;
         is=max((cbbox[0]-AMRD_g_bbox[0])/AMRD_g_dx[0]+0.5,0);
         ie=min((cbbox[1]-AMRD_g_bbox[0])/AMRD_g_dx[0]+0.5,AMRD_g_shape[0]-1);
         // if child extends beyond periodic max, the above test will miss it if ghost_width>0 ( for all dims)
         if (AMRD_periodic[0] && cbbox[1]>=AMRD_base_bbox[1]) ie=AMRD_g_shape[0]-1; 
         if (is>0 || AMRD_periodic[0] || (AMRD_g_bbox[0]-AMRD_base_bbox[0])>AMRD_g_dx[0]/2) is++;
         if (ie<(AMRD_g_shape[0]-1) || AMRD_periodic[0] || (AMRD_base_bbox[1]-AMRD_g_bbox[1])>AMRD_g_dx[0]/2) ie--;
         if (AMRD_g_dim>1)
         {
            js=max((cbbox[2]-AMRD_g_bbox[2])/AMRD_g_dx[1]+0.5,0);
            je=min((cbbox[3]-AMRD_g_bbox[2])/AMRD_g_dx[1]+0.5,AMRD_g_shape[1]-1);
            if (AMRD_periodic[1] && cbbox[3]>=AMRD_base_bbox[3]) je=AMRD_g_shape[1]-1; 
            if (js>0 || AMRD_periodic[1] || (AMRD_g_bbox[2]-AMRD_base_bbox[2])>AMRD_g_dx[1]/2) js++;
            if (je<(AMRD_g_shape[1]-1) || AMRD_periodic[1] || (AMRD_base_bbox[3]-AMRD_g_bbox[3])>AMRD_g_dx[1]/2) je--;
            if (AMRD_g_dim>2)
            {
               ks=max((cbbox[4]-AMRD_g_bbox[4])/AMRD_g_dx[2]+0.5,0);
               ke=min((cbbox[5]-AMRD_g_bbox[4])/AMRD_g_dx[2]+0.5,AMRD_g_shape[2]-1);
               if (AMRD_periodic[2] && cbbox[5]>=AMRD_base_bbox[5]) ke=AMRD_g_shape[2]-1; 
               if (ks>0 || AMRD_periodic[2] || (AMRD_g_bbox[4]-AMRD_base_bbox[4])>AMRD_g_dx[2]/2) ks++;
               if (ke<(AMRD_g_shape[2]-1) || AMRD_periodic[2] || (AMRD_base_bbox[5]-AMRD_g_bbox[5])>AMRD_g_dx[2]/2) ke--;
            }
         }
         if (is>=0 && js>=0 && ks>=0 && is<=ie && js<=je && ks<=ke &&
             ie<AMRD_g_shape[0] && je<AMRD_g_shape[1] && ke<AMRD_g_shape[2])
         {
            for (i=is; i<=ie; i++)
            {
               for (j=js; j<=je; j++)
               {
                  for (k=ks; k<=ke; k++)
                  {
                     mask[i+j*AMRD_g_shape[0]+k*AMRD_g_shape[0]*AMRD_g_shape[1]]=AMRD_CMASK_ON;
                  }
               }
            }
         }
         valid=PAMR_next_g();
      }
      PAMR_pop_iter(); 
      valid=PAMR_next_g();
   }
}

//=============================================================================
// sets the mg_cmask(if MGH) or cmask (if AMRH) of local grids to ON in 
// the interior of the grid (excluding the excised region)
//=============================================================================
void set_cmask_bdy(int L,int hier)
{
   int valid,is,js,ks,ie,je,ke,i,j,k,ind; 
   int is_c,js_c,ks_c,ie_c,je_c,ke_c;
   real cbbox[2*PAMR_MAX_DIM];
   real *mask,*ex_mask;
   real *mask_c,*ex_mask_c;

   valid=PAMR_init_s_iter(L,hier,0);

   while(valid)
   {
      if (hier==PAMR_MGH) 
      {
         mg_ldptr(); 
         mask=AMRD_cmask_mg;
         ex_mask=AMRD_chr_mg;
         ex_mask_c=AMRD_chr_mg_c;
      }
      else 
      {
         ev_ldptr();
         mask=AMRD_cmask;
         mask_c=AMRD_cmask_c;
         ex_mask=AMRD_chr;
         ex_mask_c=AMRD_chr_c;
      }
      for (j=0; j<AMRD_g_size; j++) mask[j]=AMRD_CMASK_OFF;
      if (hier==PAMR_AMRH) {
	for (j=0; j<AMRD_g_size_c; j++) mask_c[j]=AMRD_CMASK_OFF;
      }
      is=js=ks=je=ke=0;
      is_c=js_c=ks_c=je_c=ke_c=0;
      ie=AMRD_g_shape[0]-1;
      ie_c=AMRD_g_shape_c[0]-1;
      if ((AMRD_g_bbox[0]-AMRD_base_bbox[0])>AMRD_g_dx[0]/2 || AMRD_periodic[0]) {
	is += AMRD_AMR_bdy_width;
	is_c += AMRD_AMR_bdy_width_c;
      }
      if ((AMRD_base_bbox[1]-AMRD_g_bbox[1])>AMRD_g_dx[0]/2 || AMRD_periodic[0]) {
	ie -= AMRD_AMR_bdy_width;
	ie_c -= AMRD_AMR_bdy_width_c;
      }
      if (AMRD_g_dim>1)
      {
         je=AMRD_g_shape[1]-1;
	 je_c=AMRD_g_shape_c[1]-1;
         if ((AMRD_g_bbox[2]-AMRD_base_bbox[2])>AMRD_g_dx[1]/2 || AMRD_periodic[1]) {
	   js += AMRD_AMR_bdy_width;
	   js_c += AMRD_AMR_bdy_width_c;
	 }
         if ((AMRD_base_bbox[3]-AMRD_g_bbox[3])>AMRD_g_dx[1]/2 || AMRD_periodic[1]) {
	   je -= AMRD_AMR_bdy_width;
	   je_c -= AMRD_AMR_bdy_width_c;
	 }
         if (AMRD_g_dim>2)
         {
            ke=AMRD_g_shape[2]-1;
	    ke_c=AMRD_g_shape_c[2]-1;
            if ((AMRD_g_bbox[4]-AMRD_base_bbox[4])>AMRD_g_dx[2]/2 || AMRD_periodic[2]) {
	      ks += AMRD_AMR_bdy_width;
	      ks_c += AMRD_AMR_bdy_width_c;
	    }
            if ((AMRD_base_bbox[5]-AMRD_g_bbox[5])>AMRD_g_dx[2]/2 || AMRD_periodic[2]) {
	      ke -= AMRD_AMR_bdy_width;
	      ke_c -= AMRD_AMR_bdy_width_c;
	    }
         }
      }

      if (is>=0 && js>=0 && ks>=0 && is<=ie && js<=je && ks<=ke &&
          ie<AMRD_g_shape[0] && je<AMRD_g_shape[1] && ke<AMRD_g_shape[2])
      {
         for (i=is; i<=ie; i++)
         {
            for (j=js; j<=je; j++)
            {
               for (k=ks; k<=ke; k++)
               {
                  ind=i+j*AMRD_g_shape[0]+k*AMRD_g_shape[0]*AMRD_g_shape[1];
                  mask[ind]=AMRD_CMASK_ON;
               }
            }
         }
      }

      // Deal with the cell-centered mask, but only if this is the AMR hierarchy.
      if (hier==PAMR_AMRH && is_c>=0 && js_c>=0 && ks_c>=0 && is_c<=ie_c && js_c<=je_c && ks_c<=ke_c &&
          ie_c<AMRD_g_shape_c[0] && je_c<AMRD_g_shape_c[1] && ke_c<AMRD_g_shape_c[2]) {
	for (i=is_c; i<=ie_c; i++) {
	  for (j=js_c; j<=je_c; j++) {
	    for (k=ks_c; k<=ke_c; k++) {
	      ind=i+j*AMRD_g_shape_c[0]+k*AMRD_g_shape_c[0]*AMRD_g_shape_c[1];
	      mask_c[ind]=AMRD_CMASK_ON;
	    }
	  }
	}
      }
      
      valid=PAMR_next_g();
   }
}


//=============================================================================
// sets the mg_cmask(if MGH) or cmask (if AMRH) of local grids to ON in 
// the interior of the grid (excluding the excised region)
// this forces the boundary width equal to one vertex, zero cells.
// shouldn't matter for cells because this set_cmask should only be used
// for multigrid stuff.
//=============================================================================
void set_cmask_bdy_w1(int L,int hier)
{
   int valid,is,js,ks,ie,je,ke,i,j,k,ind; 
   int is_c,js_c,ks_c,ie_c,je_c,ke_c;
   real cbbox[2*PAMR_MAX_DIM];
   real *mask,*ex_mask;
   real *mask_c,*ex_mask_c;

   valid=PAMR_init_s_iter(L,hier,0);

   while(valid)
   {
      if (hier==PAMR_MGH) 
      {
         mg_ldptr(); 
         mask=AMRD_cmask_mg;
         ex_mask=AMRD_chr_mg;
         ex_mask_c=AMRD_chr_mg_c;
      }
      else 
      {
         ev_ldptr();
         mask=AMRD_cmask;
         mask_c=AMRD_cmask_c;
         ex_mask=AMRD_chr;
         ex_mask_c=AMRD_chr_c;
      }
      for (j=0; j<AMRD_g_size; j++) mask[j]=AMRD_CMASK_OFF;
      if (hier==PAMR_AMRH) {
	for (j=0; j<AMRD_g_size_c; j++) mask_c[j]=AMRD_CMASK_OFF;
      }
      is=js=ks=je=ke=0;
      is_c=js_c=ks_c=je_c=ke_c=0;
      ie=AMRD_g_shape[0]-1;
      ie_c=AMRD_g_shape_c[0]-1;
      if ((AMRD_g_bbox[0]-AMRD_base_bbox[0])>AMRD_g_dx[0]/2 || AMRD_periodic[0]) is += 1;
      if ((AMRD_base_bbox[1]-AMRD_g_bbox[1])>AMRD_g_dx[0]/2 || AMRD_periodic[0]) ie -= 1;
      if (AMRD_g_dim>1)
      {
         je=AMRD_g_shape[1]-1;
	 je_c=AMRD_g_shape_c[1]-1;
         if ((AMRD_g_bbox[2]-AMRD_base_bbox[2])>AMRD_g_dx[1]/2 || AMRD_periodic[1]) js += 1;
         if ((AMRD_base_bbox[3]-AMRD_g_bbox[3])>AMRD_g_dx[1]/2 || AMRD_periodic[1]) je -= 1;
         if (AMRD_g_dim>2)
         {
            ke=AMRD_g_shape[2]-1;
	    ke_c=AMRD_g_shape_c[2]-1;
            if ((AMRD_g_bbox[4]-AMRD_base_bbox[4])>AMRD_g_dx[2]/2 || AMRD_periodic[2]) ks += 1;
            if ((AMRD_base_bbox[5]-AMRD_g_bbox[5])>AMRD_g_dx[2]/2 || AMRD_periodic[2]) ke -= 1;
         }
      }

      if (is>=0 && js>=0 && ks>=0 && is<=ie && js<=je && ks<=ke &&
          ie<AMRD_g_shape[0] && je<AMRD_g_shape[1] && ke<AMRD_g_shape[2])
      {
         for (i=is; i<=ie; i++)
         {
            for (j=js; j<=je; j++)
            {
               for (k=ks; k<=ke; k++)
               {
                  ind=i+j*AMRD_g_shape[0]+k*AMRD_g_shape[0]*AMRD_g_shape[1];
                  mask[ind]=AMRD_CMASK_ON;
               }
            }
         }
      }

      // Deal with the cell-centered mask, but only if this is the AMR hierarchy.
      if (hier==PAMR_AMRH && is_c>=0 && js_c>=0 && ks_c>=0 && is_c<=ie_c && js_c<=je_c && ks_c<=ke_c &&
          ie_c<AMRD_g_shape_c[0] && je_c<AMRD_g_shape_c[1] && ke_c<AMRD_g_shape_c[2]) {
	for (i=is_c; i<=ie_c; i++) {
	  for (j=js_c; j<=je_c; j++) {
	    for (k=ks_c; k<=ke_c; k++) {
	      ind=i+j*AMRD_g_shape_c[0]+k*AMRD_g_shape_c[0]*AMRD_g_shape_c[1];
	      mask_c[ind]=AMRD_CMASK_ON;
	    }
	  }
	}
      }
      
      valid=PAMR_next_g();
   }
}

//=============================================================================
// this function extends the AMRD_CMASK_OFF character of boundary zones by
// num points from edges ...this is so that after a sync, local grids can
// know that a boundary with AMRD_CMASK_ON is at least num points interior
// (only works as long as num < minimum grid width)
//=============================================================================
#define AMRD_CMASK_TMP_OFF -100
void extend_cmask_bdy(int L, int hier, int num)
{
   int valid,ip,jp,kp,im,jm,km,i,j,k; 
   real *mask;

   valid=PAMR_init_s_iter(L,hier,0);

   while(valid)
   {
      if (hier==PAMR_MGH) 
      {
         mg_ldptr(); 
         mask=AMRD_cmask_mg;
      }
      else 
      {
         ev_ldptr();
         mask=AMRD_cmask;
      }
      for (i=0; i<AMRD_g_Nx; i++)
      {
         ip=min(AMRD_g_Nx-1,i+num);
         im=max(0,i-num);
         for (j=0; j<AMRD_g_Ny; j++)
         {
            jp=min(AMRD_g_Ny-1,j+num);
            jm=max(0,j-num);
            for (k=0; k<AMRD_g_Nz; k++)
            {
               kp=min(AMRD_g_Nz-1,k+num);
               km=max(0,k-num);
               if (mask[i+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_ON)
               {
                  if ( ((j==0 || j==(AMRD_g_Ny-1) || k==0 || k==(AMRD_g_Nz-1)) && 
                        (mask[ip+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF || 
                         mask[im+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF)) ||
                       ((i==0 || i==(AMRD_g_Nx-1) || k==0 || k==(AMRD_g_Nz-1)) && 
                        (mask[i+jp*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF || 
                         mask[i+jm*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF)) ||
                       ((i==0 || i==(AMRD_g_Nx-1) || j==0 || j==(AMRD_g_Ny-1)) && 
                        (mask[i+j*AMRD_g_Nx+kp*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF || 
                         mask[i+j*AMRD_g_Nx+km*AMRD_g_Nx*AMRD_g_Ny]==AMRD_CMASK_OFF)) )
                     mask[i+j*AMRD_g_Nx+k*AMRD_g_Nx*AMRD_g_Ny]=AMRD_CMASK_TMP_OFF;
               }
            }
         }
      }
               
      for (i=0; i<AMRD_g_size; i++) if (mask[i]==AMRD_CMASK_TMP_OFF) mask[i]=AMRD_CMASK_OFF;
      valid=PAMR_next_g();
   }
}

//=============================================================================
// returns 1 if any grid in level L is coarsest
//=============================================================================
int is_coarsest(int L)
{
   int valid,coarsest;

   valid=PAMR_init_s_iter(L,PAMR_MGH,1); // need to loop over all, for sync'ing

   while(valid)
   {
      PAMR_get_g_coarsest(&coarsest);
      if (coarsest) return 1;
      valid=PAMR_next_g();
   }

   return 0;
}

//=============================================================================
// sets comm=1 for coarsest grids at level L
//=============================================================================
void set_comm_coarsest(int L)
{
   int valid,coarsest;

   valid=PAMR_init_s_iter(L,PAMR_MGH,0);

   while(valid)
   {
      PAMR_get_g_coarsest(&coarsest);
      PAMR_set_g_comm(coarsest);
      valid=PAMR_next_g();
   }
}

//=============================================================================
// sets comm=1 for all grids at level L
//=============================================================================
void reset_comm(int L)
{
   int valid;

   valid=PAMR_init_s_iter(L,PAMR_MGH,0);

   while(valid)
   {
      PAMR_set_g_comm(1);
      valid=PAMR_next_g();
   }
}

//=============================================================================
// temporarily turns sync'ing on for residual
// (remember to PAMR_thaw... afterwards)
//=============================================================================
void set_res_sync(void)
{
   int i,num;

   num=AMRD_num_elliptic_vars; if (AMRD_is_t0) num+=AMRD_num_elliptic_vars_t0;

   PAMR_freeze_tf_bits();
   PAMR_clear_tf_bits();

   for (i=0; i<num; i++) PAMR_set_tf_bit(AMRD_MG_res_gfn[i],PAMR_SYNC);
}

//=============================================================================
// synchronizes/injects/interpolates all user MG constrained variables
// and elliptic_t0 vars
//=============================================================================
void sync_cnst_data(int L)
{
   int i;

   PAMR_freeze_tf_bits();
   PAMR_clear_tf_bits();

   for (i=0; i<AMRD_num_MG_cnst_data_vars; i++) PAMR_set_tf_bit(AMRD_MG_cnst_data_gfn[i],PAMR_SYNC);
   PAMR_sync(L,0,PAMR_MGH,0);

   PAMR_thaw_tf_bits();
}

void inject_cnst_data(int Lf)
{
   int i;

   //--------------------------------------------------------
   // Before injecting, volume weight the hydro variables.
   //--------------------------------------------------------

   if(AMRD_num_inject_wavg_vars){ 
	apply_wavg(0, Lf-1, PAMR_MGH);
   	apply_wavg(0, Lf, PAMR_MGH);
   }
   PAMR_freeze_tf_bits();
   PAMR_clear_tf_bits();

   for (i=0; i<AMRD_num_MG_cnst_data_vars; i++) {
     if (PAMR_gfn_var_type(AMRD_MG_cnst_data_gfn[i])==PAMR_VERTEX_CENTERED) {
       PAMR_set_tf_bit(AMRD_MG_cnst_data_gfn[i],PAMR_STRAIGHT_INJECT);
     } else {
   	if(AMRD_num_inject_wavg_vars) PAMR_set_tf_bit(AMRD_MG_cnst_data_gfn[i],PAMR_NN_ADD);
   	else PAMR_set_tf_bit(AMRD_MG_cnst_data_gfn[i],PAMR_FIRST_ORDER_CONS);
     }
   }
   PAMR_inject(Lf,0,PAMR_MGH);
   PAMR_thaw_tf_bits();

   if(AMRD_num_inject_wavg_vars){ 
	apply_wavg(1, Lf-1, PAMR_MGH);
   	apply_wavg(1, Lf, PAMR_MGH);
   }	
}

void interp_cnst_data(int Lc)
{
   int i;

   PAMR_freeze_tf_bits();
   PAMR_clear_tf_bits();

   for (i=0; i<AMRD_num_MG_cnst_data_vars; i++) {
     if (PAMR_gfn_var_type(AMRD_MG_cnst_data_gfn[i])==PAMR_VERTEX_CENTERED) {
       PAMR_set_tf_bit(AMRD_MG_cnst_data_gfn[i],PAMR_FOURTH_ORDER);
     } else {
       PAMR_set_tf_bit(AMRD_MG_cnst_data_gfn[i],PAMR_FIRST_ORDER_CONS);
     }
   }
   PAMR_interp(Lc,0,PAMR_MGH); 

   PAMR_thaw_tf_bits();
}

void sync_elliptic_vars_t0_data(int L)
{
   int i;

   PAMR_freeze_tf_bits();
   PAMR_clear_tf_bits();

   for (i=0; i<AMRD_num_elliptic_vars_t0; i++) PAMR_set_tf_bit(AMRD_elliptic_vars_t0_gfn[i],PAMR_SYNC);
   PAMR_sync(L,0,PAMR_MGH,0);

   PAMR_thaw_tf_bits();
}

void inject_elliptic_vars_t0_data(int Lf)
{
   int i;

   PAMR_freeze_tf_bits();
   PAMR_clear_tf_bits();

   for (i=0; i<AMRD_num_elliptic_vars_t0; i++) PAMR_set_tf_bit(AMRD_elliptic_vars_t0_gfn[i],PAMR_STRAIGHT_INJECT);
   PAMR_inject(Lf,0,PAMR_MGH);

   PAMR_thaw_tf_bits();
}

void interp_elliptic_vars_t0_data(int Lc)
{
   int i;

   PAMR_freeze_tf_bits();
   PAMR_clear_tf_bits();

   for (i=0; i<AMRD_num_elliptic_vars_t0; i++) PAMR_set_tf_bit(AMRD_elliptic_vars_t0_gfn[i],PAMR_FOURTH_ORDER);
   PAMR_interp(Lc,0,PAMR_MGH);

   PAMR_thaw_tf_bits();
}

//=============================================================================
// debug utility
//
// which = 0 : all
//       = 1 : MG only
//       = 2 : MG vars only
//=============================================================================
void mg_dump(int L, char *tag, int giter, int which)
{
   real t;
   int i,num;
   char buf[512];
   static real lcount=-1,pgiter=-1;
   real t_lf;

   num=AMRD_num_elliptic_vars; if (AMRD_is_t0) num+=AMRD_num_elliptic_vars_t0;

   if (!AMRD_MG_DV_trace) return;
   t_lf=PAMR_get_time(PAMR_get_max_lev(PAMR_AMRH));
   if (AMRD_MG_DV_trace_t_on>t_lf || AMRD_MG_DV_trace_t_off<t_lf) return;

   if (giter==0 && pgiter!=giter) lcount++;
   pgiter=giter;

   t=lcount+(giter/(AMRD_MG_max_iter+2.0));

   for (i=0; i<num; i++) 
   {
      if (i<AMRD_num_elliptic_vars) 
      {
         PAMR_save_gfn(AMRD_elliptic_vars[i],PAMR_MGH,0,L,t,tag,"");
      }
      else 
      {
         PAMR_save_gfn(AMRD_elliptic_vars_t0[i-AMRD_num_elliptic_vars],PAMR_MGH,0,L,t,tag,"");
      }
   }

   if (which<2)
   {
      for (i=0; i<num; i++) 
      {
         if (i<AMRD_num_elliptic_vars)
         {
	   strcpy(buf,AMRD_elliptic_vars[i]);  
	   AMRD_append_tag(buf,"_lop",AMRD_v_tag,AMRD_c_tag);
	   PAMR_save_gfn(buf,PAMR_MGH,0,L,t,tag,"");

 	   strcpy(buf,AMRD_elliptic_vars[i]);  
	   AMRD_append_tag(buf,"_rhs",AMRD_v_tag,AMRD_c_tag);
	   PAMR_save_gfn(buf,PAMR_MGH,0,L,t,tag,"");

 	   strcpy(buf,AMRD_elliptic_vars[i]);  
	   AMRD_append_tag(buf,"_res",AMRD_v_tag,AMRD_c_tag);
	   PAMR_save_gfn(buf,PAMR_MGH,0,L,t,tag,"");
         }
         else
         {
 	   strcpy(buf,AMRD_elliptic_vars_t0[i-AMRD_num_elliptic_vars]);  
	   AMRD_append_tag(buf,"_lop",AMRD_v_tag,AMRD_c_tag);
	   PAMR_save_gfn(buf,PAMR_MGH,0,L,t,tag,"");
	   
 	   strcpy(buf,AMRD_elliptic_vars_t0[i-AMRD_num_elliptic_vars]);  
	   AMRD_append_tag(buf,"_rhs",AMRD_v_tag,AMRD_c_tag);
	   PAMR_save_gfn(buf,PAMR_MGH,0,L,t,tag,"");

 	   strcpy(buf,AMRD_elliptic_vars_t0[i-AMRD_num_elliptic_vars]);  
	   AMRD_append_tag(buf,"_res",AMRD_v_tag,AMRD_c_tag);
	   PAMR_save_gfn(buf,PAMR_MGH,0,L,t,tag,"");
         }
      }
      if (which<1)
      {
         i=0; 
         while(i<AMRD_num_hyperbolic_vars) PAMR_save_gfn(AMRD_hyperbolic_vars[i++],PAMR_MGH,0,L,t,tag,"");
         PAMR_save_gfn("cmask",PAMR_MGH,0,L,t,tag,"");
         i=0; 
         while(i<AMRD_num_MGH_work_vars) PAMR_save_gfn(AMRD_MGH_work_vars[i++],PAMR_MGH,0,L,t,tag,"");
         i=0; 
         while(i<AMRD_num_AMRH_work_in_MGH_vars) PAMR_save_gfn(AMRD_AMRH_work_in_MGH_vars[i++],PAMR_MGH,0,L,t,tag,"");
      }
   }
}

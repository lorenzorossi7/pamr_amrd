//=============================================================================
// io.c --- reading/writing parameters, files, etc.
//=============================================================================

#include <stdio.h>
#include "io_w.h"
#include "globals_w.h"
#include <bbhutil.h>
#include <stdlib.h>
#include <string.h>
#include "util_w.h"
#include "cls_w.h"
#include "evolve_w.h"

int in_list(char *v, char **list, int num)
{
   int i=0;

   while(i<num)
   {
      if (!(strcmp(v,list[i]))) return 1;
      i++;
   }

   return 0;
}

//=============================================================================
// The following reads parameters associated with amrd, from pfile. At the
// same time, the context is initialized.
//
// First, read all parameters required for a restart. If a restart is
// specified, we return then, else read variable definitions, etc.
//=============================================================================
#define MAX_INIT_DEPTH 20
void init_context(char *pfile)
{
   int i,j,rho_sp[PAMR_MAX_LEVS],rho_tm[PAMR_MAX_LEVS],k,ii,inum[2*PAMR_MAX_DIM];
   int ghost_width[3]={2,2,2},min_width[3]={3,3,3},min_cwidth[3]={3,3,3};
   int interp_buffer=2,gdm_method,stage,sgfn,offset,fcsgfn; 
   int gdm_grid_by_grid,gdm_align,gdm_no_overlap,mod;
   int in_amrh,in_mgh,num_tl,amr_inject,amr_interp,amr_bdy_interp,amr_sync,mg_inject;
   int mg_interp,mg_sync,mg_noinj_to_amr,regrid_transfer;
   int amr_c_to_v,amr_v_to_c;
   int phys_bdy_type[6]={PAMR_UNKNOWN,PAMR_UNKNOWN,PAMR_UNKNOWN,PAMR_UNKNOWN,PAMR_UNKNOWN,PAMR_UNKNOWN};
   real x,lambda;
   char buffer[512],buffer_c[512];
   char chr_name_v[128],chr_name_c[128],*chr_v,*chr_c;
   char wavg_name_c[128], *wavg_c; 
   char *var,*p;
   real init_bbox[4*MAX_INIT_DEPTH*PAMR_MAX_DIM],AMRD_MG_w0;
   int init_lev[2*MAX_INIT_DEPTH],num,lev_num,init_depth;

   int num_amr_inject_v;
   char **amr_inject_v;
   int num_amr_hw_inject_v;
   char **amr_hw_inject_v;
   int num_amr_interp1_v;
   char **amr_interp1_v;
   int num_amr_interp2_v;
   char **amr_interp2_v;
   int num_amr_interp4_v;
   char **amr_interp4_v;
   int num_amr_sync_v;
   char **amr_sync_v;
   int num_amr_transfer2_v;
   char **amr_transfer2_v;
   int num_amr_transfer1_v;
   char **amr_transfer1_v;
   int num_amr_transfer4_v;
   char **amr_transfer4_v;
   int num_amr_c_to_v2_v;
   char **amr_c_to_v2_v;
   int num_amr_v_to_c2_v;
   char **amr_v_to_c2_v;
   int num_amr_bdy_interp1_v;
   char **amr_bdy_interp1_v;
   int num_amr_bdy_interp2_v;
   char **amr_bdy_interp2_v;

   int AMRD_num_TRE_vars;
   char **AMRD_TRE_vars;

   int num_mg_interp2_v;
   char **mg_interp2_v;
   int num_mg_hw_restr_v;
   char **mg_hw_restr_v;
   int num_mg_fw_restr_v;
   char **mg_fw_restr_v;
   int num_mg_hw_restr_bdy_si_v; // hwr in interior, straight injection on boundary
   char **mg_hw_restr_bdy_si_v;
   int num_mg_sync_v;
   char **mg_sync_v;

   int num_even_vars_x0min,num_even_vars_x1min,num_even_vars_x2min;
   int num_even_vars_x0max,num_even_vars_x1max,num_even_vars_x2max;
   char **even_vars_x0min,**even_vars_x1min,**even_vars_x2min;
   char **even_vars_x0max,**even_vars_x1max,**even_vars_x2max;

   int num_odd_vars_x0min,num_odd_vars_x1min,num_odd_vars_x2min;
   int num_odd_vars_x0max,num_odd_vars_x1max,num_odd_vars_x2max;
   char **odd_vars_x0min,**odd_vars_x1min,**odd_vars_x2min;
   char **odd_vars_x0max,**odd_vars_x1max,**odd_vars_x2max;

   //------------------------------------------------------------------------
   // Set defaults (some initialized in var. declarations above)
   //------------------------------------------------------------------------
   AMRD_using_cc_tre=-1; // temp value ... to check whether TRE vars aren't being mixed
   AMRD_echo_params=1;
   AMRD_ID_DV_trace=0;
   AMRD_regrid_trace=0;
   AMRD_MG_trace=0;
   AMRD_MG_DV_trace=-1;
   AMRD_MG_DV_trace_t_on=1e-15;
   AMRD_MG_DV_trace_t_off=1e-15;
   AMRD_evo_trace=0;
   AMRD_evo_DV_trace=0;
   AMRD_evo_DV_trace_t_on=1e-15;
   AMRD_evo_DV_trace_t_off=1e-15;
   AMRD_base_shape[0]=AMRD_base_shape[1]=AMRD_base_shape[2]=3;
   AMRD_base_bbox[0]=AMRD_base_bbox[2]=AMRD_base_bbox[4]=-1;
   AMRD_base_bbox[1]=AMRD_base_bbox[3]=AMRD_base_bbox[5]=1;
   AMRD_max_lev=1;
   for (i=0; i<PAMR_MAX_LEVS; i++) rho_sp[i]=rho_tm[i]=2;
   AMRD_t0=0;
   AMRD_steps=1;
   AMRD_evo_max_iter=50;
   AMRD_evo_min_iter=1;
   AMRD_MG_start_iter=1;
   AMRD_MG_max_iter=50;
   AMRD_MG_min_iter=1;
   AMRD_MG_max_citer=5000;
   AMRD_MG_pre_swp=3;
   AMRD_MG_pst_swp=3;
   AMRD_np1_initial_guess=0;
   AMRD_evo_tol=0;
   AMRD_MG_tol=0;
   AMRD_MG_crtol=1.0e-3;
   AMRD_MG_w0=1.0;
   AMRD_MG_extrap_method=AMRD_MG_EXTRAP_2ND_FROM_ETM1;
   AMRD_MG_eps_c=1.0;
   AMRD_MG_reinterp_bdy=0;
   lambda=1.0;
   AMRD_id_method=0;
   AMRD_id_pl_method=0;
   AMRD_id_user_mg_pl=0;
   AMRD_id_pl_steps=1;
   AMRD_id_pl_lambda=0.1;
   AMRD_TRE_max=0;
   AMRD_TRE_norm=0;
   AMRD_TRE_buffer=0;
   AMRD_TRE_ibc_buffer=0;
   AMRD_TRE_ibc_a_buffer=0;
   AMRD_TRE_exc_buffer=0;
   AMRD_TRE_exc_buffer_lmin=0;
   AMRD_TRE_sgpbh=0;
   AMRD_TRE_ibcp_buffer=0;
   AMRD_cls_method=CLS_SIMPLE;
   AMRD_cls_merge_dist=0;
   AMRD_cls_align_mode=CLS_ALIGN_SHRINK;
   gdm_method=0;
   gdm_grid_by_grid=0;
   gdm_align=0;
   gdm_no_overlap=0;
   AMRD_save_tag=0;
   AMRD_save_all_finer=0;
   AMRD_save_iter_mod=0;
   AMRD_c_tag=0;
   AMRD_v_tag=0;
   // AMRD_AMR_bdy_width=2 makes more sense for consistency with AMRD_AMR_bdy_width_c,
   // however for consistency with legacy codes we will place the onus on
   // new codes to manually set this to 2 (or whatever).
   AMRD_AMR_bdy_width = 1; 
   AMRD_AMR_bdy_width_c = 1;
   AMRD_resid_in_evo = 1;
   AMRD_global_var_norm_floor=1;
   AMRD_calc_global_var_norms=0;
   AMRD_global_var_norm_type=AMRD_GLOBAL_VAR_INF_NORM;
   for (i=0; i<PAMR_MAX_GFNS; i++) AMRD_global_var_norms[i]=1;
   if (!(AMRD_save_ivec0=(int *)malloc(sizeof(int)*(AMRD_MAX_IVEC_SIZE+1))))
      AMRD_stop("init_context ... out of memory\n","");
   sget_ivec_param("ivec := 0","ivec", AMRD_save_ivec0, AMRD_MAX_IVEC_SIZE);
   for (i=0; i<PAMR_MAX_LEVS; i++) 
   {
      if (!(AMRD_save_ivec[i]=(int *)malloc(sizeof(int)*(AMRD_MAX_IVEC_SIZE+1))))
         AMRD_stop("init_context ... out of memory\n","");
      sget_ivec_param("ivec := 0","ivec", AMRD_save_ivec[i], AMRD_MAX_IVEC_SIZE);
      AMRD_lsteps[i]=0;
      AMRD_tre_valid[i]=0;
      //      AMRD_TRE_max_hydro[i]=1.e+100;
   }
   AMRD_skip_frg=1;
   AMRD_regrid_interval=4;
   AMRD_regrid_min_lev=1;
   AMRD_periodic[0]=0;
   AMRD_periodic[1]=0;
   AMRD_periodic[2]=0;

   AMRD_regrid_script_name=0;
   AMRD_regrid_script=0;

   AMRD_num_f_tre_vars=0;
   for (i=0; i<AMRD_MAX_VARS; i++) AMRD_f_tre_vars[i]=0;

   AMRD_eps_diss=0;
   AMRD_diss_bdy=0;
   AMRD_diss_stride=1;
   AMRD_diss_use_6th_order=0;
   AMRD_repop_diss_bdy=0;
   AMRD_diss_freq=1;
   AMRD_rg_eps_diss=0;
   AMRD_mg_tre_gfn=0;
   AMRD_evo_ssc=1;
   AMRD_ic_n=2;
   AMRD_interp_AMR_bdy_order=4;
   AMRD_max_t_interp_order=2;
   AMRD_re_interp_width=0;
   AMRD_re_interp_width_c=0;
   AMRD_c_to_v_in_tstep=0;
   AMRD_v_to_c_in_tstep=0;

   AMRD_ex=1;
   AMRD_do_ex=0;

   AMRD_cp_restart=0;
   AMRD_cp_delta_t_hrs=0;
   AMRD_cp_first_delta_t_hrs=0;
   AMRD_cp_save_fname=0;
   AMRD_cp_restore_fname=0;

   if (my_rank==0) printf("Reading parameters from file %s\n\n",pfile);

   AMRD_int_param(pfile,"echo_params",&AMRD_echo_params,1);

   i=0;
   AMRD_int_param(pfile,"pamr_trace_lev",&i,1);
   PAMR_set_trace_lev(i);

   AMRD_int_param(pfile,"dim",&AMRD_dim,1);

   AMRD_int_param(pfile,"max_lev",&AMRD_max_lev,1);

   AMRD_int_param(pfile,"save_iter_mod",&AMRD_save_iter_mod,1);
   AMRD_int_param(pfile,"save_all_finer",&AMRD_save_all_finer,1);

   if (AMRD_save_iter_mod>1)
   printf("\nWARNING\n amrd: currently AMRD_save_iter_mod > 1 will not work correctly with non-trivial ivec's\nWARNING\n");

   if (AMRD_dim<1 || AMRD_dim>3) AMRD_stop("amrd: AMRD_dim<1 or AMRD_dim>3 not supported\n",0);

   AMRD_int_param(pfile,"MG_trace",&AMRD_MG_trace,1);
   AMRD_int_param(pfile,"MG_DV_trace",&AMRD_MG_DV_trace,1);
   AMRD_int_param(pfile,"regrid_trace",&AMRD_regrid_trace,1);
   AMRD_real_param(pfile,"MG_DV_trace_t_on",&AMRD_MG_DV_trace_t_on,1);
   AMRD_real_param(pfile,"MG_DV_trace_t_off",&AMRD_MG_DV_trace_t_off,1);
   AMRD_int_param(pfile,"ID_DV_trace",&AMRD_ID_DV_trace,1);
   AMRD_int_param(pfile,"evo_DV_trace",&AMRD_evo_DV_trace,1);
   AMRD_real_param(pfile,"evo_DV_trace_t_on",&AMRD_evo_DV_trace_t_on,1);
   AMRD_real_param(pfile,"evo_DV_trace_t_off",&AMRD_evo_DV_trace_t_off,1);
   AMRD_int_param(pfile,"evo_trace",&AMRD_evo_trace,1);

   AMRD_int_param(pfile,"steps",&AMRD_steps,1);
   AMRD_int_param(pfile,"evo_max_iter",&AMRD_evo_max_iter,1);
   AMRD_int_param(pfile,"evo_min_iter",&AMRD_evo_min_iter,1);
   AMRD_int_param(pfile,"MG_start_iter",&AMRD_MG_start_iter,1);
   AMRD_int_param(pfile,"MG_max_iter",&AMRD_MG_max_iter,1);
   AMRD_int_param(pfile,"MG_min_iter",&AMRD_MG_min_iter,1);
   AMRD_int_param(pfile,"MG_max_citer",&AMRD_MG_max_citer,1);
   AMRD_int_param(pfile,"MG_pre_swp",&AMRD_MG_pre_swp,1);
   AMRD_int_param(pfile,"MG_pst_swp",&AMRD_MG_pst_swp,1);

   AMRD_int_param(pfile,"do_ex",&AMRD_do_ex,1);
   if (AMRD_do_ex<0 && my_rank==0) printf("NOTE: AMRD_do_ex < 0 ... using chr as an eps array!\n");
   AMRD_real_param(pfile,"ex",&AMRD_ex,1);

   AMRD_real_param(pfile,"evo_tol",&AMRD_evo_tol,1);
   AMRD_real_param(pfile,"MG_tol",&AMRD_MG_tol,1);
   AMRD_real_param(pfile,"MG_crtol",&AMRD_MG_crtol,1);
   AMRD_real_param(pfile,"MG_w0",&AMRD_MG_w0,1);
   AMRD_real_param(pfile,"MG_w0",&AMRD_MG_w0,1);
   AMRD_MG_w0_r=AMRD_MG_w0_i=AMRD_MG_w0;
   AMRD_real_param(pfile,"MG_w0_r",&AMRD_MG_w0_r,1);
   AMRD_real_param(pfile,"MG_w0_i",&AMRD_MG_w0_i,1);
   AMRD_int_param(pfile,"evo_ssc",&AMRD_evo_ssc,1);
   AMRD_int_param(pfile,"MG_extrap_method",&AMRD_MG_extrap_method,1);
   AMRD_real_param(pfile,"MG_eps_c",&AMRD_MG_eps_c,1);
   AMRD_int_param(pfile,"MG_reinterp_bdy",&AMRD_MG_reinterp_bdy,1);
   AMRD_real_param(pfile,"TRE_max",&AMRD_TRE_max,1);
   AMRD_int_param(pfile,"TRE_norm",&AMRD_TRE_norm,1);
   AMRD_int_param(pfile,"TRE_buffer",&AMRD_TRE_buffer,1);
   AMRD_int_param(pfile,"TRE_ibc_buffer",&AMRD_TRE_ibc_buffer,1);
   AMRD_int_param(pfile,"TRE_ibc_a_buffer",&AMRD_TRE_ibc_a_buffer,1);
   AMRD_int_param(pfile,"TRE_exc_buffer",&AMRD_TRE_exc_buffer,1);
   AMRD_int_param(pfile,"TRE_exc_buffer_lmin",&AMRD_TRE_exc_buffer_lmin,1);
   AMRD_int_param(pfile,"TRE_sgpbh",&AMRD_TRE_sgpbh,1);
   AMRD_int_param(pfile,"TRE_ibcp_buffer",&AMRD_TRE_ibcp_buffer,1);
   AMRD_int_param(pfile,"regrid_interval",&AMRD_regrid_interval,1);
   AMRD_int_param(pfile,"regrid_min_lev",&AMRD_regrid_min_lev,1);
   if (AMRD_regrid_min_lev<1) AMRD_stop("Error ... regrid_min_lev>=1","");
   AMRD_int_param(pfile,"cls_method",&AMRD_cls_method,1);
   AMRD_int_param(pfile,"cls_merge_dist",&AMRD_cls_merge_dist,1);
   AMRD_int_param(pfile,"cls_align_mode",&AMRD_cls_align_mode,1);
   AMRD_int_param(pfile,"calc_global_var_norms",&AMRD_calc_global_var_norms,1);
   AMRD_int_param(pfile,"global_var_norm_type",&AMRD_global_var_norm_type,1);
   AMRD_real_param(pfile,"global_var_norm_floor",&AMRD_global_var_norm_floor,1);

   AMRD_int_param(pfile,"magic_cookie",&AMRD_magic_cookie,1);

   AMRD_int_param(pfile,"id_method",&AMRD_id_method,1);
   AMRD_int_param(pfile,"id_user_mg_pl",&AMRD_id_user_mg_pl,1);
   AMRD_int_param(pfile,"id_pl_method",&AMRD_id_pl_method,1);
   AMRD_int_param(pfile,"id_pl_steps",&AMRD_id_pl_steps,1);
   if (AMRD_id_pl_steps<1) AMRD_stop("AMRD_id_pl_steps must be >=1\n","");
   AMRD_real_param(pfile,"id_pl_lambda",&AMRD_id_pl_lambda,1);

   AMRD_t_interp_substeps=0; AMRD_num_t_interp_substeps=0;
   AMRD_real_param_v(pfile,"t_interp_substeps",&AMRD_t_interp_substeps,&AMRD_num_t_interp_substeps);

   AMRD_num_hyperbolic_vars=0; AMRD_hyperbolic_vars=0;
   AMRD_str_param_v(pfile,"hyperbolic_vars",&AMRD_hyperbolic_vars,&AMRD_num_hyperbolic_vars);
   if (AMRD_num_hyperbolic_vars>=AMRD_MAX_VARS) 
      AMRD_stop("number of hyperbolic_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_num_elliptic_vars=0; AMRD_elliptic_vars=0;
   AMRD_str_param_v(pfile,"elliptic_vars",&AMRD_elliptic_vars,&AMRD_num_elliptic_vars);
   if (AMRD_num_elliptic_vars>=AMRD_MAX_VARS) 
      AMRD_stop("number of elliptic_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

/*    if (AMRD_num_elliptic_vars && AMRD_num_t_interp_substeps) */
/*       AMRD_stop("AMRD FATAL WARNING: NOT YET TESTED WITH COMBINATION OF AMRD_num_elliptic_vars and " */
/*                 "t_interp_substeps.\n" */
/*                 "Remove this line from the code and proceed with caution\n",""); */

   AMRD_num_elliptic_vars_t0=0; AMRD_elliptic_vars_t0=0;
   AMRD_str_param_v(pfile,"elliptic_vars_t0",&AMRD_elliptic_vars_t0,&AMRD_num_elliptic_vars_t0);
   if (AMRD_num_elliptic_vars_t0>=AMRD_MAX_VARS) 
      AMRD_stop("number of elliptic_vars_t0 is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_num_AMRH_work_vars=0; AMRD_AMRH_work_vars=0;
   AMRD_str_param_v(pfile,"AMRH_work_vars",&AMRD_AMRH_work_vars,&AMRD_num_AMRH_work_vars);
   if (AMRD_num_AMRH_work_vars>=AMRD_MAX_VARS) 
      AMRD_stop("number of AMRH_work_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_num_AMRH_work_in_MGH_vars=0; AMRD_AMRH_work_in_MGH_vars=0;
   AMRD_str_param_v(pfile,"AMRH_work_in_MGH_vars",&AMRD_AMRH_work_in_MGH_vars,&AMRD_num_AMRH_work_in_MGH_vars);
   if (AMRD_num_AMRH_work_in_MGH_vars>=AMRD_MAX_VARS) 
      AMRD_stop("number of AMRH_work_in_MGH_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_num_fc_vars=0; AMRD_fc_vars=0;
   AMRD_str_param_v(pfile,"fc_vars",&AMRD_fc_vars,&AMRD_num_fc_vars);
   if (AMRD_num_fc_vars>=AMRD_MAX_VARS) 
      AMRD_stop("number of fc_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_num_past_bdy_interp_vars=0; AMRD_past_bdy_interp_vars=0;
   AMRD_str_param_v(pfile,"past_bdy_interp_vars",&AMRD_past_bdy_interp_vars,&AMRD_num_past_bdy_interp_vars);
   if (AMRD_num_past_bdy_interp_vars>=AMRD_MAX_VARS) 
      AMRD_stop("number of past_bdy_interp_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_num_work_repop_vars=0; AMRD_work_repop_vars=0;
   AMRD_str_param_v(pfile,"work_repop_vars",&AMRD_work_repop_vars,&AMRD_num_work_repop_vars);
   if (AMRD_num_work_repop_vars>=AMRD_MAX_VARS) 
      AMRD_stop("number of work_repop_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_num_AMRH_work_inject_vars=0; AMRD_AMRH_work_inject_vars=0;
   AMRD_str_param_v(pfile,"AMRH_work_inject_vars",&AMRD_AMRH_work_inject_vars,&AMRD_num_AMRH_work_inject_vars);
   if (AMRD_num_AMRH_work_inject_vars>=AMRD_MAX_VARS) 
      AMRD_stop("number of AMRH_work_inject_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_num_MGH_work_vars=0; AMRD_MGH_work_vars=0;
   AMRD_str_param_v(pfile,"MGH_work_vars",&AMRD_MGH_work_vars,&AMRD_num_MGH_work_vars);
   if (AMRD_num_MGH_work_vars>=AMRD_MAX_VARS) 
      AMRD_stop("number of MGH_work_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_num_MG_cnst_data_vars=0; AMRD_MG_cnst_data_vars=0;
   AMRD_str_param_v(pfile,"MG_cnst_data_vars",&AMRD_MG_cnst_data_vars,&AMRD_num_MG_cnst_data_vars); 
   if (AMRD_num_MG_cnst_data_vars>=AMRD_MAX_VARS) 
      AMRD_stop("num_MG_cnst_data_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_num_free_data_vars=0; AMRD_free_data_vars=0;
   AMRD_str_param_v(pfile,"free_data_vars",&AMRD_free_data_vars,&AMRD_num_free_data_vars); 
   if (AMRD_num_free_data_vars>=AMRD_MAX_VARS) 
      AMRD_stop("num_free_data_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_num_ex_repop1_vars=0; AMRD_ex_repop1_vars=0;
   AMRD_str_param_v(pfile,"ex_repop1_vars",&AMRD_ex_repop1_vars,&AMRD_num_ex_repop1_vars); 
   if (AMRD_num_ex_repop1_vars>=AMRD_MAX_VARS) 
      AMRD_stop("num_ex_repop1_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_num_ex_repop2_vars=0; AMRD_ex_repop2_vars=0;
   AMRD_str_param_v(pfile,"ex_repop2_vars",&AMRD_ex_repop2_vars,&AMRD_num_ex_repop2_vars); 
   if (AMRD_num_ex_repop2_vars>=AMRD_MAX_VARS) 
      AMRD_stop("num_ex_repop2_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_num_ex_repop3_vars=0; AMRD_ex_repop3_vars=0;
   AMRD_str_param_v(pfile,"ex_repop3_vars",&AMRD_ex_repop3_vars,&AMRD_num_ex_repop3_vars); 
   if (AMRD_num_ex_repop3_vars>=AMRD_MAX_VARS) 
      AMRD_stop("num_ex_repop3_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_num_ex_repop4_vars=0; AMRD_ex_repop4_vars=0;
   AMRD_str_param_v(pfile,"ex_repop4_vars",&AMRD_ex_repop4_vars,&AMRD_num_ex_repop4_vars); 
   if (AMRD_num_ex_repop4_vars>=AMRD_MAX_VARS) 
      AMRD_stop("num_ex_repop4_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_num_tnp1_diss_vars=0; AMRD_tnp1_diss_vars=0;
   AMRD_str_param_v(pfile,"tnp1_diss_vars",&AMRD_tnp1_diss_vars,&AMRD_num_tnp1_diss_vars); 
   if (AMRD_num_tnp1_diss_vars>=AMRD_MAX_VARS) 
      AMRD_stop("num_tnp1_diss_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_diss_all_past=0;
   AMRD_int_param(pfile,"diss_freq",&AMRD_diss_freq,1);
   if (AMRD_diss_freq>1) AMRD_diss_all_past=1;
   AMRD_int_param(pfile,"diss_all_past",&AMRD_diss_all_past,1);
   
   AMRD_num_tn_diss_vars=0; AMRD_tn_diss_vars=0;
   AMRD_str_param_v(pfile,"tn_diss_vars",&AMRD_tn_diss_vars,&AMRD_num_tn_diss_vars); 
   if (AMRD_num_tn_diss_vars>=AMRD_MAX_VARS) 
      AMRD_stop("num_tn_diss_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 
   
   int num_unused_vars=0; char** unused_vars=0;
   AMRD_str_param_v(pfile,"unused_vars",&unused_vars,&num_unused_vars); 
   if (num_unused_vars>=AMRD_MAX_VARS) 
      AMRD_stop("num_unused_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 
   if (num_unused_vars>0) printf ("WARNING!:  unused_vars options currently does not work with check pointing.\n"); 

   AMRD_real_param(pfile,"eps_diss",&AMRD_eps_diss,1);
   AMRD_tn_eps_diss=AMRD_eps_diss;
   AMRD_tnp1_eps_diss=AMRD_eps_diss;

   AMRD_real_param(pfile,"tn_eps_diss",&AMRD_tn_eps_diss,1);
   AMRD_real_param(pfile,"tnp1_eps_diss",&AMRD_tnp1_eps_diss,1);

   AMRD_int_param(pfile,"diss_bdy",&AMRD_diss_bdy,1);
   AMRD_repop_diss_bdy=AMRD_diss_bdy;
   AMRD_int_param(pfile,"repop_diss_bdy",&AMRD_repop_diss_bdy,1);

   AMRD_int_param(pfile,"diss_use_6th_order",&AMRD_diss_use_6th_order,1);
   AMRD_int_param(pfile,"diss_stride",&AMRD_diss_stride,1);

   AMRD_num_tnp1_liipb_vars=0; AMRD_tnp1_liipb_vars=0;
   AMRD_str_param_v(pfile,"tnp1_liipb_vars",&AMRD_tnp1_liipb_vars,&AMRD_num_tnp1_liipb_vars); 
   if (AMRD_num_tnp1_liipb_vars>=AMRD_MAX_VARS) 
      AMRD_stop("num_tnp1_liipb_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_num_tnp1_liiab_vars=0; AMRD_tnp1_liiab_vars=0;
   AMRD_str_param_v(pfile,"tnp1_liiab_vars",&AMRD_tnp1_liiab_vars,&AMRD_num_tnp1_liiab_vars); 
   if (AMRD_num_tnp1_liiab_vars>=AMRD_MAX_VARS) 
      AMRD_stop("num_tnp1_liiab_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_num_tnp1_liibb_vars=0; AMRD_tnp1_liibb_vars=0;
   AMRD_str_param_v(pfile,"tnp1_liibb_vars",&AMRD_tnp1_liibb_vars,&AMRD_num_tnp1_liibb_vars); 
   if (AMRD_num_tnp1_liibb_vars>=AMRD_MAX_VARS) 
      AMRD_stop("num_tnp1_liibb_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_num_rg_diss_vars=0; AMRD_rg_diss_vars=0;
   AMRD_str_param_v(pfile,"rg_diss_vars",&AMRD_rg_diss_vars,&AMRD_num_rg_diss_vars); 
   if (AMRD_num_rg_diss_vars>=AMRD_MAX_VARS) 
      AMRD_stop("num_rg_diss_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_real_param(pfile,"rg_eps_diss",&AMRD_rg_eps_diss,1);
   
   AMRD_num_interp_AMR_bdy_vars=0; AMRD_interp_AMR_bdy_vars=0;
   AMRD_str_param_v(pfile,"interp_AMR_bdy_vars",&AMRD_interp_AMR_bdy_vars,&AMRD_num_interp_AMR_bdy_vars); 
   if (AMRD_num_interp_AMR_bdy_vars>=AMRD_MAX_VARS) 
      AMRD_stop("num_interp_AMR_bdy_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_int_param(pfile,"interp_AMR_bdy_order",&AMRD_interp_AMR_bdy_order,1);
   AMRD_int_param(pfile,"max_t_interp_order",&AMRD_max_t_interp_order,1);

   AMRD_int_param(pfile,"re_interp_width",&AMRD_re_interp_width,1);
   AMRD_int_param(pfile,"re_interp_width_c",&AMRD_re_interp_width_c,1);

   AMRD_int_param(pfile,"v_to_c_in_tstep",&AMRD_v_to_c_in_tstep,1);
   AMRD_int_param(pfile,"c_to_v_in_tstep",&AMRD_c_to_v_in_tstep,1);

   AMRD_num_TRE_vars=0; AMRD_TRE_vars=0;
   AMRD_str_param_v(pfile,"TRE_vars",&AMRD_TRE_vars,&AMRD_num_TRE_vars); 
   if (AMRD_num_TRE_vars>AMRD_MAX_VARS) AMRD_stop("num_TRE_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_int_param(pfile,"regrid_script",&AMRD_regrid_script,1);
   if (AMRD_regrid_script) AMRD_str_param(pfile,"regrid_script_name",&AMRD_regrid_script_name,1);

   AMRD_str_param(pfile,"c_tag",&AMRD_c_tag,1); 
   AMRD_str_param(pfile,"v_tag",&AMRD_v_tag,1); 
   AMRD_int_param(pfile,"AMR_bdy_width",  &AMRD_AMR_bdy_width,1);
   AMRD_int_param(pfile,"AMR_bdy_width_c",&AMRD_AMR_bdy_width_c,1);

   AMRD_int_param(pfile,"resid_in_evo",&AMRD_resid_in_evo,1);

   AMRD_str_param(pfile,"save_tag",&AMRD_save_tag,1); 
   AMRD_ivec_param(pfile,"save_ivec",AMRD_save_ivec0,AMRD_MAX_IVEC_SIZE);
   for (i=0; i<AMRD_max_lev; i++)
   {
      sprintf(buffer,"save_ivec_%i",i+1);
      AMRD_ivec_param(pfile,buffer,AMRD_save_ivec[i],AMRD_MAX_IVEC_SIZE);
   }

   AMRD_num_save_mg_vars=0; AMRD_save_mg_vars=0;
   AMRD_str_param_v(pfile,"save_mg_vars",&AMRD_save_mg_vars,&AMRD_num_save_mg_vars);
   if (AMRD_num_save_mg_vars>=AMRD_MAX_VARS) AMRD_stop("num_save_mg_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   for (j=0; j<AMRD_MAX_TIMES; j++)
   {
      AMRD_num_save_n_vars[j]=0; AMRD_save_n_vars[j]=0;
      sprintf(buffer,"save_%i_vars",j+1);
      AMRD_str_param_v(pfile,buffer,&(AMRD_save_n_vars[j]),&(AMRD_num_save_n_vars[j])); 
      if (AMRD_num_save_n_vars[j]>=AMRD_MAX_VARS) AMRD_stop(buffer," is greater than maximum allowed ... change AMRD_MAX_VARS and recompile"); 
   }
   
   num_amr_inject_v=0; amr_inject_v=0;
   AMRD_str_param_v(pfile,"amr_inject",&amr_inject_v,&num_amr_inject_v); 
   if (num_amr_inject_v>=AMRD_MAX_VARS) AMRD_stop("num_amr_inject is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 
   
   num_amr_hw_inject_v=0; amr_hw_inject_v=0;
   AMRD_str_param_v(pfile,"amr_hw_inject",&amr_hw_inject_v,&num_amr_hw_inject_v); 
   if (num_amr_hw_inject_v>=AMRD_MAX_VARS) AMRD_stop("num_amr_hw_inject is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 
   
   AMRD_num_inject_wavg_vars=0; AMRD_inject_wavg_vars=0;
   AMRD_str_param_v(pfile,"amr_inject_wavg",&AMRD_inject_wavg_vars,&AMRD_num_inject_wavg_vars); 
   if (AMRD_num_inject_wavg_vars>=AMRD_MAX_VARS) AMRD_stop("AMRD_num_inject_wavg_vars is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_amr_interp1_v=0; amr_interp1_v=0;
   AMRD_str_param_v(pfile,"amr_interp1",&amr_interp1_v,&num_amr_interp1_v); 
   if (num_amr_interp1_v>=AMRD_MAX_VARS) AMRD_stop("num_amr_interp1 is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_amr_interp2_v=0; amr_interp2_v=0;
   AMRD_str_param_v(pfile,"amr_interp2",&amr_interp2_v,&num_amr_interp2_v); 
   if (num_amr_interp2_v>=AMRD_MAX_VARS) AMRD_stop("num_amr_interp2 is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_amr_interp4_v=0; amr_interp4_v=0;
   AMRD_str_param_v(pfile,"amr_interp4",&amr_interp4_v,&num_amr_interp4_v); 
   if (num_amr_interp4_v>=AMRD_MAX_VARS) AMRD_stop("num_amr_interp4 is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_amr_sync_v=0; amr_sync_v=0;
   AMRD_str_param_v(pfile,"amr_sync",&amr_sync_v,&num_amr_sync_v); 
   if (num_amr_sync_v>=AMRD_MAX_VARS) AMRD_stop("num_amr_sync is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_amr_transfer1_v=0; amr_transfer1_v=0;
   AMRD_str_param_v(pfile,"amr_transfer1",&amr_transfer1_v,&num_amr_transfer1_v); 
   if (num_amr_transfer1_v>=AMRD_MAX_VARS) AMRD_stop("num_amr_transfer1 is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 
   
   num_amr_transfer2_v=0; amr_transfer2_v=0;
   AMRD_str_param_v(pfile,"amr_transfer2",&amr_transfer2_v,&num_amr_transfer2_v); 
   if (num_amr_transfer2_v>=AMRD_MAX_VARS) AMRD_stop("num_amr_transfer2 is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_amr_transfer4_v=0; amr_transfer4_v=0;
   AMRD_str_param_v(pfile,"amr_transfer4",&amr_transfer4_v,&num_amr_transfer4_v); 
   if (num_amr_transfer4_v>=AMRD_MAX_VARS) AMRD_stop("num_amr_transfer4 is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_amr_c_to_v2_v=0; amr_c_to_v2_v=0;
   AMRD_str_param_v(pfile,"amr_c_to_v2",&amr_c_to_v2_v,&num_amr_c_to_v2_v); 
   if (num_amr_c_to_v2_v>=AMRD_MAX_VARS) AMRD_stop("num_amr_c_to_v2 is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_amr_v_to_c2_v=0; amr_v_to_c2_v=0;
   AMRD_str_param_v(pfile,"amr_v_to_c2",&amr_v_to_c2_v,&num_amr_v_to_c2_v); 
   if (num_amr_v_to_c2_v>=AMRD_MAX_VARS) AMRD_stop("num_amr_v_to_c2 is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_amr_bdy_interp1_v=0; amr_bdy_interp1_v=0;
   AMRD_str_param_v(pfile,"amr_bdy_interp1",&amr_bdy_interp1_v,&num_amr_bdy_interp1_v); 
   if (num_amr_bdy_interp1_v>=AMRD_MAX_VARS) AMRD_stop("num_amr_bdy_interp1 is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_amr_bdy_interp2_v=0; amr_bdy_interp2_v=0;
   AMRD_str_param_v(pfile,"amr_bdy_interp2",&amr_bdy_interp2_v,&num_amr_bdy_interp2_v); 
   if (num_amr_bdy_interp2_v>=AMRD_MAX_VARS) AMRD_stop("num_amr_bdy_interp2 is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_mg_interp2_v=0; mg_interp2_v=0;
   AMRD_str_param_v(pfile,"mg_interp2",&mg_interp2_v,&num_mg_interp2_v); 
   if (num_mg_interp2_v>=AMRD_MAX_VARS) AMRD_stop("num_mg_interp2 is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_mg_hw_restr_v=0; mg_hw_restr_v=0;
   AMRD_str_param_v(pfile,"mg_hw_restr",&mg_hw_restr_v,&num_mg_hw_restr_v); 
   if (num_mg_hw_restr_v>=AMRD_MAX_VARS) AMRD_stop("num_mg_hw_restr is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_mg_fw_restr_v=0; mg_fw_restr_v=0;
   AMRD_str_param_v(pfile,"mg_fw_restr",&mg_fw_restr_v,&num_mg_fw_restr_v); 
   if (num_mg_fw_restr_v>=AMRD_MAX_VARS) AMRD_stop("num_mg_fw_restr is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_mg_hw_restr_bdy_si_v=0; mg_hw_restr_bdy_si_v=0;  // hwr in the interior, straight injection on the boundary
   AMRD_str_param_v(pfile,"mg_hw_restr_bdy_si",&mg_hw_restr_bdy_si_v,&num_mg_hw_restr_bdy_si_v); 
   if (num_mg_hw_restr_bdy_si_v>=AMRD_MAX_VARS) AMRD_stop("num_mg_hw_restr_bdy_si is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_mg_sync_v=0; mg_sync_v=0;
   AMRD_str_param_v(pfile,"mg_sync",&mg_sync_v,&num_mg_sync_v); 
   if (num_mg_sync_v>=AMRD_MAX_VARS) AMRD_stop("num_mg_sync is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   AMRD_int_param(pfile,"base_shape",AMRD_base_shape,AMRD_dim);
   AMRD_real_param(pfile,"base_bbox",AMRD_base_bbox,2*AMRD_dim);

   AMRD_real_param(pfile,"cp_delta_t_hrs",&AMRD_cp_delta_t_hrs,1);
   if (AMRD_cp_delta_t_hrs>0) AMRD_str_param(pfile,"cp_save_fname",&AMRD_cp_save_fname,1); 
   AMRD_int_param(pfile,"cp_restart",&AMRD_cp_restart,1);
   if (AMRD_cp_restart) AMRD_str_param(pfile,"cp_restore_fname",&AMRD_cp_restore_fname,1);
   AMRD_cp_first_delta_t_hrs=AMRD_cp_delta_t_hrs;
   AMRD_real_param(pfile,"cp_first_delta_t_hrs",&AMRD_cp_first_delta_t_hrs,1);
   if (AMRD_cp_first_delta_t_hrs<=0) AMRD_cp_first_delta_t_hrs=AMRD_cp_delta_t_hrs;

   AMRD_int_param(pfile,"np1_initial_guess",&AMRD_np1_initial_guess,1);

   if (!(pamr_context=PAMR_init_context(0,0,AMRD_dim,AMRD_base_shape,AMRD_base_bbox,AMRD_cp_restore_fname,AMRD_v_tag,AMRD_c_tag)))
      AMRD_stop("PAMR_init_context failed\n",0);

   AMRD_int_param(pfile,"num_evo_tl",&AMRD_num_evo_tl,1);

   AMRD_int_param(pfile,"periodic",AMRD_periodic,AMRD_dim);

   // =========================================================================
   // parameters after this point are all i.c. or PAMR ones that should
   // have been saved and restored by pamr. ... except AMRD_f_tre_vars... idiot
   // =========================================================================
   if (AMRD_cp_restart) {
     i=0;j=0;
     while(i<AMRD_num_hyperbolic_vars) {
       if (in_list(AMRD_hyperbolic_vars[i],AMRD_TRE_vars,AMRD_num_TRE_vars)) {
	 strcpy(buffer,AMRD_hyperbolic_vars[i]); 
	 AMRD_append_tag(buffer,"_tre",AMRD_v_tag,AMRD_c_tag);
	 AMRD_f_tre_vars[j++]=strdup(buffer);
       }
       i++;
     }
     if (!(amrd_do_cp(AMRD_cp_restore_fname,AMRD_CP_RESTORE))) AMRD_stop("amrd_do_cp failed\n","");
     return;
   }

   PAMR_set_periodic_bdy(AMRD_periodic);

   AMRD_real_param(pfile,"t0",&AMRD_t0,1);

   AMRD_int_param(pfile,"ic_n",&AMRD_ic_n,1);
   if (AMRD_ic_n!=2) AMRD_stop("amrd: AMRD_ic_n!=2 not yet supported\n","");
   AMRD_int_param(pfile,"skip_frg",&AMRD_skip_frg,1);

   init_depth=2;
   AMRD_int_param(pfile,"init_depth",&init_depth,1);
   if (init_depth < 2 || init_depth > MAX_INIT_DEPTH) AMRD_stop("amrd: invalid init_depth","");

   //--------------------------------------------------------------------------
   // "HISTORIC" NOTES:
   //
   // PAMR currently cannot guarantee that interpolation from a 
   // parent level will give *exact* results on all corner/edge points
   // of the child, especially when the edges of a child aren't aligned
   // with those of it's parent, which *does* occur when grids are split
   // up for parallelization, even though the sequential AMR hierarhcy is
   // always aligned. (I'm not entirely sure why there occasionally are tiny
   // differences with the simple stencils I've tested so far ... if the parent is
   // in sync, then the edge alignment shouldn't matter. There may
   // be a 'mild' bug in PAMR, such as a round-off mismatch choosing which
   // parent points to interpolate from, etc.). 
   // And, as all my current MG experience suggests, functions need to
   // be *exactly* in sync for the vcycle to converge. Therefore,
   // for now turn off PAMR's bdy_width 'feature' ... this should have
   // a neglible effect on comminication speeds, as this is a 
   // 'd-2' dimensional affect on the pieces of grids communicated
   // during sync'ing, which communicates d-1 dimensional segments.
   //
   // Addendum: we need AMR_bdy_width>0 for PAMR_bdy_interp() to
   //           work ... recent experiments indicate that it may be
   //           OK simply to keep MG_bdy_width=0. If not, we will need
   //           to modify either the manner in which PAMR does 
   //           PAMR_bdy_interp (e.g., send an argument specify the boundary
   //           width, thus ignoring the internal parameter)
   //           or the manner in which AMR_bdy_width is interpetid (i.e.,
   //           keep this flag seperate from ghost_width, and when sync'ing
   //           only use ghostwidth)
   //
   // Addendum b: removed AMR_bdy_width/MG_bdy_width parameters from PAMR,
   //             and added argument to PAMR_bdy_interp.
   //--------------------------------------------------------------------------
   
   PAMR_set_interp_buffer(interp_buffer);

   AMRD_int_param(pfile,"rho_sp",rho_sp,1);
   AMRD_int_param(pfile,"rho_tm",rho_tm,1);
   for (i=1; i<PAMR_MAX_LEVS; i++) { rho_sp[i]=rho_sp[i-1]; rho_tm[i]=rho_tm[i-1]; }
   AMRD_int_param(pfile,"rho_sp_all",rho_sp,AMRD_max_lev);
   AMRD_int_param(pfile,"rho_tm_all",rho_tm,AMRD_max_lev);

   PAMR_set_rho(rho_sp,rho_tm,AMRD_max_lev);

   AMRD_real_param(pfile,"lambda",&lambda,1);
   PAMR_set_lambda(lambda);

   AMRD_int_param(pfile,"ghost_width",ghost_width,AMRD_dim);
   PAMR_set_ghost_width(ghost_width);

   AMRD_int_param(pfile,"gdm_grid_by_grid",&gdm_grid_by_grid,1);
   if (gdm_grid_by_grid) gdm_method|=PAMR_GDM_GRID_BY_GRID;

   AMRD_int_param(pfile,"gdm_align",&gdm_align,1);
   if (gdm_align) gdm_method|=PAMR_GDM_ALIGN;

   AMRD_int_param(pfile,"gdm_no_overlap",&gdm_no_overlap,1);
   if (gdm_no_overlap) gdm_method|=PAMR_GDM_NO_OVERLAP;

   PAMR_set_gdm(gdm_method);

   AMRD_int_param(pfile,"min_width",min_width,AMRD_dim);
   PAMR_set_min_width(min_width);

   AMRD_int_param(pfile,"min_mg_cwidth",min_cwidth,AMRD_dim);
   PAMR_set_MG_coarse_width(min_cwidth);

   num_even_vars_x0min=0; even_vars_x0min=0;
   AMRD_str_param_v(pfile,"even_vars_x0min",&even_vars_x0min,&num_even_vars_x0min); 
   if (num_even_vars_x0min>=AMRD_MAX_VARS) 
      AMRD_stop("num_even_vars_x0min is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_even_vars_x1min=0; even_vars_x1min=0;
   AMRD_str_param_v(pfile,"even_vars_x1min",&even_vars_x1min,&num_even_vars_x1min); 
   if (num_even_vars_x1min>=AMRD_MAX_VARS) 
      AMRD_stop("num_even_vars_x1min is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_even_vars_x2min=0; even_vars_x2min=0;
   AMRD_str_param_v(pfile,"even_vars_x2min",&even_vars_x2min,&num_even_vars_x2min); 
   if (num_even_vars_x2min>=AMRD_MAX_VARS) 
      AMRD_stop("num_even_vars_x2min is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_even_vars_x0max=0; even_vars_x0max=0;
   AMRD_str_param_v(pfile,"even_vars_x0max",&even_vars_x0max,&num_even_vars_x0max); 
   if (num_even_vars_x0max>=AMRD_MAX_VARS) 
      AMRD_stop("num_even_vars_x0max is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_even_vars_x1max=0; even_vars_x1max=0;
   AMRD_str_param_v(pfile,"even_vars_x1max",&even_vars_x1max,&num_even_vars_x1max); 
   if (num_even_vars_x1max>=AMRD_MAX_VARS) 
      AMRD_stop("num_even_vars_x1max is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_even_vars_x2max=0; even_vars_x2max=0;
   AMRD_str_param_v(pfile,"even_vars_x2max",&even_vars_x2max,&num_even_vars_x2max); 
   if (num_even_vars_x2max>=AMRD_MAX_VARS) 
      AMRD_stop("num_even_vars_x2max is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_odd_vars_x0min=0; odd_vars_x0min=0;
   AMRD_str_param_v(pfile,"odd_vars_x0min",&odd_vars_x0min,&num_odd_vars_x0min); 
   if (num_odd_vars_x0min>=AMRD_MAX_VARS) 
      AMRD_stop("num_odd_vars_x0min is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_odd_vars_x1min=0; odd_vars_x1min=0;
   AMRD_str_param_v(pfile,"odd_vars_x1min",&odd_vars_x1min,&num_odd_vars_x1min); 
   if (num_odd_vars_x1min>=AMRD_MAX_VARS) 
      AMRD_stop("num_odd_vars_x1min is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_odd_vars_x2min=0; odd_vars_x2min=0;
   AMRD_str_param_v(pfile,"odd_vars_x2min",&odd_vars_x2min,&num_odd_vars_x2min); 
   if (num_odd_vars_x2min>=AMRD_MAX_VARS) 
      AMRD_stop("num_odd_vars_x2min is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_odd_vars_x0max=0; odd_vars_x0max=0;
   AMRD_str_param_v(pfile,"odd_vars_x0max",&odd_vars_x0max,&num_odd_vars_x0max); 
   if (num_odd_vars_x0max>=AMRD_MAX_VARS) 
      AMRD_stop("num_odd_vars_x0max is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_odd_vars_x1max=0; odd_vars_x1max=0;
   AMRD_str_param_v(pfile,"odd_vars_x1max",&odd_vars_x1max,&num_odd_vars_x1max); 
   if (num_odd_vars_x1max>=AMRD_MAX_VARS) 
      AMRD_stop("num_odd_vars_x1max is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   num_odd_vars_x2max=0; odd_vars_x2max=0;
   AMRD_str_param_v(pfile,"odd_vars_x2max",&odd_vars_x2max,&num_odd_vars_x2max); 
   if (num_odd_vars_x2max>=AMRD_MAX_VARS) 
      AMRD_stop("num_odd_vars_x2max is greater than maximum allowed ... change AMRD_MAX_VARS and recompile",0); 

   //--------------------------------------------------------------------------
   // now define all of the variables
   //--------------------------------------------------------------------------
   i=0;
   int i_AMRD_fcs_vars;
   int i_AMRD_past_bdy_interp_vars;
   int i_AMRD_work_repop_vars;
   int i_AMRD_AMRH_work_inject_vars;
   int i_AMRD_inject_wavg_vars;
   i_AMRD_fcs_vars = 0;
   i_AMRD_past_bdy_interp_vars = 0;
   i_AMRD_work_repop_vars = 0;
   i_AMRD_AMRH_work_inject_vars = 0;
   i_AMRD_inject_wavg_vars=0;
   if (AMRD_echo_params && my_rank==0) printf("\nDefining hyperbolic variables:");
   if (AMRD_echo_params && my_rank==0) printf("\n==============================\n");
   if (AMRD_num_hyperbolic_vars==0) var=0; else var=AMRD_hyperbolic_vars[i];
   stage=0;
   if (!var)
   {
      stage++;
      if (AMRD_num_elliptic_vars==0) var=0; else var=AMRD_elliptic_vars[i];
      if (!var)
      {
         stage++;
         if (AMRD_num_AMRH_work_vars==0) var=0; else var=AMRD_AMRH_work_vars[i];
         if (!var)
         {
            stage++;
            if (AMRD_num_MGH_work_vars==0) var=0; else var=AMRD_MGH_work_vars[i];
            if  (!var)
            {
               stage++;
               if (AMRD_num_elliptic_vars_t0==0) var=0; else var=AMRD_elliptic_vars_t0[i];
               if (!var) AMRD_stop("Error ... no variables defined","");
            }
         }
      }
   }

   while(var)
   {
      amr_inject=PAMR_NO_INJECT;
      amr_interp=PAMR_NO_INTERP;
      amr_bdy_interp=PAMR_NO_INTERP;
      amr_sync=PAMR_NO_SYNC;
      amr_c_to_v = PAMR_C_TO_V_NO_TRANSFER;
      amr_v_to_c = PAMR_V_TO_C_NO_TRANSFER;
      mg_inject=PAMR_NO_INJECT;
      mg_interp=PAMR_NO_INTERP;
      mg_sync=PAMR_NO_SYNC;
      switch(stage)
      {
         case 0: // hyperbolic
            in_amrh=1; in_mgh=1; num_tl=AMRD_num_evo_tl;
            mg_noinj_to_amr=1;
            break;
         case 1: // elliptic
            in_amrh=1; in_mgh=1; num_tl=AMRD_num_evo_tl;
            mg_noinj_to_amr=0;
            break;
         case 2: // AMR work
            in_amrh=1; in_mgh=0; num_tl=1;
	    mg_noinj_to_amr=1;
            if (in_list(var,AMRD_AMRH_work_in_MGH_vars,AMRD_num_AMRH_work_in_MGH_vars)) in_mgh=1;
            break;
         case 3: // MG work
            in_amrh=0; in_mgh=1; num_tl=0;
            mg_noinj_to_amr=1;
            break;
         case 4: // elliptic_t0
            in_amrh=0; in_mgh=1; num_tl=0;
            mg_noinj_to_amr=0;
            break;
      }

      if (in_list(var,amr_sync_v,num_amr_sync_v)) amr_sync=PAMR_SYNC;
      if (in_list(var,mg_sync_v,num_mg_sync_v)) mg_sync=PAMR_SYNC;

      if (PAMR_var_type(var)==PAMR_VERTEX_CENTERED) {
	// Set default transfer properties for VC variables:
	regrid_transfer=PAMR_SECOND_ORDER;
	if (in_list(var,amr_interp2_v,num_amr_interp2_v)) amr_interp=amr_bdy_interp=PAMR_SECOND_ORDER;
	if (in_list(var,amr_interp4_v,num_amr_interp4_v)) amr_interp=amr_bdy_interp=PAMR_FOURTH_ORDER;
	if (in_list(var,amr_inject_v,num_amr_inject_v)) amr_inject=PAMR_STRAIGHT_INJECT;
	if (in_list(var,amr_hw_inject_v,num_amr_hw_inject_v)) amr_inject=PAMR_HW_RESTR;
	if (in_list(var,mg_hw_restr_v,num_mg_hw_restr_v)) mg_inject=PAMR_HW_RESTR;
	if (in_list(var,mg_fw_restr_v,num_mg_fw_restr_v)) mg_inject=PAMR_FW_RESTR;
	if (in_list(var,mg_hw_restr_bdy_si_v,num_mg_hw_restr_bdy_si_v)) mg_inject=PAMR_HW_RESTR_BDY_SI; // hwr in interior, straight inject on boundary
	if (in_list(var,amr_transfer2_v,num_amr_transfer2_v)) regrid_transfer=PAMR_SECOND_ORDER;
	if (in_list(var,amr_transfer4_v,num_amr_transfer4_v)) regrid_transfer=PAMR_FOURTH_ORDER;
	if (in_list(var,amr_bdy_interp2_v,num_amr_bdy_interp2_v)) amr_bdy_interp=PAMR_SECOND_ORDER;;
      }	else {
	// Set default transfer properties for CC variables:
	regrid_transfer=PAMR_MC;
	if (in_list(var,amr_interp1_v,num_amr_interp1_v)) amr_interp=amr_bdy_interp=PAMR_FIRST_ORDER_CONS;
	if (in_list(var,amr_interp2_v,num_amr_interp2_v)) amr_interp=amr_bdy_interp=PAMR_MC;
	if (in_list(var,amr_inject_v,num_amr_inject_v)) amr_inject=PAMR_NN_AVERAGE;
	if (in_list(var,AMRD_inject_wavg_vars,AMRD_num_inject_wavg_vars)) amr_inject=PAMR_NN_ADD;
	if (in_list(var,amr_transfer1_v,num_amr_transfer1_v)) regrid_transfer=PAMR_FIRST_ORDER_CONS;
	if (in_list(var,amr_transfer2_v,num_amr_transfer2_v)) regrid_transfer=PAMR_MC;
	// For variables which undergo interpolation, we use first order conservative
	// as the default, since something must be defined, and extensive is not a 
	// valid choice for the amr boundaries.
	if (amr_interp != PAMR_NO_INTERP) amr_bdy_interp=PAMR_FIRST_ORDER_CONS;
	if (in_list(var,amr_bdy_interp1_v,num_amr_bdy_interp1_v)) amr_bdy_interp=PAMR_FIRST_ORDER_CONS;
	if (in_list(var,amr_bdy_interp2_v,num_amr_bdy_interp2_v)) amr_bdy_interp=PAMR_MC;
      }

      if (in_list(var,amr_c_to_v2_v,num_amr_c_to_v2_v)) amr_c_to_v=PAMR_C_TO_V_SECOND_ORDER;
      if (in_list(var,amr_v_to_c2_v,num_amr_v_to_c2_v)) amr_v_to_c=PAMR_V_TO_C_SECOND_ORDER;

      phys_bdy_type[0]=phys_bdy_type[1]=phys_bdy_type[2]=phys_bdy_type[3]=phys_bdy_type[4]=phys_bdy_type[5]=PAMR_UNKNOWN;
      if (in_list(var,even_vars_x0min,num_even_vars_x0min)) phys_bdy_type[0]=PAMR_EVEN;
      if (in_list(var,even_vars_x1min,num_even_vars_x1min)) phys_bdy_type[2]=PAMR_EVEN;
      if (in_list(var,even_vars_x2min,num_even_vars_x2min)) phys_bdy_type[4]=PAMR_EVEN;
      if (in_list(var,even_vars_x0max,num_even_vars_x0max)) phys_bdy_type[1]=PAMR_EVEN;
      if (in_list(var,even_vars_x1max,num_even_vars_x1max)) phys_bdy_type[3]=PAMR_EVEN;
      if (in_list(var,even_vars_x2max,num_even_vars_x2max)) phys_bdy_type[5]=PAMR_EVEN;
      if (in_list(var,odd_vars_x0min,num_odd_vars_x0min)) phys_bdy_type[0]=PAMR_ODD;
      if (in_list(var,odd_vars_x1min,num_odd_vars_x1min)) phys_bdy_type[2]=PAMR_ODD;
      if (in_list(var,odd_vars_x2min,num_odd_vars_x2min)) phys_bdy_type[4]=PAMR_ODD;
      if (in_list(var,odd_vars_x0max,num_odd_vars_x0max)) phys_bdy_type[1]=PAMR_ODD;
      if (in_list(var,odd_vars_x1max,num_odd_vars_x1max)) phys_bdy_type[3]=PAMR_ODD;
      if (in_list(var,odd_vars_x2max,num_odd_vars_x2max)) phys_bdy_type[5]=PAMR_ODD;

      sgfn=PAMR_def_var_full(var,in_amrh,in_mgh,num_tl,amr_inject,amr_interp,
                             amr_bdy_interp,amr_sync,mg_inject,mg_interp,
                             mg_sync,mg_noinj_to_amr,regrid_transfer,amr_c_to_v,
			     amr_v_to_c,phys_bdy_type);

/*       if (AMRD_echo_params && my_rank==0) printf("%s: in_amrh=%i, in_mgh=%i, num_tl=%i, amr_inject=%i\n",var,in_amrh,in_mgh,num_tl,amr_inject); */
/*       if (AMRD_echo_params && my_rank==0) printf("\tamr_interp=%i, amr_bdy_interp=%i, amr_sync=%i, mg_inject=%i\n",amr_interp,amr_bdy_interp,amr_sync,mg_inject); */
/*       if (AMRD_echo_params && my_rank==0) printf("\tmg_interp=%i, mg_sync=%i, mg_noinj_to_amr=%i, regrid_transfer=%i\n",mg_interp,mg_sync,mg_noinj_to_amr,regrid_transfer); */
/*       if (AMRD_echo_params && my_rank==0) printf("\tphys_bdy_type=[%i,%i,%i,%i,%i,%i]\n\n", */
/*          phys_bdy_type[0],phys_bdy_type[1],phys_bdy_type[2],phys_bdy_type[3],phys_bdy_type[4],phys_bdy_type[5]); */

      if(in_list(var,AMRD_inject_wavg_vars,AMRD_num_inject_wavg_vars)){
        if (PAMR_var_type(var)==PAMR_VERTEX_CENTERED) AMRD_stop("ERROR! AMRD does not support inject_wavg for vertex-centered variables",0); 
      	if(stage==2) AMRD_inject_wavg_gfn[i_AMRD_inject_wavg_vars] = sgfn;
      	else if(stage==0 && num_tl>=2) AMRD_inject_wavg_gfn[i_AMRD_inject_wavg_vars] = sgfn+1; //wavg inject time level 2 for hyperbolic
	else AMRD_stop("ERROR! Can only wavg inject hyberbolic or work variables",0); 
        i_AMRD_inject_wavg_vars++;
      }

      if (stage==2) {
	// This will be a list of work variables which will need to be 
	// repopulated when doing excision.  This was introduced for
	// hydrodynamics. (Branson)
	if (in_list(var,AMRD_work_repop_vars,AMRD_num_work_repop_vars)) {
	  AMRD_work_repop_gfn[i_AMRD_work_repop_vars] = sgfn;
	  i_AMRD_work_repop_vars++;
	}

	// This will be a list of work variables which will need to be 
	// injected along with the evolution variables.  This was introduced for
	// hydrodynamics. (Branson)
	if (in_list(var,AMRD_AMRH_work_inject_vars,AMRD_num_AMRH_work_inject_vars)) {
	  AMRD_AMRH_work_inject_gfn[i_AMRD_AMRH_work_inject_vars] = sgfn;
	  AMRD_AMRH_work_inject_op[i_AMRD_AMRH_work_inject_vars] = amr_inject;
	  i_AMRD_AMRH_work_inject_vars++;
	}
      }

      phys_bdy_type[0]=phys_bdy_type[1]=phys_bdy_type[2]=phys_bdy_type[3]=phys_bdy_type[4]=phys_bdy_type[5]=PAMR_UNKNOWN;
      if (stage==0)
      {
         for (j=0; j<num_tl; j++) AMRD_AMR_f_gfn[j][i]=sgfn+j;

	 if (in_list(var,AMRD_past_bdy_interp_vars,AMRD_num_past_bdy_interp_vars)) {
	   for (j=0; j<num_tl; j++) AMRD_past_bdy_interp_vars_gfn[j][i_AMRD_past_bdy_interp_vars] = sgfn+j;
	   AMRD_past_bdy_interp_op[i_AMRD_past_bdy_interp_vars] = amr_bdy_interp;
	   i_AMRD_past_bdy_interp_vars++;
	 }

	 if (in_list(var,AMRD_fc_vars,AMRD_num_fc_vars)) {
	   // You need to create the appropriate fcs variable here.
            in_amrh=1; 
            num_tl=1;
            mg_noinj_to_amr=1;
	    in_mgh=0;
	    strcpy(buffer,var);
	    AMRD_append_tag(buffer,"_fcs",AMRD_v_tag,AMRD_c_tag);
	    // Branson commented out.
	    //            if (AMRD_echo_params && my_rank==0) printf("%s:\n\n",buffer);
	    amr_inject=PAMR_NO_INJECT;
	    amr_interp=PAMR_NO_INTERP;
	    amr_bdy_interp=PAMR_NO_INTERP;
            mg_inject=PAMR_NO_INJECT;
            mg_interp=PAMR_NO_INTERP;
	    mg_sync=PAMR_NO_SYNC;
	    amr_c_to_v=PAMR_C_TO_V_NO_TRANSFER;
	    amr_v_to_c=PAMR_V_TO_C_NO_TRANSFER;
	    regrid_transfer=PAMR_NO_INTERP;

            AMRD_fcs_gfn[i_AMRD_fcs_vars]=PAMR_def_var_full(buffer,in_amrh,in_mgh,num_tl,amr_inject,amr_interp,
							    amr_bdy_interp,amr_sync,mg_inject,mg_interp,
							    mg_sync,mg_noinj_to_amr,regrid_transfer,amr_c_to_v,
							    amr_v_to_c,phys_bdy_type);

	    // Note that this little counter is being used to index the gfn array.
	    i_AMRD_fcs_vars++;
         }

         if (stage==0 && in_list(var,AMRD_TRE_vars,AMRD_num_TRE_vars))
         {
            amr_inject=PAMR_NO_INJECT;
            amr_interp=PAMR_NO_INTERP;
            amr_bdy_interp=PAMR_NO_INTERP;
            amr_sync=PAMR_NO_SYNC;
            mg_inject=PAMR_NO_INJECT;
            mg_interp=PAMR_NO_INTERP;
	    amr_c_to_v=PAMR_C_TO_V_NO_TRANSFER;
	    amr_v_to_c=PAMR_V_TO_C_NO_TRANSFER;
            mg_sync=PAMR_NO_SYNC;
	    if (PAMR_var_type(var)==PAMR_VERTEX_CENTERED) {
	      regrid_transfer=PAMR_SECOND_ORDER; // important for check-pointing between different # of nodes.
              if (AMRD_using_cc_tre==1) AMRD_stop("Error ... at present all TRE_vars must be of the same type (cell/vertex centered)","");
	      AMRD_using_cc_tre = 0;
	    } else {
	      regrid_transfer=PAMR_MC; // important for check-pointing between different # of nodes.
              if (AMRD_using_cc_tre==0) AMRD_stop("Error ... at present all TRE_vars must be of the same type (cell/vertex centered)","");
	      AMRD_using_cc_tre = 1;
	    }
            in_amrh=1; 
            num_tl=1;
            mg_noinj_to_amr=1;
            //-----------------------------------------------------------------
            // we use the first tre variable to store a norm
            // of the MG TREs, to calculate initial data via method 1
            //-----------------------------------------------------------------
            if (AMRD_num_f_tre_vars==0) in_mgh=1; else in_mgh=0;
	    strcpy(buffer,var); 
	    AMRD_append_tag(buffer,"_tre",AMRD_v_tag,AMRD_c_tag);
            AMRD_f_tre_vars[AMRD_num_f_tre_vars]=strdup(buffer);

            AMRD_f_tre_fn[AMRD_num_f_tre_vars]=i;
            AMRD_f_tre_gfn[AMRD_num_f_tre_vars++]=PAMR_def_var_full(buffer,in_amrh,in_mgh,num_tl,amr_inject,amr_interp,
                                        amr_bdy_interp,amr_sync,mg_inject,mg_interp,
                                        mg_sync,mg_noinj_to_amr,regrid_transfer,amr_c_to_v,amr_v_to_c,phys_bdy_type);
            if (in_mgh) AMRD_mg_tre_gfn=AMRD_f_tre_gfn[AMRD_num_f_tre_vars-1]+1;

            if (AMRD_echo_params && my_rank==0) printf("%s:\n\n",buffer);
         }

      }

      if (stage==1 || stage==4)
      {
         if (stage==1) offset=0; else offset=AMRD_num_elliptic_vars;
         AMRD_MG_f_gfn[i+offset]=PAMR_get_gfn(var,PAMR_MGH,0);

         if (stage==1) // no AMR correspondence for _t0 vars
            for (k=1; k<=AMRD_num_evo_tl; k++) 
               AMRD_AMR_mgf_gfn[k-1][i]=AMRD_MG_f_gfn[i]-AMRD_num_evo_tl+k-1;
	 strcpy(buffer,var);
	 AMRD_append_tag(buffer,"_res",AMRD_v_tag,AMRD_c_tag);
         amr_inject=PAMR_NO_INJECT;
         amr_interp=PAMR_NO_INTERP;
         amr_bdy_interp=PAMR_NO_INTERP;
         amr_sync=PAMR_NO_SYNC;
         mg_inject=PAMR_HW_RESTR;
         mg_interp=PAMR_NO_INTERP;
         if (in_list(buffer,mg_fw_restr_v,num_mg_fw_restr_v)) mg_inject=PAMR_FW_RESTR;
         if (in_list(var,mg_interp2_v,num_mg_interp2_v)) mg_interp=PAMR_SECOND_ORDER;
         mg_sync=PAMR_NO_SYNC;
         regrid_transfer=PAMR_NO_TRANSFER;
         in_amrh=0; in_mgh=1; num_tl=0;
         mg_noinj_to_amr=0;

         AMRD_MG_res_gfn[i+offset]=
            PAMR_def_var_full(buffer,in_amrh,in_mgh,num_tl,amr_inject,amr_interp,
                              amr_bdy_interp,amr_sync,mg_inject,mg_interp,
                              mg_sync,mg_noinj_to_amr,regrid_transfer,amr_c_to_v,amr_v_to_c,phys_bdy_type);
         if (AMRD_echo_params && my_rank==0) printf("%s: mg_inject=%i, mg_interp (for %s !)=%i\n\n",buffer,mg_inject,var,mg_interp);

	 strcpy(buffer,var);
	 AMRD_append_tag(buffer,"_lop",AMRD_v_tag,AMRD_c_tag);
         mg_inject=PAMR_NO_INJECT;
         mg_interp=PAMR_NO_INTERP;
         AMRD_MG_lop_gfn[i+offset]=
            PAMR_def_var_full(buffer,in_amrh,in_mgh,num_tl,amr_inject,amr_interp,
                              amr_bdy_interp,amr_sync,mg_inject,mg_interp,
                              mg_sync,mg_noinj_to_amr,regrid_transfer,amr_c_to_v,amr_v_to_c,
			      phys_bdy_type);
         if (AMRD_echo_params && my_rank==0) printf("%s:\n\n",buffer);

	 strcpy(buffer,var);
	 AMRD_append_tag(buffer,"_rhs",AMRD_v_tag,AMRD_c_tag);
         mg_inject=PAMR_NO_INJECT;
         mg_interp=PAMR_NO_INTERP;
         AMRD_MG_rhs_gfn[i+offset]=
            PAMR_def_var_full(buffer,in_amrh,in_mgh,num_tl,amr_inject,amr_interp,
                              amr_bdy_interp,amr_sync,mg_inject,mg_interp,
                              mg_sync,mg_noinj_to_amr,regrid_transfer,
			      amr_c_to_v,amr_v_to_c,phys_bdy_type);
         if (AMRD_echo_params && my_rank==0) printf("%s:\n\n",buffer);

         if (stage==1)
         {
	   strcpy(buffer,var);
	   AMRD_append_tag(buffer,"_brs",AMRD_v_tag,AMRD_c_tag);
	    in_amrh=1; in_mgh=0; num_tl=1;
	    amr_c_to_v = PAMR_C_TO_V_NO_TRANSFER;
	    amr_v_to_c = PAMR_V_TO_C_NO_TRANSFER;
            AMRD_MG_brs_gfn[i]=
               PAMR_def_var_full(buffer,in_amrh,in_mgh,num_tl,amr_inject,amr_interp,
				 amr_bdy_interp,amr_sync,mg_inject,mg_interp,
				 mg_sync,mg_noinj_to_amr,regrid_transfer,amr_c_to_v,
				 amr_v_to_c,phys_bdy_type);
            if (AMRD_echo_params && my_rank==0) printf("%s:\n\n",buffer);

	    strcpy(buffer,var);
	    AMRD_append_tag(buffer,"_extrap_tm1",AMRD_v_tag,AMRD_c_tag);
            amr_inject=PAMR_STRAIGHT_INJECT;
            amr_interp=PAMR_FOURTH_ORDER;
            if (in_list(buffer,amr_interp2_v,num_amr_interp2_v)) amr_interp=amr_bdy_interp=PAMR_SECOND_ORDER;
            regrid_transfer=PAMR_FOURTH_ORDER;
            if (in_list(buffer,amr_transfer2_v,num_amr_transfer2_v)) regrid_transfer=PAMR_SECOND_ORDER;
            in_amrh=1; in_mgh=0; num_tl=1;
            mg_noinj_to_amr=1;
            AMRD_AMR_mgf_extrap_tm1_gfn[i]=
               PAMR_def_var_full(buffer,in_amrh,in_mgh,num_tl,amr_inject,amr_interp,
				 amr_bdy_interp,amr_sync,mg_inject,mg_interp,
				 mg_sync,mg_noinj_to_amr,regrid_transfer,amr_c_to_v,amr_v_to_c,phys_bdy_type);
            if (AMRD_echo_params && my_rank==0) printf("%s (in AMRH): amr_interp=%i, regrid_transfer=%i\n\n",buffer,amr_interp,regrid_transfer);

	    strcpy(buffer,var);
	    AMRD_append_tag(buffer,"_extrap_tm2",AMRD_v_tag,AMRD_c_tag);
            amr_inject=PAMR_STRAIGHT_INJECT;
            amr_interp=PAMR_FOURTH_ORDER;
            if (in_list(buffer,amr_interp2_v,num_amr_interp2_v)) amr_interp=amr_bdy_interp=PAMR_SECOND_ORDER;
            regrid_transfer=PAMR_FOURTH_ORDER;
            if (in_list(buffer,amr_transfer2_v,num_amr_transfer2_v)) regrid_transfer=PAMR_SECOND_ORDER;
            AMRD_AMR_mgf_extrap_tm2_gfn[i]=
               PAMR_def_var_full(buffer,in_amrh,in_mgh,num_tl,amr_inject,amr_interp,
                              amr_bdy_interp,amr_sync,mg_inject,mg_interp,
                              mg_sync,mg_noinj_to_amr,regrid_transfer,amr_c_to_v,amr_v_to_c,phys_bdy_type);
            if (AMRD_echo_params && my_rank==0) printf("%s (in AMRH): amr_interp=%i, regrid_transfer=%i\n\n",buffer,amr_interp,regrid_transfer);
         } 
      }

      i++;
      switch(stage)
      {
         case 0: if (i>=AMRD_num_hyperbolic_vars) var=0; else var=AMRD_hyperbolic_vars[i]; break;
         case 1: if (i>=AMRD_num_elliptic_vars) var=0; else var=AMRD_elliptic_vars[i]; break;
         case 2: if (i>=AMRD_num_AMRH_work_vars) var=0; else var=AMRD_AMRH_work_vars[i]; break;
         case 3: if (i>=AMRD_num_MGH_work_vars) var=0; else var=AMRD_MGH_work_vars[i]; break;
         case 4: if (i>=AMRD_num_elliptic_vars_t0) var=0; else var=AMRD_elliptic_vars_t0[i]; break;
      }

      while(!var && stage<5) 
      {
         stage++;
         i=0;
         switch(stage)
         {
            case 1: if (i>=AMRD_num_elliptic_vars) var=0; else var=AMRD_elliptic_vars[i]; 
                    if (AMRD_echo_params && my_rank==0) printf("\nDefining elliptic variables:");
                    if (AMRD_echo_params && my_rank==0) printf("\n===========================\n");
                    break;
            case 2: if (i>=AMRD_num_AMRH_work_vars) var=0; else var=AMRD_AMRH_work_vars[i]; 
                    if (AMRD_echo_params && my_rank==0) printf("\nDefining AMR work variables (%i):",AMRD_num_AMRH_work_vars);
                    if (AMRD_echo_params && my_rank==0) printf("\n================================\n");
                    break;
            case 3: if (i>=AMRD_num_MGH_work_vars) var=0; else  var=AMRD_MGH_work_vars[i]; 
                    if (AMRD_echo_params && my_rank==0) printf("\nDefining MG work variables (%i):",AMRD_num_MGH_work_vars);
                    if (AMRD_echo_params && my_rank==0) printf("\n===============================\n");
                    break;
            case 4: if (i>=AMRD_num_elliptic_vars_t0) var=0; else var=AMRD_elliptic_vars_t0[i]; 
                    if (AMRD_echo_params && my_rank==0) printf("\nDefining elliptic_t0 variables:");
                    if (AMRD_echo_params && my_rank==0) printf("\n===============================\n");
                    break;
         }
      }
   }

   if (AMRD_echo_params && my_rank==0) printf("\nDefining misc. variables:");
   if (AMRD_echo_params && my_rank==0) printf("\n=========================\n");

   amr_inject=PAMR_NO_INJECT;
   amr_interp=PAMR_NO_INTERP;
   amr_bdy_interp=PAMR_NO_INTERP;
   amr_sync=PAMR_NO_SYNC;
   mg_inject=PAMR_NO_INJECT;
   mg_interp=PAMR_NO_INTERP;
   mg_sync=PAMR_NO_SYNC;
   regrid_transfer=PAMR_NO_TRANSFER;
   in_amrh=1; in_mgh=1; num_tl=1;
   mg_noinj_to_amr=1;
   // We won't allow these masks to undergo c_to_v and v_to_c transformations.
   amr_c_to_v = PAMR_C_TO_V_NO_TRANSFER;
   amr_v_to_c = PAMR_V_TO_C_NO_TRANSFER;

   strcpy(buffer,"cmask");
   if (AMRD_v_tag) strcat(buffer,AMRD_v_tag);
   AMRD_cmask_gfn=PAMR_def_var_full(buffer,in_amrh,in_mgh,num_tl,amr_inject,amr_interp,
                               amr_bdy_interp,amr_sync,mg_inject,mg_interp,
                               mg_sync,mg_noinj_to_amr,regrid_transfer,amr_c_to_v,amr_v_to_c,phys_bdy_type);
   AMRD_cmask_mg_gfn=AMRD_cmask_gfn+1;
   if (AMRD_echo_params && my_rank==0) printf("%s:\n\n",buffer);

   if (!AMRD_v_tag && !AMRD_c_tag) strcpy(buffer,"__cmask_c"); // no cc-var's, but defining dummy cmask for simplicity
   else strcpy(buffer,"cmask");
   if (AMRD_c_tag) strcat(buffer,AMRD_c_tag);
   AMRD_cmask_c_gfn=PAMR_def_var_full(buffer,in_amrh,in_mgh,num_tl,amr_inject,amr_interp,
		   	          amr_bdy_interp,amr_sync,mg_inject,mg_interp,
                                  mg_sync,mg_noinj_to_amr,regrid_transfer,PAMR_C_TO_V_NO_TRANSFER,PAMR_V_TO_C_NO_TRANSFER,phys_bdy_type);
   if (AMRD_echo_params && my_rank==0) printf("%s:\n\n",buffer);

   // The mask for flux correction shall have the same properties as the other masks.
   if (AMRD_c_tag || AMRD_v_tag)
   {
      strcpy(buffer,"fc_mask");
      if (AMRD_c_tag) strcat(buffer,AMRD_c_tag);
      AMRD_fc_mask_gfn=PAMR_def_var_full(buffer,in_amrh,in_mgh,num_tl,amr_inject,amr_interp,
				amr_bdy_interp,amr_sync,mg_inject,mg_interp,
				mg_sync,mg_noinj_to_amr,regrid_transfer,PAMR_C_TO_V_NO_TRANSFER,PAMR_V_TO_C_NO_TRANSFER,phys_bdy_type);
      if (AMRD_echo_params && my_rank==0) printf("%s:\n\n",buffer);
   }
   else AMRD_fc_mask_gfn=0;

   AMRD_chr_gfn=0; AMRD_chr_mg_gfn=0;
   AMRD_chr_c_gfn=0; AMRD_chr_mg_c_gfn=0;
   if (AMRD_do_ex)
   {
      chr_v=chr_name_v; // not saving the names in global storage, only the gfn's ... can do so later if needs be
      chr_c=chr_name_c;
      if (AMRD_v_tag) sprintf(chr_v,"chr%s",AMRD_v_tag);
      else strcpy(chr_v,"chr");

      if (AMRD_c_tag) sprintf(chr_c,"chr%s",AMRD_c_tag); 
      else if (AMRD_v_tag && !AMRD_c_tag) strcpy(chr_c,"chr");
      else strcpy(chr_c,"__AMRD_chr_c"); // no cc-var's, but defining dummy AMRD_chr_c for simplicity

      AMRD_chr_gfn=PAMR_def_var_full(chr_v,in_amrh,in_mgh,num_tl,amr_inject,amr_interp,
                            amr_bdy_interp,amr_sync,mg_inject,mg_interp,
                            mg_sync,mg_noinj_to_amr,regrid_transfer,amr_c_to_v,amr_v_to_c,phys_bdy_type);
      AMRD_chr_mg_gfn=AMRD_chr_gfn+1;
      if (AMRD_echo_params && my_rank==0) printf("%s:\n\n",chr_v);

      AMRD_chr_c_gfn=PAMR_def_var_full(chr_c,in_amrh,in_mgh,num_tl,amr_inject,amr_interp,
                               amr_bdy_interp,amr_sync,mg_inject,mg_interp,
                               mg_sync,mg_noinj_to_amr,regrid_transfer,amr_c_to_v,amr_v_to_c,phys_bdy_type);
      AMRD_chr_mg_c_gfn=AMRD_chr_c_gfn+1;
      if (AMRD_echo_params && my_rank==0) printf("%s:\n\n",chr_c);
   }


   AMRD_wavg_c_gfn=0; AMRD_wavg_mg_c_gfn=0;
   if (AMRD_num_inject_wavg_vars)
   {
      wavg_c=wavg_name_c; // not saving the names in global storage, only the gfn's ... can do so later if needs be

      if (AMRD_c_tag) sprintf(wavg_c,"wavg%s",AMRD_c_tag); 
      else if (AMRD_v_tag && !AMRD_c_tag) strcpy(wavg_c,"wavg");
      else strcpy(wavg_c,"__AMRD_wavg_c"); // no cc-var's, but defining dummy AMRD_wavg_c for simplicity


      AMRD_wavg_c_gfn=PAMR_def_var_full(wavg_c,in_amrh,in_mgh,num_tl,amr_inject,amr_interp,
                               amr_bdy_interp,amr_sync,mg_inject,mg_interp,
                               mg_sync,mg_noinj_to_amr,regrid_transfer,amr_c_to_v,amr_v_to_c,phys_bdy_type);
      AMRD_wavg_mg_c_gfn=AMRD_wavg_c_gfn+1;
      if (AMRD_echo_params && my_rank==0) printf("%s:\n\n",wavg_c);
   }

   in_mgh=0;

   strcpy(buffer,"AMRD_tre");
   if (AMRD_v_tag) strcat(buffer,AMRD_v_tag);
   AMRD_tre_gfn=PAMR_def_var_full(buffer,in_amrh,in_mgh,num_tl,amr_inject,amr_interp,
                             amr_bdy_interp,amr_sync,mg_inject,mg_interp,
                             mg_sync,mg_noinj_to_amr,regrid_transfer,amr_c_to_v,amr_v_to_c,phys_bdy_type);
   if (AMRD_echo_params && my_rank==0) printf("%s:\n\n",buffer);

   if (!AMRD_v_tag && !AMRD_c_tag) strcpy(buffer,"__AMRD_tre_c"); // no cc-var's, but defining dummy AMRD_tre_c for simplicity
   else strcpy(buffer,"AMRD_tre");
   if (AMRD_c_tag) strcat(buffer,AMRD_c_tag);
   AMRD_tre_c_gfn=PAMR_def_var_full(buffer,in_amrh,in_mgh,num_tl,amr_inject,amr_interp,
				    amr_bdy_interp,amr_sync,mg_inject,mg_interp,
				    mg_sync,mg_noinj_to_amr,regrid_transfer,amr_c_to_v,amr_v_to_c,phys_bdy_type);
   if (AMRD_echo_params && my_rank==0) printf("%s:\n\n",buffer);

   strcpy(buffer,"amrd_w1");
   if (AMRD_v_tag) strcat(buffer,AMRD_v_tag);
   AMRD_w1_gfn=PAMR_def_var_full(buffer,in_amrh,in_mgh,num_tl,amr_inject,amr_interp,
                                 amr_bdy_interp,amr_sync,mg_inject,mg_interp,
                                 mg_sync,mg_noinj_to_amr,regrid_transfer,amr_c_to_v,amr_v_to_c,phys_bdy_type);
   if (AMRD_echo_params && my_rank==0) printf("%s:\n\n",buffer);

   if (!AMRD_v_tag && !AMRD_c_tag) strcpy(buffer,"__AMRD_w1_c"); // no cc-var's, but defining dummy AMRD_w1_c for simplicity
   else strcpy(buffer,"amrd_w1");
   if (AMRD_c_tag) strcat(buffer,AMRD_c_tag);
   AMRD_w1_c_gfn=PAMR_def_var_full(buffer,in_amrh,in_mgh,num_tl,amr_inject,amr_interp,
				   amr_bdy_interp,amr_sync,mg_inject,mg_interp,
				   mg_sync,mg_noinj_to_amr,regrid_transfer,amr_c_to_v,amr_v_to_c,phys_bdy_type);
   if (AMRD_echo_params && my_rank==0) printf("%s:\n\n",buffer);

   for (i=0; i<AMRD_num_ex_repop1_vars; i++) AMRD_ex_repop1_sgfn[i]=PAMR_get_gfn(AMRD_ex_repop1_vars[i],PAMR_AMRH,1);
   for (i=0; i<AMRD_num_ex_repop2_vars; i++) AMRD_ex_repop2_sgfn[i]=PAMR_get_gfn(AMRD_ex_repop2_vars[i],PAMR_AMRH,1);
   for (i=0; i<AMRD_num_ex_repop3_vars; i++) AMRD_ex_repop3_sgfn[i]=PAMR_get_gfn(AMRD_ex_repop3_vars[i],PAMR_AMRH,1);
   for (i=0; i<AMRD_num_ex_repop4_vars; i++) AMRD_ex_repop4_sgfn[i]=PAMR_get_gfn(AMRD_ex_repop4_vars[i],PAMR_AMRH,1);

   for (i=0; i<AMRD_num_tnp1_diss_vars; i++)
   {
      AMRD_tnp1_diss_gfn[i]=PAMR_get_gfn(AMRD_tnp1_diss_vars[i],PAMR_AMRH,1);
      PAMR_get_var_attribs(AMRD_tnp1_diss_vars[i],inum,inum,inum,inum,
                           inum,inum,inum,inum,inum,inum,inum,inum,inum,inum,&AMRD_tnp1_diss_pbt[i][0]);
   }
   if (AMRD_num_evo_tl>2 && AMRD_diss_all_past) mod=AMRD_num_evo_tl-1; else mod=1;

   for (i=0; i<(AMRD_num_tn_diss_vars*mod); i++) 
   {
      if (AMRD_num_evo_tl<2 || (  !(in_list(AMRD_tn_diss_vars[i%AMRD_num_tn_diss_vars],AMRD_hyperbolic_vars,AMRD_num_hyperbolic_vars))
                                &&!(in_list(AMRD_tn_diss_vars[i%AMRD_num_tn_diss_vars],AMRD_elliptic_vars,AMRD_num_elliptic_vars)) ) )
      {
         printf("specified <%s> as a tn_diss_var when either\nAMRD_num_evo_tl<2 or variable is not a hyperbolic/elliptic variable\n",
                AMRD_tn_diss_vars[i%AMRD_num_tn_diss_vars]);
         AMRD_stop("","");
      }
      AMRD_tn_diss_gfn[i]=PAMR_get_gfn(AMRD_tn_diss_vars[i%AMRD_num_tn_diss_vars],PAMR_AMRH,2+(int)(i/AMRD_num_tn_diss_vars));
      PAMR_get_var_attribs(AMRD_tn_diss_vars[i%AMRD_num_tn_diss_vars],inum,inum,inum,inum,
                           inum,inum,inum,inum,inum,inum,inum,inum,inum,inum,&AMRD_tn_diss_pbt[i][0]);
   }
   AMRD_num_tn_diss_vars*=mod;

   for (i=0; i<AMRD_num_tnp1_liipb_vars; i++) AMRD_tnp1_liipb_gfn[i]=PAMR_get_gfn(AMRD_tnp1_liipb_vars[i],PAMR_AMRH,1);
   for (i=0; i<AMRD_num_tnp1_liiab_vars; i++) AMRD_tnp1_liiab_gfn[i]=PAMR_get_gfn(AMRD_tnp1_liiab_vars[i],PAMR_AMRH,1);
   for (i=0; i<AMRD_num_tnp1_liibb_vars; i++) AMRD_tnp1_liibb_gfn[i]=PAMR_get_gfn(AMRD_tnp1_liibb_vars[i],PAMR_AMRH,1);
   for (i=0,ii=0; i<AMRD_num_rg_diss_vars; i++) 
   {
      AMRD_rg_diss_gfn[ii]=PAMR_get_gfn(AMRD_rg_diss_vars[i],PAMR_AMRH,1);
      PAMR_get_var_attribs(AMRD_rg_diss_vars[i],inum,inum,&num_tl,inum,
                           inum,inum,inum,inum,inum,inum,inum,inum,inum,inum,&AMRD_rg_diss_pbt[ii][0]);
      ii++;

      //-----------------------------------------------------------------------
      // PG compiler seems to have a bug ... the following statement
      //
      // rg_diss_gfn[ii++]=rg_diss_gfn[ii-1]+1;
      //
      // is not the same as the line below ... it seems to increment ii
      // before make the assignment
      //-----------------------------------------------------------------------
      while(--num_tl)
      {
         AMRD_rg_diss_gfn[ii]=AMRD_rg_diss_gfn[ii-1]+1; 
         for (j=0; j<2*PAMR_MAX_DIM; j++) AMRD_rg_diss_pbt[ii][j]=AMRD_rg_diss_pbt[ii-1][j];
         ii++;
      }
   }
   AMRD_tnum_rg_diss_vars=ii;
   AMRD_num_MG_inject_wavg_vars=0;
   for (i=0; i<AMRD_num_interp_AMR_bdy_vars; i++) AMRD_interp_AMR_bdy_f_gfn[i]=PAMR_get_gfn(AMRD_interp_AMR_bdy_vars[i],PAMR_AMRH,1);
   for (i=0; i<AMRD_num_MG_cnst_data_vars; i++){
	AMRD_MG_cnst_data_gfn[i]=PAMR_get_gfn(AMRD_MG_cnst_data_vars[i],PAMR_MGH,0);
        if( PAMR_var_type(AMRD_MG_cnst_data_vars[i])==PAMR_CELL_CENTERED && AMRD_num_inject_wavg_vars) {
        	AMRD_MG_inject_wavg_gfn[AMRD_num_MG_inject_wavg_vars] = AMRD_MG_cnst_data_gfn[i];
		AMRD_num_MG_inject_wavg_vars++; 
	}
   }
   /*
   Get the gfns for free data variables which will be synced/interpolated/injected
   after the set free data hook function is called.  Note that we are *not* saving this 
   list in the check-point file since this should be needed only at the initial time.
   */
   for (i=0; i<AMRD_num_free_data_vars; i++){
        AMRD_free_data_gfn[i]=PAMR_get_gfn(AMRD_free_data_vars[i],PAMR_AMRH,1);
        //Free data will be set on the AMRD_ic_n time-level for variables with multiple time levels.
        if (in_list(AMRD_free_data_vars[i],AMRD_hyperbolic_vars,AMRD_num_hyperbolic_vars) ||
            in_list(AMRD_free_data_vars[i],AMRD_elliptic_vars,AMRD_num_elliptic_vars)   ) {
             AMRD_free_data_gfn[i] += (AMRD_ic_n -1); 
        }
   }
   for (i=0; i<AMRD_num_elliptic_vars_t0; i++) AMRD_elliptic_vars_t0_gfn[i]=PAMR_get_gfn(AMRD_elliptic_vars_t0[i],PAMR_MGH,0);

   //==========================================================================
   // If max_lev>1, then assume AMR, and initialize base+shadow level;
   // else, only initialize base level.
   //==========================================================================
   init_lev[0]=1;
   init_lev[1]=2;
   for (i=0; i<(AMRD_dim*4*init_depth); i++) init_bbox[i]=0;
   for (i=0; i<AMRD_dim; i++)
   {
      init_bbox[2*i]=init_bbox[2*(i+AMRD_dim)]=AMRD_base_bbox[2*i];
      init_bbox[2*i+1]=init_bbox[2*(i+AMRD_dim)+1]=AMRD_base_bbox[2*i+1];
   }
   if (AMRD_max_lev>1)
   {
      if (AMRD_echo_params && my_rank==0) printf("\nAMR on: initializing base level plus shadow, up to level %i\n\n",init_depth);
      lev_num=2;
      num=2;
      while (lev_num<init_depth)
      {
         lev_num++;
         num++;
         init_lev[num-1]=lev_num;
         sprintf(buffer,"init_bbox_%i",lev_num);
         init_bbox[(num-1)*2*AMRD_dim]=-1e6;
         AMRD_real_param(pfile,buffer,&init_bbox[(num-1)*2*AMRD_dim],2*AMRD_dim);
         if (init_bbox[(num-1)*2*AMRD_dim]==-1e6) AMRD_stop("init_depth>2 specified, but following variable not found:",buffer);
         num++;
         init_lev[num-1]=lev_num;
         sprintf(buffer,"init_bbox_b_%i",lev_num);
         init_bbox[(num-1)*2*AMRD_dim]=-1e6;
         AMRD_real_param(pfile,buffer,&init_bbox[(num-1)*2*AMRD_dim],2*AMRD_dim);
         if (init_bbox[(num-1)*2*AMRD_dim]==-1e6) num--; // _b is optional
      }
   }
   else
   {
      if (AMRD_echo_params && my_rank==0) printf("\nAMR off: initializing base level\n\n");
      lev_num=num=1;
   }
   
   //Turn off unused variables
   for (i=0; i<num_unused_vars; i++) PAMR_disable_var(unused_vars[i]);
   
   if (!(PAMR_compose_hierarchy(1,lev_num,num,init_lev,init_bbox,AMRD_t0))) AMRD_stop("PAMR_compose_hierarchy failed ...",0);

   for (i=1; i<=lev_num; i++) call_app(app_AMRH_var_clear,i,PAMR_AMRH);

   if (AMRD_do_ex) PAMR_excision_on(chr_v,chr_c,app_fill_ex_mask,AMRD_ex,1);

   if (AMRD_echo_params && my_rank==0) printf("\ninit_context() complete\n");
   if (AMRD_echo_params && my_rank==0) printf("\n%s\n",line_break);

   if (AMRD_TRE_vars)            free(AMRD_TRE_vars); 
   if (mg_interp2_v)             free(mg_interp2_v); 
   if (mg_hw_restr_v)            free(mg_hw_restr_v); 
   if (mg_fw_restr_v)            free(mg_fw_restr_v);
   if (mg_hw_restr_bdy_si_v)     free(mg_hw_restr_bdy_si_v);
   if (mg_sync_v)                free(mg_sync_v);
   if (amr_inject_v)             free(amr_inject_v);
   if (amr_interp1_v)            free(amr_interp1_v);
   if (amr_interp2_v)            free(amr_interp2_v);
   if (amr_interp4_v)            free(amr_interp4_v);
   if (amr_sync_v)               free(amr_sync_v);
   if (amr_transfer2_v)          free(amr_transfer2_v);
   if (amr_transfer1_v)          free(amr_transfer1_v);
   if (amr_c_to_v2_v)            free(amr_c_to_v2_v);
   if (amr_v_to_c2_v)            free(amr_v_to_c2_v);
   if (amr_bdy_interp1_v)        free(amr_bdy_interp1_v);
   if (amr_bdy_interp2_v)        free(amr_bdy_interp2_v);
   if (amr_transfer4_v)          free(amr_transfer4_v);
   if (even_vars_x0min)          free(even_vars_x0min);
   if (even_vars_x1min)          free(even_vars_x1min);
   if (even_vars_x2min)          free(even_vars_x2min);
   if (even_vars_x0max)          free(even_vars_x0max);
   if (even_vars_x1max)          free(even_vars_x1max);
   if (even_vars_x2max)          free(even_vars_x2max);
   if (odd_vars_x0min)           free(odd_vars_x0min);
   if (odd_vars_x1min)           free(odd_vars_x1min);
   if (odd_vars_x2min)           free(odd_vars_x2min);
   if (odd_vars_x0max)           free(odd_vars_x0max);
   if (odd_vars_x1max)           free(odd_vars_x1max);
   if (odd_vars_x2max)           free(odd_vars_x2max);
   if (unused_vars)              free(unused_vars);

   return;
}

//=============================================================================
// save0 saves *all* save_n_vars[] 
//=============================================================================
void save0(int iter)
{
   int Lf,L,tl,i,ltrace=0,itag,do_close=0;
   char ltag[100];

   IFL printf("save0 ...\n");

   if (AMRD_save_iter_mod>0) 
   {
      itag=((int)(iter/AMRD_save_iter_mod))*AMRD_save_iter_mod;
      sprintf(ltag,"_%i",itag);
      if (iter % AMRD_save_iter_mod == (AMRD_save_iter_mod-1)) do_close=1;
   }
   else ltag[0]=0;

   Lf=PAMR_get_max_lev(PAMR_AMRH);
   for (L=1; L<=Lf; L++) call_app(app_pre_io_calc,L,PAMR_AMRH);

   for (tl=1; tl<=AMRD_num_evo_tl; tl++)
   {
      i=0; 
      while(i<AMRD_num_save_n_vars[tl-1]) 
      {
         if (do_close) PAMR_save_gfn_close(AMRD_save_n_vars[tl-1][i++],PAMR_AMRH,tl,-1,-1.0,AMRD_save_tag,ltag);
         else PAMR_save_gfn(AMRD_save_n_vars[tl-1][i++],PAMR_AMRH,tl,-1,-1.0,AMRD_save_tag,ltag);
      }
   }

   IFL printf("save0 ... done\n");
}

//=============================================================================
// saveL saves level lev save_n_vars[], and all finer levels
// if (AMRD_save_all_finer)
//=============================================================================
void saveL(int lev,int iter)
{
   int Lf,L,tl,i,ltrace=0,itag,do_close=0;
   char ltag[100];

   IFL printf("saveL ...\n");

   if (AMRD_save_iter_mod>0) 
   {
      itag=((int)(iter/AMRD_save_iter_mod))*AMRD_save_iter_mod;
      sprintf(ltag,"_L%i_%i",lev,itag);
      if (iter % AMRD_save_iter_mod == (AMRD_save_iter_mod-1)) do_close=1;
   }
   else sprintf(ltag,"_L%i",lev);

   if (AMRD_save_all_finer) Lf=PAMR_get_max_lev(PAMR_AMRH); else Lf=lev;
   for (L=lev; L<=Lf; L++) call_app(app_pre_io_calc,L,PAMR_AMRH);

   for (L=lev; L<=Lf; L++)
   {
      for (tl=1; tl<=AMRD_num_evo_tl; tl++)
      {
         i=0; 
         while(i<AMRD_num_save_n_vars[tl-1]) 
         {
            if (do_close) PAMR_save_gfn_close(AMRD_save_n_vars[tl-1][i++],PAMR_AMRH,tl,L,-1.0,AMRD_save_tag,ltag);
            else PAMR_save_gfn(AMRD_save_n_vars[tl-1][i++],PAMR_AMRH,tl,L,-1.0,AMRD_save_tag,ltag);
         }
      }
   }

   IFL printf("saveL ... done\n");
}

//=================================================================================
// if (dir==AMRD_CP_SAVE)
//   copy n bytes from p to *q,
//   increments *q, and adds number of bytes to tot_size
// else if (dir==PAMR_CP_DIR_RESTORE)
//   same, but direction of copy is reversed. 
//=================================================================================
#define AMRD_CP_MAX_SIZE 1400000
void amrd_copy_block(char *p, char **q, int n, int dir, int *tot_size)
{
   char *p0;
   int n0;

   if (n==0) return;
   
   if ((*tot_size+n) > (AMRD_CP_MAX_SIZE)) 
      AMRD_stop("amrd_copy_block: error ... AMRD_CP_MAX_SIZE too small\n","");
   *tot_size+=n;
   
   n0=n;

   p0=p;

   if (dir==AMRD_CP_SAVE) while(n0--) *(*q)++=*p0++;
   else while(n0--) *p0++=*(*q)++;

   return;
}

//=============================================================================
// if dir == AMRD_CP_SAVE
// then saves all needed global and user parameters to cp_file_AMRD_my_rank.sdf 
//
// if dir == AMRD_CP_RESTORE
// then either restores from that file or cp_file_AMRD_0.sdf
//=============================================================================
int amrd_do_cp(char *cp_file, int dir)
{
   int size0,size,ret=0,i,id_tag=-1,ver=AMRD_CP_VERSION,j;
   int t_AMRD_MAX_VARS,t_AMRD_MAX_TIMES,t_PAMR_MAX_GFNS,t_PAMR_MAX_DIM,t_PAMR_MAX_LEVS;
   char *data,*data0=0,file_name[256];

   if (dir==AMRD_CP_RESTORE)
   {
      sprintf(file_name,"%s_AMRD_%i",cp_file,my_rank);
      if (!(gft_read_shape(file_name,1,&size0)))
      {
         sprintf(file_name,"%s_AMRD_0",cp_file);
         if (!(gft_read_shape(file_name,1,&size0))) { printf("amrd_do_cp: error opening file %s\n",file_name); goto fin;}
      }
      size0*=sizeof(double);
      if (!(data0=data=malloc(size0))) { printf("amrd_do_cp: error ... out of memory\n"); goto fin;}
      if (!(gft_read_brief(file_name,1,(double *)data))) { printf("amrd_do_cp: error reading file %s\n",file_name); goto fin;}
   }
   else
   {
      size0=AMRD_CP_MAX_SIZE+AMRD_cp_data_size;
      if (!(data0=malloc(size0))) {printf("amrd_do_cp: error ... out of memory\n"); return 0;}
   }

   data=data0; size=0;

   // in old-cp stuff, no version info nor variable size indicators were saved
   // first variable was lsteps (i.e. >0), so use a -1 tag to indicate old vs new

   if (dir==AMRD_CP_SAVE || (dir==AMRD_CP_RESTORE && ((int *)data)[0]==id_tag))
   {
      t_AMRD_MAX_VARS=AMRD_MAX_VARS;
      t_AMRD_MAX_TIMES=AMRD_MAX_TIMES;
      t_PAMR_MAX_GFNS=PAMR_MAX_GFNS;
      t_PAMR_MAX_LEVS=PAMR_MAX_LEVS;
      t_PAMR_MAX_DIM=PAMR_MAX_DIM;
      amrd_copy_block((char *)&id_tag,&data,sizeof(int),dir,&size);
      amrd_copy_block((char *)&ver,&data,sizeof(int),dir,&size);
      amrd_copy_block((char *)&t_AMRD_MAX_VARS,&data,sizeof(int),dir,&size);
      amrd_copy_block((char *)&t_AMRD_MAX_TIMES,&data,sizeof(int),dir,&size);
      amrd_copy_block((char *)&t_PAMR_MAX_GFNS,&data,sizeof(int),dir,&size);
      amrd_copy_block((char *)&t_PAMR_MAX_LEVS,&data,sizeof(int),dir,&size);
      amrd_copy_block((char *)&t_PAMR_MAX_DIM,&data,sizeof(int),dir,&size);
      if (ver!=AMRD_CP_VERSION && !my_rank) printf("WARNING: AMRD_CP_VERSION mismatch\n");

      // cannot yet handle the situation where the size of arrays have decreased.
      // (or when MAX_DIM changes)
     
      if (dir==AMRD_CP_RESTORE)
      {
         if (t_AMRD_MAX_VARS>AMRD_MAX_VARS)
         {
            if (!my_rank) printf("ERROR: AMRD_MAX_VARS in cp file (%i) is greater than current value (%i)\n",
                                 t_AMRD_MAX_VARS,AMRD_MAX_VARS);
            AMRD_stop("","");
         }
         if (t_AMRD_MAX_TIMES>AMRD_MAX_TIMES)
         {
            if (!my_rank) printf("ERROR: AMRD_MAX_TIMES in cp file (%i) is greater than current value (%i)\n",
                                 t_AMRD_MAX_TIMES,AMRD_MAX_TIMES);
            AMRD_stop("","");
         }
         if (t_PAMR_MAX_GFNS>PAMR_MAX_GFNS)
         {
            if (!my_rank) printf("ERROR: PAMR_MAX_GFNS in cp file (%i) is greater than current value (%i)\n",
                                 t_PAMR_MAX_GFNS,PAMR_MAX_GFNS);
            AMRD_stop("","");
         }
         if (t_PAMR_MAX_LEVS>PAMR_MAX_LEVS)
         {
            if (!my_rank) printf("ERROR: PAMR_MAX_LEVS in cp file (%i) is greater than current value (%i)\n",
                                 t_PAMR_MAX_LEVS,PAMR_MAX_LEVS);
            AMRD_stop("","");
         }
         if (t_PAMR_MAX_DIM!=PAMR_MAX_DIM)
         {
            if (!my_rank) printf("ERROR: PAMR_MAX_DIM in cp file (%i) is different from current value (%i)\n",
                                 t_PAMR_MAX_DIM,PAMR_MAX_DIM);
            AMRD_stop("","");
         }
      }
   }
   else
   {
      AMRD_stop("ERROR: amrd_do_cp ... invalid id_tag or 'dir'","");
   }

   // now just copy all needed variables to/from data
   // many variables are read from the parameter file, and in principle
   // we can recalculate the gfn's from this. However, for simplicity
   // save them here. At some point it may be good to save *all* variables
   // that should not be changed upon a restore here, but for now keep it simple.
   
   amrd_copy_block((char *)AMRD_lsteps,&data,sizeof(int)*t_PAMR_MAX_LEVS,dir,&size);
   amrd_copy_block((char *)AMRD_MG_res_gfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   amrd_copy_block((char *)AMRD_MG_lop_gfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   amrd_copy_block((char *)AMRD_MG_brs_gfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   amrd_copy_block((char *)AMRD_MG_rhs_gfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   amrd_copy_block((char *)AMRD_MG_f_gfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   if (ver>1) amrd_copy_block((char *)AMRD_fcs_gfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   amrd_copy_block((char *)AMRD_tn_diss_pbt,&data,sizeof(int)*t_AMRD_MAX_VARS*2*t_PAMR_MAX_DIM,dir,&size);
   amrd_copy_block((char *)AMRD_tnp1_diss_pbt,&data,sizeof(int)*t_AMRD_MAX_VARS*2*t_PAMR_MAX_DIM,dir,&size);
   amrd_copy_block((char *)AMRD_rg_diss_pbt,&data,sizeof(int)*t_AMRD_MAX_VARS*2*t_PAMR_MAX_DIM,dir,&size);
   amrd_copy_block((char *)&AMRD_cmask_gfn,&data,sizeof(int)*1,dir,&size);
   amrd_copy_block((char *)&AMRD_cmask_mg_gfn,&data,sizeof(int)*1,dir,&size);
   amrd_copy_block((char *)&AMRD_cmask_c_gfn,&data,sizeof(int)*1,dir,&size);
   amrd_copy_block((char *)&AMRD_fc_mask_gfn,&data,sizeof(int)*1,dir,&size);
   amrd_copy_block((char *)&AMRD_tre_gfn,&data,sizeof(int)*1,dir,&size);
   amrd_copy_block((char *)&AMRD_mg_tre_gfn,&data,sizeof(int)*1,dir,&size);
   amrd_copy_block((char *)&AMRD_chr_gfn,&data,sizeof(int)*1,dir,&size);
   amrd_copy_block((char *)&AMRD_chr_mg_gfn,&data,sizeof(int)*1,dir,&size);
   if (ver>1) amrd_copy_block((char *)&AMRD_chr_c_gfn,&data,sizeof(int)*1,dir,&size);
   if (ver>1) amrd_copy_block((char *)&AMRD_chr_mg_c_gfn,&data,sizeof(int)*1,dir,&size);
   if (ver>2) amrd_copy_block((char *)&AMRD_wavg_c_gfn,&data,sizeof(int)*1,dir,&size);
   if (ver>2) amrd_copy_block((char *)&AMRD_wavg_mg_c_gfn,&data,sizeof(int)*1,dir,&size);
   // need to do the following two like this, due to the manner in which c-stores arrays,
   // and AMRD_MAX_VARS might have changed
   for (j=0; j<t_AMRD_MAX_TIMES; j++) 
      amrd_copy_block((char *)AMRD_AMR_f_gfn[j],&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   for (j=0; j<t_AMRD_MAX_TIMES; j++) 
      amrd_copy_block((char *)AMRD_AMR_mgf_gfn[j],&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   if (ver>1)
   {
      for (j=0; j<t_AMRD_MAX_TIMES; j++) 
        amrd_copy_block((char *)AMRD_past_bdy_interp_vars_gfn[j],&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
      amrd_copy_block((char *)AMRD_past_bdy_interp_op,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   }
   amrd_copy_block((char *)AMRD_AMR_mgf_extrap_tm1_gfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   amrd_copy_block((char *)AMRD_AMR_mgf_extrap_tm2_gfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   amrd_copy_block((char *)AMRD_tnp1_diss_gfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   amrd_copy_block((char *)AMRD_tn_diss_gfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   amrd_copy_block((char *)AMRD_ex_repop1_sgfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   amrd_copy_block((char *)AMRD_ex_repop2_sgfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   amrd_copy_block((char *)AMRD_ex_repop3_sgfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   amrd_copy_block((char *)AMRD_ex_repop4_sgfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   amrd_copy_block((char *)AMRD_tnp1_liipb_gfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   amrd_copy_block((char *)AMRD_tnp1_liiab_gfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   amrd_copy_block((char *)AMRD_tnp1_liibb_gfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   amrd_copy_block((char *)AMRD_global_var_norms,&data,sizeof(real)*t_PAMR_MAX_GFNS,dir,&size);
   amrd_copy_block((char *)AMRD_rg_diss_gfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   amrd_copy_block((char *)AMRD_interp_AMR_bdy_f_gfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   amrd_copy_block((char *)AMRD_f_tre_gfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   amrd_copy_block((char *)AMRD_f_tre_fn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   amrd_copy_block((char *)&AMRD_w1_gfn,&data,sizeof(int)*1,dir,&size);
   if (ver>1) amrd_copy_block((char *)&AMRD_w1_c_gfn,&data,sizeof(int)*1,dir,&size);
   amrd_copy_block((char *)AMRD_tre_valid,&data,sizeof(int)*t_PAMR_MAX_LEVS,dir,&size);
   amrd_copy_block((char *)AMRD_MG_cnst_data_gfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   amrd_copy_block((char *)&AMRD_state,&data,sizeof(int)*1,dir,&size);
   amrd_copy_block((char *)&AMRD_is_t0,&data,sizeof(int)*1,dir,&size);
   amrd_copy_block((char *)&AMRD_top,&data,sizeof(int)*1,dir,&size);
   amrd_copy_block((char *)AMRD_stack,&data,sizeof(int)*2*t_PAMR_MAX_LEVS,dir,&size);
   amrd_copy_block((char *)&rgs_t,&data,sizeof(real)*1,dir,&size);
   amrd_copy_block((char *)&rgs_L1,&data,sizeof(int)*1,dir,&size);
   amrd_copy_block((char *)&rgs_Lf,&data,sizeof(int)*1,dir,&size);
   amrd_copy_block((char *)&AMRD_num_tn_diss_vars,&data,sizeof(int)*1,dir,&size);
   amrd_copy_block((char *)&AMRD_num_f_tre_vars,&data,sizeof(int)*1,dir,&size);
   amrd_copy_block((char *)&AMRD_tnum_rg_diss_vars,&data,sizeof(int)*1,dir,&size);
   amrd_copy_block((char *)&AMRD_curr_cp_tag,&data,1,dir,&size);
   if (ver>1) amrd_copy_block((char *)&AMRD_using_cc_tre,&data,sizeof(int)*1,dir,&size);
   if (ver>1) amrd_copy_block((char *)AMRD_AMRH_work_inject_gfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   if (ver>1) amrd_copy_block((char *)AMRD_AMRH_work_inject_op,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   if (ver>2) amrd_copy_block((char *)AMRD_work_repop_gfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   if (ver>2) amrd_copy_block((char *)AMRD_inject_wavg_gfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   if (ver>2) amrd_copy_block((char *)AMRD_MG_inject_wavg_gfn,&data,sizeof(int)*t_AMRD_MAX_VARS,dir,&size);
   if (ver>2) amrd_copy_block((char *)&AMRD_num_MG_inject_wavg_vars,&data,sizeof(int)*1,dir,&size);

   //align to double
   data+=(sizeof(double)-size%(sizeof(double)));
   size+=(sizeof(double)-size%(sizeof(double)));

   app_user_cp(dir,data);

   size+=AMRD_cp_data_size;
   data+=AMRD_cp_data_size;

   size+=(sizeof(double)-size%(sizeof(double )));
   data+=(sizeof(double)-size%(sizeof(double )));

   if (dir==AMRD_CP_SAVE)
   {
      sprintf(file_name,"%s_AMRD_%i",cp_file,my_rank);
      size/=sizeof(double);
      if (!(gft_out_brief(file_name,0.0,&size,1,(double *)data0))) { printf("amrd_do_cp: error saving file %s\n",file_name); goto fin;}
      // this is important for the integrity of a saved file if the job terminates abnormally
      if (dir==AMRD_CP_SAVE) gft_close(file_name);
   }

   ret=1;

fin:
   if (data0) free(data0);
   return ret; 
}

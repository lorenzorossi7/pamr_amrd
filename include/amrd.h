#ifndef _AMRD_H
#define _AMRD_H

#include "pamr.h"

/*============================================================================= */
/* amrd.h - adaptive mesh refinement driver                                     */
/*                                                                              */
/*  Copyright 2003/4 F.Pretorius                                                */
/*                                                                              */
/*============================================================================= */

/*============================================================================= */
/* ... API for application                                                      */
/*============================================================================= */

void amrd(int argc, char **argv, 
          int (*app_id_f)(void),
          void (*app_var_pre_init_f)(char *pfile),
          void (*app_var_post_init_f)(char *pfile),
          void (*app_AMRH_var_clear_f)(void),
          void (*app_free_data_f)(void),
          void (*app_t0_cnst_data_f)(void),
          real (*app_evo_residual_f)(void),
          real (*app_MG_residual_f)(void),
          void (*app_evolve_f)(int iter),
          real (*app_MG_relax_f)(void),
          void (*app_L_op_f)(void),
          void (*app_pre_io_calc_f)(void),
          void (*app_scale_tre_f)(void),
          void (*app_post_regrid_f)(void),
          void (*app_post_tstep_f)(int L),
          void (*app_fill_ex_mask_f)(real *mask, int dim, int *shape, real *bbox, real excised),
          void (*app_fill_bh_bboxes_f)(real *bbox, int *num, int max_num));

/*============================================================================= */
/* some utility routines                                                        */
/*============================================================================= */

real *AMRD_get_global_norms(void);
void AMRD_ev_ldptr(void);

/*============================================================================= */
/* define a post-tstep hook                                                     */
/*============================================================================= */
extern void (*app_pre_tstep)(int L);
void amrd_set_app_pre_tstep_hook(void (*app_pre_tstep_f)(int L));

/*============================================================================= */
/* define a hook to initialize elliptic variables that are only define at the   */
/* initial time                                                                 */
/*============================================================================= */
extern void (*app_elliptic_vars_t0_init)(void);
void amrd_set_elliptic_vars_t0_init(void (*app_elliptic_vars_t0_init)(void));

/*============================================================================= */
/* defines a check-pointing hook                                                */
/*============================================================================= */
#define AMRD_CP_SAVE 0
#define AMRD_CP_RESTORE 1
extern void (*app_user_cp)(int save_restore, char *data);
void amrd_set_app_user_cp_hook(void (*app_pre_tstep_f)(int save_restore, char *data), int cp_data_size);
extern int AMRD_cp_data_size,AMRD_cp_restart;
extern real AMRD_cp_delta_t_hrs,AMRD_cp_first_delta_t_hrs;
extern char *AMRD_cp_save_fname;
extern char *AMRD_cp_restore_fname;
extern char AMRD_curr_cp_tag;

/*============================================================================= */
/* error handler                                                                */
/*============================================================================= */
void AMRD_stop(char *err1, char *err2);

/*============================================================================= */
/* parameter input                                                              */
/*============================================================================= */
void AMRD_int_param(char *pfile, char *name, int *var, int size);
void AMRD_real_param(char *pfile, char *name, real *var, int size);
void AMRD_str_param(char *pfile, char *name, char **var, int size);
void AMRD_ivec_param(char *pfile, char *name, int *var, int size);

/*============================================================================= */
/* Some of the following set of variables should more appropriately be          */
/* considered variables internal to amrd. However, since I'm defining them as   */
/* external globals, I'm including the definitions here to highlight            */
/* potential conflicts with application specific global variable names          */
/*============================================================================= */

/*============================================================================= */
/* Prototypes for required functions that the user must specify. Grids are      */
/* looped over using the PAMR iterator mechanism, hence grid-based hooks        */
/* can read their arguments via PAMR_get_g_...()                                */
/*============================================================================= */

/*============================================================================= */
/* Returns 0 to use default mechanism, or is expected to calculate              */
/* the correct initial hierarhcy and return 1:                                  */
/*============================================================================= */
extern int (*app_id)(void);

/*============================================================================= */
/* Sets custom parameters, variables, etc. Split up into two segments,          */
/* one called before the pamr context is initialized and standard               */
/* parameters are read, and the other afterwards.                               */
/* Argument is parameter file                                                   */
/*============================================================================= */
extern void (*app_var_pre_init)(char *pfile);
extern void (*app_var_post_init)(char *pfile);

/*============================================================================= */
/* Sets all AMR/MG variables to their 'zero' values:                            */
/*============================================================================= */
extern void (*app_AMRH_var_clear)(void);

/*============================================================================= */
/* Initial data for free fields:                                                */
/*============================================================================= */
extern void (*app_free_data)(void);

/*============================================================================= */
/* Initial constraint data --- called after each MG iteration:                  */
/*============================================================================= */
extern void (*app_t0_cnst_data)(void);

/*============================================================================= */
/* Called prior to saving variables to disk                                     */
/*============================================================================= */
extern void (*app_pre_io_calc)(void);

/*============================================================================= */
/* Returns some norm of the residual for the evolution variables,               */
/* called after an evolution iteration:                                         */
/*============================================================================= */
extern real (*app_evo_residual)(void);

/*============================================================================= */
/* Returns some norm of the residual for the MG variables, *AND*                */
/* stores the point-wise residual in "f_res" for each MG variable "f" (for      */
/* computing new RHS's)                                                         */
/*============================================================================= */
extern real (*app_MG_residual)(void);

/*============================================================================= */
/* Performs 1 iteration of the evolution equations                              */
/*============================================================================= */
extern void (*app_evolve)(int iter);

/*============================================================================= */
/* Performs 1 relaxation sweep of the MG equations, and returns an estimate     */
/* of the norm of the residual.                                                 */
/*============================================================================= */
extern real (*app_MG_relax)(void);

/*============================================================================= */
/* Computes the differential operator L acting on the current grid,             */
/* in a region specified by the grid function "cmask". Stores the result        */
/* in "f_lop" for each MG variable "f"                                          */
/*============================================================================= */
extern void (*app_L_op)(void);

/*============================================================================= */
/* Called after calculating the TRE for all variables                           */
/*============================================================================= */
extern void (*app_scale_tre)(void);

/*============================================================================= */
/* Called after regridding (to re-initialize time-independent grid functions)   */
/*============================================================================= */
extern void (*app_post_regrid)(void);

/*============================================================================= */
/* Called after tstep (no grid-iteration)                                       */
/*============================================================================= */
extern void (*app_post_tstep)(int L);

/*============================================================================= */
/* passed to PAMR's PAMR_excision_on() function                                 */
/*============================================================================= */
extern void (*app_fill_ex_mask)(real *mask, int dim, int *shape, real *bbox, real excised);

/*============================================================================= */
/* used if TRE_sgpbh is on, returns a list of bounding boxes surrounding the    */
/* holes (up to a maximum of max_num)                                           */
/*============================================================================= */
extern void (*app_fill_bh_bboxes)(real *bbox, int *num, int max_num);

/*----------------------------------------------------------------------------- */

/*============================================================================= */
/* Some run-time variables read from parameter file --- others, in particular   */
/* those managed by libpamr, are not stored as globals; the current state       */
/* should be queried/set via the pamr functions.                                */
/*============================================================================= */

/*============================================================================= */
/* Variable name lists:                                                         */
/* The user needs to declare at least the hyperbolic and elliptic               */
/* variables --- these will then be defined automatically.                      */
/* Work variables can also be declared via these lists.                         */
/*                                                                              */
/* For each elliptic variable "f", the following additional MG                  */
/* variables will be created (note that the rhs's are MG generated only ...     */
/* users must place there rhs's on the lhs):                                    */
/*                                                                              */
/* f_res                                                                        */
/* f_lop                                                                        */
/* f_rhs                                                                        */
/* f_brs                                                                        */
/*                                                                              */
/* and the following AMR work variables (currently only second order extrapolation)*/
/*                                                                              */
/* f_extrap_tm1                                                                 */
/* f_extrap_tm2                                                                 */
/*                                                                              */
/* For each TRE_vars variable (in parameter file)                               */
/*                                                                              */
/* f_tre                                                                        */
/*                                                                              */
/* other variables:                                                             */
/*                                                                              */
/* cmask : child-mask. In vcycle, used to set region where L_op should act.     */
/* tre : net, flagged tre                                                       */
/*                                                                              */
/* AMRD_w1, AMRD_w2, ... : work arrays                                          */
/*============================================================================= */

#define AMRD_MAX_CP_TIMES 256
#define AMRD_CP_VERSION 1
#define AMRD_MAX_IVEC_SIZE 50000
#define AMRD_MAX_VARS 256
#define AMRD_MAX_TIMES 5
#define AMRD_MAX_VNAME_SIZE 32

extern int AMRD_num_hyperbolic_vars;
extern int AMRD_num_elliptic_vars_t0;
extern int AMRD_num_elliptic_vars;
extern int AMRD_num_AMRH_work_vars;
extern int AMRD_num_MGH_work_vars;
extern int AMRD_num_AMRH_work_in_MGH_vars;

extern char *AMRD_hyperbolic_vars[AMRD_MAX_VARS];
extern char *AMRD_elliptic_vars_t0[AMRD_MAX_VARS];
extern char *AMRD_elliptic_vars[AMRD_MAX_VARS];
extern char *AMRD_AMRH_work_vars[AMRD_MAX_VARS];
extern char *AMRD_MGH_work_vars[AMRD_MAX_VARS];
extern char *AMRD_AMRH_work_in_MGH_vars[AMRD_MAX_VARS];

extern char *app_name;

/*============================================================================= */
/* basic options                                                                */
/*============================================================================= */

extern int AMRD_dim;     
extern int AMRD_num_evo_tl; 
extern int AMRD_ic_n; /* which time level to initialize at t=0 */

extern int AMRD_base_shape[PAMR_MAX_DIM];
extern real AMRD_base_bbox[PAMR_MAX_DIM*2];
extern int AMRD_max_lev;
extern real AMRD_t0;

/*============================================================================= */
/* iterator and MG controls                                                     */
/*============================================================================= */
extern int AMRD_steps;     
extern int AMRD_lsteps[PAMR_MAX_LEVS];     
extern int AMRD_evo_max_iter;
extern int AMRD_evo_min_iter;
extern int AMRD_MG_max_iter;
extern int AMRD_MG_min_iter;
extern int AMRD_MG_max_citer;
extern real AMRD_evo_tol;  
extern real AMRD_MG_tol;
extern real AMRD_MG_crtol;
extern real AMRD_MG_w0_r,AMRD_MG_w0_i;
extern int AMRD_evo_ssc;
extern int AMRD_MG_pre_swp;
extern int AMRD_MG_pst_swp;
#define AMRD_NP1_IG_NONE 0
#define AMRD_NP1_IG_FIRST_ORDER 1
extern int AMRD_np1_initial_guess;

/*============================================================================= */
/* MG extrapolation options                                                     */
/*============================================================================= */
#define AMRD_MG_EXTRAP_2ND_FROM_ETM1 0
#define AMRD_MG_EXTRAP_2ND_FROM_TN 1
extern int AMRD_MG_extrap_method;
extern real AMRD_MG_eps_c;

/*============================================================================= */
/* id options                                                                   */
/*============================================================================= */
extern int AMRD_id_method;
extern int AMRD_id_pl_method,AMRD_id_pl_steps,AMRD_id_user_mg_pl;
extern real AMRD_id_pl_lambda;

/*============================================================================= */
/* TRE parameters                                                               */
/*============================================================================= */
extern real AMRD_TRE_max;
extern int AMRD_TRE_norm;
extern int AMRD_TRE_buffer;
extern int AMRD_TRE_ibc_buffer;
extern int AMRD_TRE_exc_buffer;
extern int AMRD_TRE_exc_buffer_lmin;
extern int AMRD_TRE_sgpbh;
extern int AMRD_TRE_ibcp_buffer;
extern int AMRD_regrid_interval;
extern int AMRD_regrid_min_lev;
extern int AMRD_skip_frg;

#define AMRD_NO_REGRID_SCRIPT 0
#define AMRD_READ_REGRID_SCRIPT 1
#define AMRD_WRITE_REGRID_SCRIPT 2
extern int AMRD_regrid_script;
extern char *AMRD_regrid_script_name;

/*============================================================================= */
/* I/O parameters                                                               */
/*============================================================================= */
extern char *AMRD_restart_file;
extern char *AMRD_cp_tag;
extern real AMRD_cp_times[AMRD_MAX_CP_TIMES];
extern char *AMRD_save_tag;

extern int *AMRD_save_ivec[PAMR_MAX_LEVS];
extern int *AMRD_save_ivec0;
extern char *AMRD_vt_times;
extern char *AMRD_save_n_vars[AMRD_MAX_TIMES][AMRD_MAX_VARS];
extern char *AMRD_save_mg_vars[AMRD_MAX_VARS];

extern int AMRD_save_iter_mod;

extern int AMRD_save_all_finer;

extern char AMRD_coord_names[PAMR_MAX_DIM][AMRD_MAX_VNAME_SIZE];

/*============================================================================= */
/* tracing parameters                                                           */
/*============================================================================= */
extern int AMRD_echo_params;
extern int AMRD_MG_trace;
extern int AMRD_MG_DV_trace;
extern real AMRD_MG_DV_trace_t_on;
extern real AMRD_MG_DV_trace_t_off;
extern int AMRD_evo_trace;
extern int AMRD_evo_DV_trace;
extern real AMRD_evo_DV_trace_t_on;
extern real AMRD_evo_DV_trace_t_off;
extern int AMRD_ID_DV_trace;

/*============================================================================= */
/* to allow the user to calculate normalized residuals                          */
/* (currently, global norms are updated every coarse grid time step only)       */
/*============================================================================= */
#define AMRD_GLOBAL_VAR_INF_NORM 1
#define AMRD_GLOBAL_VAR_L2_NORM 2
extern int AMRD_calc_global_var_norms;
extern int AMRD_global_var_norm_type;
extern real AMRD_global_var_norms[PAMR_MAX_GFNS];
extern real AMRD_global_var_norm_floor;

/*============================================================================= */
/* excision parameters                                                          */
/*============================================================================= */
extern real AMRD_ex;
extern int AMRD_do_ex;

/*============================================================================= */
/* a utility function that will call dmrepop3d1_ for all AMRH                   */
/* variables at all time levels, together with apropriate synchronization.      */
/* --- repeats the process n-times for repopulation of up to n zones            */
/* on the finest level                                                          */
/*============================================================================= */
void AMRD_repopulate(int n, int default_order);

extern int AMRD_num_ex_repop1_vars;
extern char *AMRD_ex_repop1_vars[AMRD_MAX_VARS];

extern int AMRD_num_ex_repop2_vars;
extern char *AMRD_ex_repop2_vars[AMRD_MAX_VARS];

extern int AMRD_num_ex_repop3_vars;
extern char *AMRD_ex_repop3_vars[AMRD_MAX_VARS];

extern int AMRD_num_ex_repop4_vars;
extern char *AMRD_ex_repop4_vars[AMRD_MAX_VARS];

/*============================================================================= */
/* clustering parameters                                                        */
/*============================================================================= */
extern int AMRD_cls_method;
extern int AMRD_cls_merge_dist;
extern int AMRD_cls_align_mode;

/*============================================================================= */
/* automatic Kreiss-Oliger filtering                                            */
/*============================================================================= */
extern int AMRD_num_tnp1_diss_vars;
extern char *AMRD_tnp1_diss_vars[AMRD_MAX_VARS];
extern int AMRD_tnp1_diss_pbt[AMRD_MAX_VARS][PAMR_MAX_DIM*2];
extern int AMRD_num_tn_diss_vars;
extern char *AMRD_tn_diss_vars[AMRD_MAX_VARS];
extern int AMRD_tn_diss_pbt[AMRD_MAX_VARS][PAMR_MAX_DIM*2];
extern real AMRD_eps_diss,AMRD_tn_eps_diss,AMRD_tnp1_eps_diss;
extern int AMRD_diss_bdy;
extern int AMRD_repop_diss_bdy;
extern int AMRD_diss_freq;
extern int AMRD_diss_all_past;
extern int AMRD_num_rg_diss_vars;
extern int AMRD_tnum_rg_diss_vars;
extern char *AMRD_rg_diss_vars[AMRD_MAX_VARS];
extern real AMRD_rg_eps_diss;

/*============================================================================= */
/* linearize points 1 in from boundaries during evolution step                  */
/* liip/a/bb : linearly interpolate i-1 physical/amr/both boundaries            */
/*============================================================================= */
extern int AMRD_num_tnp1_liipb_vars;
extern int AMRD_num_tnp1_liiab_vars;
extern int AMRD_num_tnp1_liibb_vars;
extern char *AMRD_tnp1_liipb_vars[AMRD_MAX_VARS];
extern char *AMRD_tnp1_liiab_vars[AMRD_MAX_VARS];
extern char *AMRD_tnp1_liibb_vars[AMRD_MAX_VARS];

/*============================================================================= */
/* AMR boundary 'smoothing' for elliptic functions (or others)                  */
/*============================================================================= */
extern int AMRD_num_interp_AMR_bdy_vars;
extern int AMRD_interp_AMR_bdy_order;
extern char *AMRD_interp_AMR_bdy_vars[AMRD_MAX_VARS];

/*============================================================================= */
/* other globals                                                                */
/*============================================================================= */

extern int pamr_context;
extern int my_rank;

/*============================================================================= */
/* cmask (characteristic, or child mask) defines                                */
/*============================================================================= */
#define AMRD_CMASK_ON 1.0
#define AMRD_CMASK_OFF 0.0

/*============================================================================= */
/* "internal" globals                                                           */
/*============================================================================= */
extern int AMRD_MG_res_gfn[AMRD_MAX_VARS],AMRD_MG_lop_gfn[AMRD_MAX_VARS],AMRD_MG_brs_gfn[AMRD_MAX_VARS],
           AMRD_MG_rhs_gfn[AMRD_MAX_VARS],AMRD_MG_f_gfn[AMRD_MAX_VARS];
extern real *AMRD_MG_res[AMRD_MAX_VARS],*AMRD_MG_lop[AMRD_MAX_VARS],*AMRD_MG_brs[AMRD_MAX_VARS],
            *AMRD_MG_rhs[AMRD_MAX_VARS],*AMRD_MG_f[AMRD_MAX_VARS];
extern int AMRD_cmask_gfn,AMRD_cmask_mg_gfn,AMRD_tre_gfn,AMRD_mg_tre_gfn,AMRD_w1_gfn;

extern real *AMRD_cmask,*AMRD_cmask_mg,*AMRD_tre,*AMRD_mg_tre,*AMRD_w1;

extern int AMRD_chr_gfn,AMRD_chr_mg_gfn;
extern real *AMRD_chr,*AMRD_chr_mg;

extern int AMRD_AMR_f_gfn[AMRD_MAX_TIMES][AMRD_MAX_VARS];
extern real *AMRD_AMR_f[AMRD_MAX_TIMES][AMRD_MAX_VARS];

extern int AMRD_AMR_mgf_gfn[AMRD_MAX_TIMES][AMRD_MAX_VARS];
extern real *AMRD_AMR_mgf[AMRD_MAX_TIMES][AMRD_MAX_VARS];

extern int AMRD_AMR_mgf_extrap_tm1_gfn[AMRD_MAX_VARS];
extern real *AMRD_AMR_mgf_extrap_tm1[AMRD_MAX_VARS];

extern int AMRD_AMR_mgf_extrap_tm2_gfn[AMRD_MAX_VARS];
extern real *AMRD_AMR_mgf_extrap_tm2[AMRD_MAX_VARS];

extern int AMRD_num_f_tre_vars;
extern int AMRD_f_tre_gfn[AMRD_MAX_VARS],AMRD_f_tre_fn[AMRD_MAX_VARS];
extern real *AMRD_f_tre[AMRD_MAX_VARS];
extern char *AMRD_f_tre_vars[AMRD_MAX_VARS];

extern real *AMRD_g_x[PAMR_MAX_DIM],AMRD_g_dx[PAMR_MAX_DIM];
extern int AMRD_g_shape[PAMR_MAX_DIM],AMRD_g_ghost_width[2*PAMR_MAX_DIM],
           AMRD_g_rank,AMRD_g_dim,AMRD_g_size,AMRD_g_Nx,AMRD_g_Ny,AMRD_g_Nz;
extern int AMRD_g_phys_bdy[2*PAMR_MAX_DIM];
extern real AMRD_g_bbox[2*PAMR_MAX_DIM];

extern int AMRD_periodic[PAMR_MAX_DIM];

extern int AMRD_tre_valid[PAMR_MAX_LEVS];

extern int AMRD_ex_repop1_sgfn[AMRD_MAX_VARS];
extern int AMRD_ex_repop2_sgfn[AMRD_MAX_VARS];
extern int AMRD_ex_repop3_sgfn[AMRD_MAX_VARS];
extern int AMRD_ex_repop4_sgfn[AMRD_MAX_VARS];

extern int AMRD_tnp1_diss_gfn[AMRD_MAX_VARS];
extern real *AMRD_tnp1_diss[AMRD_MAX_VARS];
extern int AMRD_tn_diss_gfn[AMRD_MAX_VARS];
extern real *AMRD_tn_diss[AMRD_MAX_VARS];

extern int AMRD_tnp1_liipb_gfn[AMRD_MAX_VARS];
extern int AMRD_tnp1_liiab_gfn[AMRD_MAX_VARS];
extern int AMRD_tnp1_liibb_gfn[AMRD_MAX_VARS];
extern real *AMRD_tnp1_liipb[AMRD_MAX_VARS];
extern real *AMRD_tnp1_liiab[AMRD_MAX_VARS];
extern real *AMRD_tnp1_liibb[AMRD_MAX_VARS];

extern int AMRD_rg_diss_gfn[AMRD_MAX_VARS*AMRD_MAX_TIMES];
extern int AMRD_rg_diss_pbt[AMRD_MAX_VARS][PAMR_MAX_DIM*2];
extern real *AMRD_rg_diss[AMRD_MAX_VARS*AMRD_MAX_TIMES];

extern int AMRD_interp_AMR_bdy_f_gfn[AMRD_MAX_VARS];
extern real *AMRD_interp_AMR_bdy_f[AMRD_MAX_VARS];

extern int AMRD_num_MG_cnst_data_vars;
extern int AMRD_MG_cnst_data_gfn[AMRD_MAX_VARS];
extern real *AMRD_MG_cnst_data[AMRD_MAX_VARS];
extern char *AMRD_MG_cnst_data_vars[AMRD_MAX_VARS];

extern int AMRD_elliptic_vars_t0_gfn[AMRD_MAX_VARS];

#define AMRD_STATE_PRE_INIT 1
#define AMRD_STATE_POST_INIT 2
#define AMRD_STATE_GENERATE_ID 3
#define AMRD_STATE_EVOLVE 4

extern int AMRD_is_t0,AMRD_state;

extern int AMRD_max_t_interp_order;

extern int AMRD_re_interp_width;

extern int AMRD_curr_MG_iter;

extern int AMRD_magic_cookie;

#endif /* _AMRD_H */
